#!/bin/bash
# extract_prot_pbc_chunked.sh
#
# Genera metadyn_prot_pbc.trr a piccoli spezzoni temporali, per restare
# sotto il limite di CPU time del login node (qui: 600s / 10 minuti).
# Ogni chunk e' un trjconv indipendente e leggero; un file di stato
# (.chunks_done) tiene traccia di quali chunk sono gia' stati generati,
# cos lo script puo' essere interrotto e rilanciato senza ripetere
# lavoro gia' fatto.
#
# USO:
#   cd 08-metadyn/metadyn_0   (o qualunque altra replica)
#   bash ../extract_prot_pbc_chunked.sh
#
# Se il processo viene ucciso a meta' (CPU time exceeded), rilancia lo
# stesso comando: riprende dal chunk successivo a quelli gia' salvati.

set -e

#BEGIN INPUT
chunk_ps=2000        # dimensione di ciascun chunk, in ps (2 ns). Riduci
                      # ulteriormente (es. 1000) se anche questo supera
                      # i 600s di CPU time sul tuo sistema.
gmx=gmx_plumed
outGrp=1              # gruppo "Protein" per trjconv
#END INPUT

if [ ! -f metadyn.tpr ] || [ ! -f metadyn.trr ]; then
  echo "-manca metadyn.tpr o metadyn.trr in questa cartella, esco"
  exit 1
fi

# Determina il tempo totale disponibile nella traiettoria grezza,
# interrogando l'ultimo frame con gmx check (operazione leggera, legge
# solo gli header dei frame, non i dati - non dovrebbe avvicinarsi al
# limite di CPU time)
total_ps=$($gmx check -f metadyn.trr 2>&1 | grep -i "Last frame" | awk '{print $NF}')
if [ -z "${total_ps}" ]; then
  echo "-impossibile determinare la durata totale da gmx check, controlla manualmente"
  exit 1
fi
echo "Traiettoria grezza: fino a ${total_ps} ps"

mkdir -p chunks
touch .chunks_done

b=0
while (( $(echo "$b < $total_ps" | bc -l) )); do
  e=$(echo "$b + $chunk_ps" | bc)
  # L'ultimo chunk puo' essere piu' corto: non superare total_ps
  if (( $(echo "$e > $total_ps" | bc -l) )); then
    e=$total_ps
  fi

  chunk_label="${b}-${e}"
  chunk_file="chunks/prot_pbc_${chunk_label}.trr"

  if grep -qx "${chunk_label}" .chunks_done 2>/dev/null; then
    echo "Chunk ${chunk_label} ps: gia' fatto, salto"
    b=$e
    continue
  fi

  echo "Chunk ${chunk_label} ps: estraggo..."
  $gmx trjconv -s metadyn.tpr -f metadyn.trr -o "${chunk_file}" \
    -pbc mol -b "${b}" -e "${e}" << STOP >& "chunks/log_${chunk_label}.out"
${outGrp}
STOP

  if [ -s "${chunk_file}" ]; then
    echo "${chunk_label}" >> .chunks_done
    echo "  -> completato (${chunk_file})"
  else
    echo "  -> ATTENZIONE: chunk ${chunk_label} sembra vuoto/fallito, non segnato come fatto"
    echo "     Controlla chunks/log_${chunk_label}.out prima di rilanciare"
    exit 1
  fi

  b=$e
done

echo ""
echo "Tutti i chunk generati. Verifica quanti sono attesi vs quanti fatti:"
n_done=$(wc -l < .chunks_done)
echo "Chunk completati: ${n_done}"

echo ""
echo "Per unire i chunk in un unico metadyn_prot_pbc.trr, quando pronto:"
echo "  bash ../concat_chunks.sh"
