#!/bin/bash
#SBATCH --job-name 06ModeProj_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

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

# BEGIN INPUT
chunk_ps=2000        # Dimensione chunk in ps (2 ns)
max_ps=80000         # Limite massimo uniformato a 80 ns (80000 ps)
gmx=gmx_plumed
outGrp=1             # Gruppo "Protein"
# END INPUT

if ! command -v $gmx &> /dev/null; then
  echo "- ERRORE: $gmx non è presente nel PATH. Carica i moduli necessari prima di eseguire."
  exit 1
fi

if [ ! -f metadyn.tpr ] || [ ! -f metadyn.trr ]; then
  echo "- Manca metadyn.tpr o metadyn.trr in $(pwd), esco"
  exit 1
fi

# Limite impostato direttamente a 80 ns (80000 ps) senza eseguire gmx check
total_ps=${max_ps}

echo "Traiettoria: elaborazione impostata fino a ${total_ps} ps (80 ns)"

mkdir -p chunks
touch .chunks_done

b=0
while (( $(echo "$b < $total_ps" | bc -l) )); do
  e=$(echo "$b + $chunk_ps" | bc)
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
    -pbc mol -b "${b}" -e "${e}" > "chunks/log_${chunk_label}.out" 2>&1 << STOP
${outGrp}
STOP
  if [ -s "${chunk_file}" ]; then
    echo "${chunk_label}" >> .chunks_done
    echo "  -> completato (${chunk_file})"
  else
    echo "  -> ATTENZIONE: chunk ${chunk_label} fallito. Controlla chunks/log_${chunk_label}.out"
    exit 1
  fi

  b=$e
done

echo "Tutti i chunk fino a ${total_ps} ps sono stati generati."
