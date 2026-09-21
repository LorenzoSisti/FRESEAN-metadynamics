#!/bin/bash
for dir in metadyn_*/; do
  i="${dir%/}"
  echo ""
  echo "=== Inizio ${i} ==="
  cd "${i}" || continue
  bash ../extract_prot_pbc_chunked.sh
  cd ..
done

echo ""
echo "Giro completato. Rilancia lo stesso script per proseguire le"
echo "repliche eventualmente interrotte a meta'."
