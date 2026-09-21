#!/bin/bash
#SBATCH --job-name=trjconv_recover
#SBATCH --partition=lrd_all_serial
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=IscrC_hDDC
#SBATCH --output=trjconv_recover_%A_%a.out
#SBATCH --array=0-5

# Estrae, per ciascuna replica metadyn_${SLURM_ARRAY_TASK_ID}, la traiettoria
# protein-only con fix PBC, troncata a 80 ns, riprendendo esattamente il passo
# di startme.sh che non e' mai stato eseguito (mdrun si e' fermato tra 80 e 90 ns
# prima di arrivare al blocco trjconv finale).
#
# lrd_all_serial e' budget-free (non consuma core-hours dal Project Account)
# ed e' pensata proprio per pre/post-processing: qui non serve la GPU perche'
# trjconv/convert-tpr sono operazioni CPU-only.
#
# Sottomissione: sbatch extract_prot_pbc_80ns.sh   (dalla cartella 08-metadyn/)

set -e

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

outGrp=1
END_TIME_PS=80000   # 80 ns espressi in ps, unita' richiesta da trjconv -e

REP_DIR="metadyn_${SLURM_ARRAY_TASK_ID}"

if [ ! -d "${REP_DIR}" ]; then
    echo "Cartella ${REP_DIR} non trovata, esco."
    exit 0
fi

cd "${REP_DIR}"

if [ ! -f metadyn.trr ]; then
    echo "[${REP_DIR}] metadyn.trr assente, salto questa replica."
    exit 0
fi

if [ -f metadyn_prot_pbc.trr ]; then
    echo "[${REP_DIR}] metadyn_prot_pbc.trr gia' presente, salto."
    exit 0
fi

if [ ! -f metadyn_prot.tpr ]; then
    echo "[${REP_DIR}] Creo metadyn_prot.tpr..."
    gmx_plumed convert-tpr -s metadyn.tpr -o metadyn_prot.tpr << STOP >& tpr-convert.out
${outGrp}
STOP
fi

echo "[${REP_DIR}] Estrazione traiettoria protein-only troncata a 80 ns..."
gmx_plumed trjconv -s metadyn.tpr -f metadyn.trr -o metadyn_prot_pbc.trr \
    -pbc mol -e ${END_TIME_PS} << STOP >& trjconv.out
${outGrp}
STOP

echo "[${REP_DIR}] Completato: metadyn_prot_pbc.trr pronto."

# Decommentare per liberare spazio una volta verificato che l'estrazione e'
# andata a buon fine (il .trr all-atom puo' pesare centinaia di GB):
# rm -f metadyn.trr
