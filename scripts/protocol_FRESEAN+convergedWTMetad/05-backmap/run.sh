#!/bin/bash
#SBATCH --job-name 05backmap_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

############################
# CHECK NUMPY INSTALLATION #
############################
# Niente mamba/conda su Leonardo (nessun modulo disponibile in nessun
# profilo verificato): usiamo il python3 di sistema, installando numpy
# per l'utente con pip3 --user se non gia' presente.

python3 -c "import numpy" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "numpy non trovato, tento l'installazione con pip3 --user..."
    pip3 install --user numpy >& pip_install_numpy.out
    python3 -c "import numpy" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo "Installazione di numpy fallita. Controlla pip_install_numpy.out"
        echo "-exiting"
        exit 1
    fi
    echo "numpy installato correttamente."
fi


############################
#        BEGIN INPUT       #
############################

#Replica scelta come riferimento (dalla correlazione di sottospazio calcolata
#in 04-FRESEAN con correlate_modes.py: R2 aveva la media piu' alta con le altre 4)
REF=2

#all-atom reference structure (must match CG reference)
inpRefPDBaa=../03-CG/R${REF}/ref.pdb

#coarse-grained zero-frequency eigenvectors
inpModesXYZ=../04-FRESEAN/R${REF}/evec_freq1_mode1-30_cg.xyz

# 20 ns unbiased trajectory (protein-only, PBC-fixed; scritta da 03-CG, non da 02-MD)
inpTRR=../03-CG/R${REF}/sample-NPT_prot_pbc.trr

# If next=1 -> start next step automatically
# If next=0 -> terminate after this step
next=0
nextDir=../06-ModeProj

# If ERROR_FLAG=0 -> no error
ERROR_FLAG=0

############################
#   CHECK FILE EXISTANCE   #
############################

files=(
${inpRefPDBaa}
${inpModesXYZ}
${inpTRR}
backmap.py
)

for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

###########################
#  BEGIN MODE CONVERSION  #
###########################

cp ${inpRefPDBaa} ref.pdb
cp ${inpModesXYZ} input-cg-modes.xyz

reflinecount=`cat ref.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms} '`
echo "Found ${reflinecount} atoms in ref.pdb"

# Convert FRESEAN modes 1-30 at zero frequency (backmap all-atom)
minmode=1
maxmode=30

for ((mode=${minmode}; mode<=${maxmode}; mode++))
do
python3 backmap.py $mode >& prep_plumed_mode${mode}.out
modelinecount=`cat evec_${mode}_aa_backmapped.xyz | awk 'BEGIN {natoms=0} $1 ~ "X" {natoms++} END {print natoms}'`
echo "Found ${modelinecount} atoms in evec_${mode}_aa_backmapped.xyz"
if [ ${modelinecount} != ${reflinecount} ]; then
        echo "Mode ${mode} (evec_${mode}_aa_backmapped.xyz) is incomplete."
        ERROR_FLAG=1
fi

done

if [ ${ERROR_FLAG} == 1 ];then
	exit
fi

for i in $(seq $minmode $maxmode); do
    cat evec_${i}_aa_backmapped.xyz
done > evec_aa_${minmode}-${maxmode}_backmapped.xyz

echo "05-backmap completato usando la replica di riferimento R${REF}"

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
    if [ -f evec_aa_${minmode}-${maxmode}_backmapped.xyz ]; then
        if [ -d ${nextDir} ]; then
            curDir=`pwd`
            cd ${nextDir}
            if [ -f run.sh ]; then
                sbatch run.sh
            fi
            cd ${curDir}
        fi
    fi
fi
