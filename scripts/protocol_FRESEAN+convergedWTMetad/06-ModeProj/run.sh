#!/bin/bash
#SBATCH --job-name 06ModeProj_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

############################
# CHECK NUMPY INSTALLATION #
############################
# Niente mamba/conda su Leonardo: usiamo il python3 di sistema,
# installando numpy per l'utente con pip3 --user se non gia' presente.

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

#Replica scelta come riferimento (stessa di 05-backmap: R2)
REF=2

#all-atom reference structure (must match CG reference)
inpRefPDBaa=../03-CG/R${REF}/ref.pdb

#coarse-grained zero-frequency eigenvectors, backmapped all-atom (da 05-backmap)
inpModesXYZ=../05-backmap/evec_aa_1-30_backmapped.xyz

#custom all-atom topology (da 03-CG, generata da fresean mtop)
inpMTOPaa=../03-CG/R${REF}/topol_prot-aa.mtop

# 20 ns unbiased trajectory (protein-only, PBC-fixed; da 03-CG, non da 02-MD)
inpTRR=../03-CG/R${REF}/sample-NPT_prot_pbc.trr

# If next=1 -> start next step automatically
# If next=0 -> terminate after this step
next=0
nextDir=../07-resample

# If ERROR_FLAG=0 -> no error
ERROR_FLAG=0

############################
#   CHECK FILE EXISTANCE   #
############################

files=(
${inpRefPDBaa}
${inpModesXYZ}
${inpMTOPaa}
${inpTRR}
prep_plumed.py
plumed-mode-projection.dat
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
cp ${inpMTOPaa} topol_prot-aa.mtop

# Convert FRESEAN modes 7 & 8 at zero frequency
mode=7
python3 prep_plumed.py input-cg-modes.xyz topol_prot-aa.mtop $mode ref.pdb >& prep_plumed_mode${mode}.out
mode=8
python3 prep_plumed.py input-cg-modes.xyz topol_prot-aa.mtop $mode ref.pdb >& prep_plumed_mode${mode}.out

###########################
#  CHECK PREP COMPLETION  #
###########################

reflinecount=`cat ref.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms} '`
mode7linecount=`cat evec_7_aa_scaled.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms}'`
mode8linecount=`cat evec_8_aa_scaled.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms}'`
echo "Found ${reflinecount} atoms in ref.pdb"
echo "Found ${mode7linecount} atoms in evec_7_aa_scaled.pdb"
echo "Found ${mode8linecount} atoms in evec_8_aa_scaled.pdb"

if [ ${mode7linecount} != ${reflinecount} ]; then
	echo "Mode 7 (evec_7_aa_scaled.pdb) is incomplete."
	ERROR_FLAG=1
fi

if [ ${mode8linecount} != ${reflinecount} ]; then
	echo "Mode 8 (evec_8_aa_scaled.pdb) is incomplete."
	ERROR_FLAG=1
fi

if [ ${ERROR_FLAG} == 1 ];then
	exit
fi

#Combine extracted modes with all-atom reference input PLUMED input
cat ref.pdb evec_7_aa_scaled.pdb evec_8_aa_scaled.pdb > plumed-mode-input.pdb
rm ref.pdb input-cg-modes.xyz evec_7_aa_scaled.pdb evec_8_aa_scaled.pdb

#Analyze displacement fluctuations along FRESEAN modes
plumed driver --mf_trr ${inpTRR} --plumed plumed-mode-projection.dat --kt 2.494339 > plumed-driver.out

#number of histogram bins
bins=200

#print standard deviations of displacement fluctuations along FRSEAN modes
python3 standard-deviation.py plumed-mode-projection.out $bins >& standard-deviation.out

echo "06-ModeProj completato usando la replica di riferimento R${REF}"

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f plumed-mode-projection.dat ]; then
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
