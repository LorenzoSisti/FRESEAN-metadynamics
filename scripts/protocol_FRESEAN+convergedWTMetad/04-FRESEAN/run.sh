#!/bin/bash
#SBATCH --job-name 04FRESEAN_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod
#SBATCH --array=1-5

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

#BEGIN INPUT

#input topology in custom format (specifica della replica R${SLURM_ARRAY_TASK_ID})
inpMTOP=../../03-CG/R${SLURM_ARRAY_TASK_ID}/topol_prot-cg.mtop

#input GMX trajectory file (specifica della replica R${SLURM_ARRAY_TASK_ID})
inpTRR=../../03-CG/R${SLURM_ARRAY_TASK_ID}/sample-NPT_prot_pbc-CG.trr

#input reference file rot trans+rot fitting (specifica della replica R${SLURM_ARRAY_TASK_ID})
inpRefGROcg=../../03-CG/R${SLURM_ARRAY_TASK_ID}/ref-cg.gro

#file di job condiviso, sta in 04-FRESEAN
inpJOB=../static.job

#number of frames in trajectory
nFrames=1000000

#length of correlation functions in frames
nCorr=100

#stem name of output files
outName=cg

#trajectory time step (picoseconds)
TRJtimestep=0.020

next=0
nextDir=../../05-backmap
#END INPUT

set -e

mkdir -p R${SLURM_ARRAY_TASK_ID}
cd R${SLURM_ARRAY_TASK_ID}

files=(
${inpMTOP}
${inpTRR}
${inpRefGROcg}
${inpJOB}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Generate input parameter file for FRESEAN mode analysis
cat << STOP >& matrix.inp
#fnTop (custom .mtop file)
${inpMTOP}
#fnCrd (CP2K position and velocity files)
${inpTRR}
#fnJob (Job file defining atom groups)
${inpJOB}
#nRead (Number of frames in the file)
${nFrames}
#analysisInterval (Analysis should be performed every X frames)
1
#fnRef (Reference file generated from genREF.exe)
${inpRefGROcg}
#alignGrp
0
#analyzeGrp
0
#wrap
0
#nCorr (correlation time in trajecotry steps)
${nCorr}
#winSigma (wavenumbers cm^-1)
10.0
#binaryMatrix
1
#doGenModes
0
#convergence (only for generalized normal modes)
1.0e-5
#maxIter (only for generalized normal modes)
100
#fnOut (Output file appender)
${outName}
STOP

#Perform FRESEAN mode analysis for input trajectory
#Here: we use 32 CPU cores --> see SBATCH flag: --cpus-per-task
export OMP_NUM_THREADS=32
fresean matrix -f matrix.inp >& matrix.out

#Diagonalize velocity cross correlation matrix at all sampled frequencies
fresean eigen -m matrix_${outName}.mmat -n $nCorr >& eigen.out

#Generate input parameter file to extract FRESEAN modes sampled at zero frequency
cat << STOP >& extract.inp
#fnEigVec (file containing eigenvectors in binary format)
evec_matrix_${outName}.mmat
#extractMode ( Mode 0 -> freqSel is in wavenumbers; Mode 1 -> freqSel is matrix index )
1
#freqSel (extract mode 0 -> Integer frequency in cm^-1;  mode 1 -> frequency index)
1
#trrFreq
${TRJtimestep}
#modeStart (First vibrational mode to read)
1
#modeEnd (Last vibrational mode to read)
30
#fnOut ( Output file appender )
cg.xyz
STOP

#Extract zero frequency modes
fresean extract -f extract.inp >& extract.out

echo "04-FRESEAN replica R${SLURM_ARRAY_TASK_ID} completata"
cd ..

#Start next part of the project if next=1
#NB: next=0 di default: 05-backmap (e prima ancora, il confronto a matrice
#    di correlazione fra le 5 repliche) deve avvenire solo dopo che TUTTE
#    e 5 le istanze del job array sono completate, non alla prima che finisce.
if [ ${next} -eq 1 ]; then
  if [ -f R${SLURM_ARRAY_TASK_ID}/evec_freq1_mode1-30_cg.xyz ]; then
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
