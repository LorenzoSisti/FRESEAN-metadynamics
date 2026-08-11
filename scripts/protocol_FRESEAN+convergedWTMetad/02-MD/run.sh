#!/bin/bash
#SBATCH --job-name 02MD_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod
#SBATCH --array=1-5

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

export OMP_NUM_THREADS=8
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#BEGIN INPUT
inpTOP=../../00-prep/topol.top
inpGRO=../../01-em+equi/equi/equi.gro
next=0
nextDir=../../03-CG
gmx=gmx_plumed
set -e
#END INPUT

mkdir -p R${SLURM_ARRAY_TASK_ID}
cd R${SLURM_ARRAY_TASK_ID}

files=(
../sample-NPT.mdp
${inpTOP}
${inpGRO}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Run the equilibrium MD simulation to sample vibrational modes
#Here: 1 GPU and 8 CPU cores per replica
$gmx grompp -f ../sample-NPT.mdp -c ${inpGRO} -p ${inpTOP} -o sample-NPT.tpr -maxwarn 1 >& grompp.out
mpirun -np 4 $gmx mdrun -v -deffnm sample-NPT -cpi -ntomp 8 -nb gpu -pme gpu -npme 1 -pin off -nstlist 500 >& mdrun.out

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f sample-NPT.gro ]; then
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
