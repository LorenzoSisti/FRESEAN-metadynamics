#!/bin/bash

#SBATCH --job-name em_equi_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=04:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

export OMP_NUM_THREADS=8
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#BEGIN INPUT
inpTOP=../00-prep/topol.top
inpGRO=../00-prep/prep.gro
next=1
nextDir=../02-MD
gmx=gmx_plumed
set -e
#END INPUT

files=(
em.mdp
equi.mdp
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

#Run energy minimization of the freshly solvated protein
if [ ! -f em/em.gro ]; then
if [ ! -d em ]; then
mkdir em
fi
cd em
$gmx grompp -f ../em.mdp -c ../${inpGRO} -p ../${inpTOP} -r ../${inpGRO} -o em.tpr >& grompp.out
mpirun -np 1 $gmx mdrun -v -deffnm em -cpi -ntomp 8 -pin on >& mdrun.out
cd ..
else
echo "em/em.gro già presente, salto EM"
fi

#Run equilibration MD simulation with position restraints on heavy protein atoms
if [ ! -d equi ]; then
mkdir equi
fi
cd equi
$gmx grompp -f ../equi.mdp -c ../em/em.gro -p ../${inpTOP} -r ../${inpGRO} -o equi.tpr >& grompp.out
mpirun -np 4 $gmx mdrun -v -deffnm equi -cpi -ntomp 8 -nb gpu -pme gpu -npme 1 -pin off -nstlist 500 >& mdrun.out
cd ..

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f equi/equi.gro ]; then
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
