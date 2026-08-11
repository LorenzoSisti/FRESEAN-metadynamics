#!/bin/bash
#SBATCH --job-name 07resample_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

export OMP_NUM_THREADS=8
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

#BEGIN INPUT
# 07-resample e' allo stesso livello di 00-prep/01-em+equi (sibling,
# non annidata di un livello extra): un solo ".." basta, non "../..".
workingDir=..
inpTOP=$(readlink -f ${workingDir}/00-prep/topol.top)
inpGRO=$(readlink -f ${workingDir}/01-em+equi/equi/equi.gro)
gmx=gmx_plumed
maxh=23.5   # margine di sicurezza sotto --time=24:00:00, per checkpoint pulito
#END INPUT

set -e

files=(
${inpTOP}
${inpGRO}
sample-states.mdp
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

# Crea la cartella solo alla primissima submission
mkdir -p run-NPT
cd run-NPT

if [ ! -f sample-states.tpr ]; then
  $gmx grompp -f ../sample-states.mdp -c ${inpGRO} -p ${inpTOP} \
    -o sample-states.tpr -maxwarn 1 >& grompp.out
fi

# Continua da checkpoint se presente (run precedente andata in timeout),
# altrimenti parti da zero
if [ -f sample-states.cpt ]; then
  mpirun -np 4 $gmx mdrun -v -deffnm sample-states -cpi sample-states.cpt \
    -ntomp 8 -nb gpu -pme gpu -npme 1 -pin off -nstlist 500 -maxh ${maxh} \
    >> mdrun.out 2>&1
else
  mpirun -np 4 $gmx mdrun -v -deffnm sample-states \
    -ntomp 8 -nb gpu -pme gpu -npme 1 -pin off -nstlist 500 -maxh ${maxh} \
    >> mdrun.out 2>&1
fi

cd ..

# Controlla se la simulazione ha raggiunto nsteps; se no, ri-sottometti
target_nsteps=$(grep -i "^nsteps" sample-states.mdp | awk -F'=' '{print $2}' | awk '{print $1}')
current_step=$(tail -50 run-NPT/sample-states.log | grep -A1 "^ *Step " | tail -1 | awk '{print $1}')

if [ -z "${current_step}" ] || [ "${current_step}" -lt "${target_nsteps}" ]; then
  echo "Simulazione incompleta (step ${current_step:-0}/${target_nsteps}), ri-sottometto..."
  sbatch --dependency=afterany:${SLURM_JOB_ID} run.sh
  exit 0
fi

echo "07-resample: 100 ns completati (${current_step}/${target_nsteps} step)"

# Pull configuration every 5 ns (solo una volta che la simulazione e' completa)
mkdir -p snapshots
$gmx trjconv -s run-NPT/sample-states.tpr -f run-NPT/sample-states.trr -pbc mol -sep -dt 5000.0 -ndec 8 -o snapshots/state_.gro <<END
0
END

echo "07-resample completato: snapshot ogni 5 ns pronti in snapshots/"
