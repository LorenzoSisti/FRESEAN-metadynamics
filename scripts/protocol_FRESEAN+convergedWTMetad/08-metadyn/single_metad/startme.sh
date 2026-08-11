#!/bin/bash
#SBATCH --job-name METAD_0
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

export OMP_NUM_THREADS=8
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# Disabilita i backup automatici di GROMACS: con l'auto-resubmit sotto,
# ogni tentativo che ritrova la simulazione gia' completa (comune, dato
# che mdrun esce quasi subito se lo step nel checkpoint >= nsteps) genera
# comunque un file di output finale, e GROMACS per default rinomina il
# precedente come backup invece di sovrascriverlo. Dopo ~99 tentativi
# questo causa un errore fatale ("Won't make more than 99 backups"),
# gia' incontrato in 07-resample. -1 disabilita del tutto i backup.
export GMX_MAXBACKUP=-1

#BEGIN INPUT
# metadyn_${i}/ e' due livelli sotto la root del protocollo
# (08-metadyn/metadyn_${i}/), quindi "../.." e' corretto qui - ma lo
# rendiamo comunque assoluto con readlink -f per robustezza futura,
# nel caso lo script venga esteso con ulteriori "cd" in seguito.
workingDir=../..
inpTOP=$(readlink -f ${workingDir}/00-prep/topol.top)
inpGRO=$(readlink -f ${workingDir}/07-resample/snapshots/state_${replica}.gro)
inpPlumedPDB=$(readlink -f ${workingDir}/06-ModeProj/plumed-mode-input.pdb)
gmx=gmx_plumed
maxh=23.5
outGrp=1
#END INPUT

set -e

files=(
${inpTOP}
${inpGRO}
${inpPlumedPDB}
metadyn.mdp
plumed-mode-metadyn.dat
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Copy PLUMED input file with reference structure and FRESEAN modes to standardized file name
cp ${inpPlumedPDB} plumed-mode-input.pdb

if [ ! -f metadyn.tpr ]; then
#Create tpr input file for WT-metadynamics simulation
$gmx grompp -f metadyn.mdp -c ${inpGRO} -p ${inpTOP} -o metadyn.tpr -maxwarn 1 >& grompp.out
fi

if [ ! -f metadyn_prot.tpr ]; then
#Create tpr file with only protein atoms for analysis
$gmx convert-tpr -s metadyn.tpr -o metadyn_prot.tpr << STOP >& tpr-convert.out
${outGrp}
STOP
fi

#Run the WT-metadynamics simulation with GROMACS & PLUMED.
#Se e' gia' completa (checkpoint con step >= nsteps), mdrun esce quasi
#subito senza ricalcolare nulla - questo rende sicuro il resubmit ripetuto.
if [ -f metadyn.cpt ]; then
  mpirun -np 1 $gmx mdrun -v -deffnm metadyn -cpi metadyn.cpt -ntomp 8 \
    -nb gpu -pme gpu -pin off -maxh ${maxh} -plumed plumed-mode-metadyn.dat \
    >> mdrun.out 2>&1
else
  mpirun -np 1 $gmx mdrun -v -deffnm metadyn -ntomp 8 \
    -nb gpu -pme gpu -pin off -maxh ${maxh} -plumed plumed-mode-metadyn.dat \
    >> mdrun.out 2>&1
fi

# Controlla se la simulazione ha raggiunto nsteps; se no, ri-sottometti.
target_nsteps=$(grep -i "^nsteps" metadyn.mdp | awk -F'=' '{print $2}' | awk '{print $1}')
current_step=$(tail -50 metadyn.log | grep -A1 "^ *Step " | tail -1 | awk '{print $1}')

if [ -z "${current_step}" ] || [ "${current_step}" -lt "${target_nsteps}" ]; then
  echo "Simulazione incompleta (step ${current_step:-0}/${target_nsteps}), ri-sottometto..."
  sbatch --job-name=METAD_REPLICA_${replica}.run --export=replica=${replica} \
    --dependency=afterany:${SLURM_JOB_ID} startme.sh
  exit 0
fi

echo "metadyn replica ${replica}: 100 ns completati (${current_step}/${target_nsteps} step)"

# kT in kJ/mol can be determined by running `plumed kT --temp 300`
kt=2.494339
#Generate free energy surface in FRESEAN mode space
plumed sum_hills --min -0.02,-0.02 --max 0.02,0.02 --bin 200,200 --hills plumed-mode-metadyn.hills --outfile plumed-mode-metadyn.fes --mintozero --kt 2.494339 --stride 5000

#Write out metadynamics trajectory for only protein atoms after PBC fix
$gmx trjconv -s metadyn.tpr -f metadyn.trr -o metadyn_prot_pbc.trr -pbc mol << STOP >& trjconv.out
${outGrp}
STOP

# Libera spazio: il .trr all-atom completo (proteina+acqua+ioni, centinaia
# di GB) non serve piu' una volta estratta la traiettoria protein-only.
rm -f metadyn.trr

echo "metadyn replica ${replica} completata: metadyn_prot_pbc.trr e plumed-mode-metadyn.fes pronti"
