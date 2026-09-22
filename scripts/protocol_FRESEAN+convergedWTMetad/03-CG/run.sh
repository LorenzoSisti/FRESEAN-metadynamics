#!/bin/bash
#SBATCH --job-name 03CG_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod
#SBATCH --array=1-5

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

#BEGIN INPUT
gmx=gmx_plumed

#GMX topology for protein (the same across all replicas, you find it in 00-prep)
inpTOPprot=../../00-prep/topol_prot.top

inpEMmdp=../../00-prep/em.mdp
inpBOXgro=../../00-prep/box.gro

#job file for fresean coarse (the same across all replicas)
inpJOB=../static.job

#GMX files from MD sampling (specific for each replica R${SLURM_ARRAY_TASK_ID})
inpTPR=../../02-MD/R${SLURM_ARRAY_TASK_ID}/sample-NPT.tpr
inpTRR=../../02-MD/R${SLURM_ARRAY_TASK_ID}/sample-NPT.trr

#output trajectory for protein with PBC fixes (in the replicas subfolders)
outTRRprotAA=sample-NPT_prot_pbc.trr

#output trajectory for protein in CG representation
outTRRprotCG=sample-NPT_prot_pbc-CG.trr

#GMX group number for protein (gmx trjconv)
outGrp=1

#time between frames in input trajectory
TRJtimestep=0.020

#number of steps in input trajectory
TRJframes=1000000

next=0
nextDir=../../04-FRESEAN
#END INPUT

set -e

mkdir -p R${SLURM_ARRAY_TASK_ID}
cd R${SLURM_ARRAY_TASK_ID}

files=(
${inpTOPprot}
${inpEMmdp}
${inpBOXgro}
${inpJOB}
${inpTPR}
${inpTRR}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Flatten the topology (resolve #include for multi-chain systems, e.g. homodimers)
#before feeding it to fresean, which does not follow #include directives itself
$gmx grompp -f ${inpEMmdp} -c ${inpBOXgro} -p ${inpTOPprot} \
  -pp topol_prot_full.top -o pp.tpr -maxwarn 10 >& grompp_pp.out

#Generate a all-atom topology file (custom format) for the protein
fresean mtop << STOP >& mtop.out
topol_prot_full.top
topol_prot_full.top
topol_prot_full.top
topol_prot-aa.mtop
STOP

#Make the protein whole (fix PBC jumps) and write out protein-only trajectory
$gmx trjconv -s ${inpTPR} -f ${inpTRR} -o ${outTRRprotAA} -pbc mol << STOP >& trjconv1.out
${outGrp}
STOP

#Generate all-atom reference structure for rot+trans fitting
$gmx trjconv -s ${inpTPR} -f ${outTRRprotAA} -o ref.gro -e 0.0 << STOP >& trjconv2.out
${outGrp}
STOP

#some postprocessing for PLUMED compatibility of ref.pdb
$gmx editconf -f ref.gro -o ref.pdb
tail -n +5 ref.pdb > tmp.pdb
head -n -2 tmp.pdb > ref.pdb
sed -i '1i REMARK TYPE=OPTIMAL' ref.pdb
sed -i '$aEND' ref.pdb
rm tmp.pdb

#Generate input parameter file for coarse-graining
#WARNING: update expected soon
cat << STOP >& coarse.inp
#fnTop
topol_prot-aa.mtop
#fnCrd
${outTRRprotAA}
#fnVel (if format xyz,crd,dcd)
#fnJob
${inpJOB}
#grp
0
#nRead
${TRJframes}
#analysisInterval
1
#fnOutTraj
tmptraj.gro
#fnOutTopol
topol_prot-cg.mtop
STOP

#Apply coarse-graining to all-atom protein trajectory
fresean coarse -f coarse.inp >& coarse.out

#Convert coarse-grained trajectory from gro (ASCII text) to trr (binary) format
$gmx trjconv -s tmptraj.gro -f tmptraj.gro -timestep ${TRJtimestep} -o ${outTRRprotCG} << STOP >& trjconv3.out
0
STOP

#Generate coarse-grained reference structure for rot+trans fitting
nAtoms=`head -n 2 tmptraj.gro | tail -n 1`
nLines=${nAtoms}
((nLines +=3))
head -n ${nLines} tmptraj.gro >& ref-cg.gro
rm tmptraj.gro

echo "03-CG replica R${SLURM_ARRAY_TASK_ID} completata: ${outTRRprotCG} pronto per 04-FRESEAN"
cd ..

#Start next part of the project if next=1
#NB: next=0 di default qui perché 04-FRESEAN deve processare tutte e 5
#    le repliche insieme (matrice di correlazione dei modi); si lancia
#    manualmente solo dopo che tutte e 5 le istanze del job array sono
#    completate con successo, non automaticamente dalla prima che finisce.
if [ ${next} -eq 1 ]; then
  if [ -f R${SLURM_ARRAY_TASK_ID}/${outTRRprotCG} ]; then
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
