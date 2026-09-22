#!/bin/bash
#SBATCH --job-name prep_DDC_apo_open
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --account=IscrC_hDDC
#SBATCH --partition=boost_usr_prod

source /leonardo_scratch/large/userexternal/lsisti00/env-plumed.sh

#BEGIN INPUT
gmx=gmx_plumed
pdb=DDC_homo_open.pdb
box_d=1.0          # nm, giving minimal buffer for editconf -d
salt=0.15          # target saline concentration (mol/L)
#END INPUT

set -e   # stop the simulation if an error occurs

#Generate first topology file (protein + crystal water/ions)
$gmx pdb2gmx -f ${pdb} -p topol_prot.top -o prot.gro -ff amber99sb-ildn -water tip3p

#Define simulation box for preiodic boundary conditions (PBC)
$gmx editconf -f prot.gro -o box.gro -c  -d 2.0 -bt triclinic

#Make copy of protein topology file and delete crystal water/ions from original
cp topol_prot.top topol.top
awk '{if($1!="SOL"&&$1!="NA"&&$1!="CL") {printf("%s\n",$0);}}' topol.top >& topol_prot.top

#Solvate the protein and update new topology file (keep crystal water/ions)
$gmx solvate -cp box.gro -cs -p -o solv.gro
rm \#topol.top.1\#

#Add ions to the solution (requires *.tpr input file)
$gmx grompp -f em.mdp -c solv.gro -p -o tmp.tpr -maxwarn 1
$gmx genion -s tmp.tpr -p -o prep.gro -neutral -conc ${salt} << STOP
SOL
STOP
rm \#topol.top.1\#
