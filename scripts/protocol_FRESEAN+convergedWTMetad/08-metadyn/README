# Step 08 - Well-tempered Metadynamics

This step executes 20 independent production replicas of Well-Tempered Metadynamics (WT-metaD) biased along vibrational modes 7 and 8 (derived from reference replica R2) as collective variables (CVs). The goal is to thoroughly sample phase space and reconstruct the free energy surface (FES) governing the conformational transition.

## Input

System topology (topol.top), initial conformational snapshots (state_${i}.gro), PLUMED CV reference structure (plumed-mode-input.pdb), and calibrated PLUMED parameter files.

## Output

Metadynamics raw simulation files, deposited bias potentials (plumed-mode-metadyn.hills), and protein-only PBC-corrected trajectories for each of the 20 replicas.

# How to execute the codes in this folder via HPC CINECA Leonardo supercomputer

run.sh in 08-metadyn/ is NOT a SLURM job: it must be executed using bash run.sh (never sbatch run.sh). It is a launcher script that runs on the login node and submits 19–20 actual jobs via sbatch—one for each replica (each containing the proper SLURM header inside startme.sh).

In plumed-mode-metadyn.dat, the SIGMA parameter for the two collective variables (modes 7 and 8) must be set to HALF of the standard deviation of the unbiased fluctuations calculated in 06-ModeProj (found in the standard-deviation.out file); it is not an arbitrary value. For the apo branch: SIGMA=4.4e-05,1.6e-04 (calculated from mode 7 std = 8.78e-05, mode 8 std = 3.21e-04). For the holo branch, it must be recalculated from scratch using the same criteria, based on the standard-deviation.out values produced by the 06-ModeProj step of the holo system.

## How to run 08-metadyn

Prerequisiti (devono essere gia' completati e verificati):
- 06-ModeProj completato per la replica scelta come riferimento (fornisce
  plumed-mode-input.pdb e le deviazioni standard in standard-deviation.out)
- 07-resample completato (fornisce i 20 snapshot in snapshots/state_0.gro
  ... state_19.gro, uno ogni 5 ns su 100 ns di dinamica non biased)

Passi:

1. Dentro single_metad/, verificare/impostare SIGMA in
   plumed-mode-metadyn.dat a META' della deviazione standard delle
   fluttuazioni non biased calcolate in 06-ModeProj (vedi nota sopra).
   NON lasciare mai il placeholder SIGMA=XXX,YYY, PLUMED fallirebbe
   al parsing.

2. Test preliminare su UNA sola replica, prima di lanciarle tutte:

     cd 08-metadyn
     cp -r single_metad metadyn_0
     cd metadyn_0
     sbatch --job-name=METAD_REPLICA_0.run --export=replica=0 startme.sh

   Verificare che grompp.out, mdrun.out e plumed-mode-metadyn.out non
   mostrino errori e che il bias (colonna metad.bias) cresca in modo
   sensato prima di procedere.

3. Solo dopo aver confermato che la replica di test funziona, lanciare
   le restanti (modificare run.sh impostando l'indice di partenza del
   ciclo a 1, dato che metadyn_0 e' gia' stata creata manualmente):

     cd 08-metadyn
     bash run.sh          # ATTENZIONE: mai "sbatch run.sh"

   run.sh NON e' un job SLURM: e' uno script lanciatore che gira sul
   login node e sottomette a sua volta 19 job veri (uno per replica)
   tramite sbatch, ciascuno con l'header SLURM corretto definito dentro
   startme.sh.

4. Ogni replica gira con auto-resubmit: se il --time allocato scade
   prima di raggiungere nsteps, lo script si ri-sottomette da solo via
   sbatch --dependency=afterany, riprendendo dal checkpoint. Non serve
   intervento manuale in condizioni normali.

5. Monitorare periodicamente lo stato con:

     squeue --me
     cindata          # spazio disco
     saldo -b          # budget ore-GPU del progetto

   ATTENZIONE: se il budget del progetto (saldo -b) supera il 100%,
   sbatch rifiuta silenziosamente i nuovi resubmit con errore "invalid
   account or expired budget", visibile solo aprendo lo slurm-*.out
   del tentativo fallito. Le repliche restano ferme al loro ultimo
   checkpoint valido (nessun dato perso, ma nessun avanzamento) finche'
   il budget non viene reintegrato.

Le sottocartelle metadyn_0 ... metadyn_19 sono escluse da git
(.gitignore) perche' contengono output di simulazione molto pesanti
(centinaia di GB a replica): sono istanze di esecuzione di
single_metad/, non varianti di codice. single_metad/ resta l'unica
fonte di verita' versionata per gli script di questo step.
