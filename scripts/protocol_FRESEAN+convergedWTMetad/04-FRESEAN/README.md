# Step 04 - FRESEAN Analysis

This step performs the core FRESEAN vibrational analysis on the coarse-grained trajectories across all 5 independent replicas (R1–R5). It computes velocity cross-correlation matrices, diagonalizes them in the frequency domain, extracts low-frequency vibrational modes, and quantitatively compares modes between replicas to evaluate reproducibility and select a reference trajectory.

## Input

Coarse-grained trajectories (.trr) and CG topologies (.mtop) from 03-CG.

## Output

Frequency-domain correlation matrices, extracted low-frequency mode files (.xyz), and inter-replica mode correlation/overlap analysis results.

# Note to myself

Classifica per correlazione media di sottospazio 3D (modi 7,8,9) con le altre 4 repliche: R2 = 0.233 (migliore), R5 = 0.228, R1 = 0.227, R3 = 0.186, R4 = 0.104 (peggiore, con una coppia R1-R4 quasi nulla: 0.004). Replica scelta come riferimento per tutti gli step successivi: R2.

