# Step 07 - Resampling

This step executes an extended 100 ns unbiased molecular dynamics simulation starting from the equilibrated structure to sample conformational space. From this trajectory, 21 distinct snapshots are extracted at regular 5 ns intervals to serve as starting configurations for downstream production metadynamics replicas.

## Input

System topology (topol.top) from 00-prep, equilibrated structure (equi.gro) from 01-em+equi, and simulation parameter file (sample-states.mdp).

## Output

100 ns unbiased MD trajectory and log files, along with 21 extracted conformational snapshots (state_0.gro to state_20.gro).


