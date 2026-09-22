# Step 06 - Mode Projection

This step isolates specific low-frequency vibrational modes (modes 7 and 8) chosen as collective variables (CVs) for downstream metadynamics. It converts these modes into PLUMED-compatible reference inputs and projects the unbiased trajectory onto them to measure natural mode fluctuations, which are then used to calibrate metadynamics parameters.

## Input

Backmapped all-atom mode trajectory (evec_aa_1-30_backmapped.xyz), reference structure (ref.pdb), topology (.mtop), and unbiased trajectory (.trr) from R2.

## Output

PLUMED CV reference structure (plumed-mode-input.pdb), mode projection output (plumed-mode-projection.out), and calculated fluctuation standard deviations (standard-deviation.out).


