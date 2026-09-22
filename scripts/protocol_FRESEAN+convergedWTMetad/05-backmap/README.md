# Step 05 - Backmapping

This step re-projects ("backmaps") the extracted coarse-grained vibrational modes back onto the full all-atom representation. To maintain efficiency, this transformation is executed exclusively on the selected reference replica (R2) rather than across all five replicas.

## Input

Reference all-atom structure (ref.pdb) from 03-CG/R2, coarse-grained mode file (.xyz) from 04-FRESEAN/R2, and the backmapping script (backmap.py).

## Output

Individual all-atom mapped mode files (evec_${mode}_aa_backmapped.xyz) and their concatenated single multi-frame file (evec_aa_1-30_backmapped.xyz).

