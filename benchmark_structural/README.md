# To generate figure 2BCDEFG
## Structural benchmark

Step 0:
---
Source data:
 - 496 PDB chain sequences with structural tandem repeats (average TM-score > 0.7) detected by CE-Symm

Files:
 - input.aa.fasta
 - answer.tsv

Step 1:
---
Run TR detecting tools and parsing the results

Notes:
 - Analysis_structural_benchmark.txt

Scripts:
 - parse_HHrepID.2.r
 - parse_MREPS.1.py
 - parse_RADAR.1.py
 - parse_TRF.1.py
 - parse_treks.3.py
 - parse_trust.4.py
	 
Step 2:
---
Compute TPR and FPR and plot

Notes:
 - Analysis_structural_benchmark.txt

Scripts:
 - eval_ROC.3.r


Show an example
---

Notes:
 - Analysis_4CNL_example.txt