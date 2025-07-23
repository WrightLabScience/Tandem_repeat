# To generate figure 2BCDEFG, S1
## Structural benchmark

Step 0:
---
Source data:
 - 1,312 PDB chain sequences with structural tandem repeats (average TM-score 0.6 - 1) detected by CE-Symm

Files:
 - TPs_v5.fas.gz
 - FPs_v5.fas.gz

Step 1:
---
Generate lookup table to convert auth_seq_id to label_seq_id
Run TR detecting tools and parsing the results

Notes:
 - Analysis_structural_benchmark.txt

Scripts:
 - parse_HHrepID.2.r
 - parse_RADAR.2.r
 - parse_treks.4.r
 - parse_trust.5.r
	 
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


Validate tandem repeats identified by protein tandem repeat finding tools
---

Notes:
 - Analysis_Tally_validation.txt

Scripts:
 - parse_HHrepID2Tally.1.r
 - parse_trust.5.r
 - parse_RADAR.2.r
 - parse_treks.4.r


Get summary of source of PDB used in test data set
---

Notes:
 - Analysis_PDB_source.txt