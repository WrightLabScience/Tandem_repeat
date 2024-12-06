# Tandem repeat units' percent identities distribution for each gene group
To generate figure 6A, 6B, S1

Step 0: Prepare source data
---
Prokaryotic and eukaryotic KEGG orthologous gene groups (Kgroups)

Input:
 - 13485 eukaryotic Kgroups
 - 9734 prokaryotic Kgroups

Notes:
 - Analysis_TR_unit_distance.txt


Step 1: Run DetectRepeats
---
Find tandem repeat in Kgroups' amino acid (aa) sequences using DetectRepeats()

Notes:
 - Analysis_TR_unit_distance.txt

Script:
 - DECIPHER_DetectRepeats.10.r

	 
Step 2: Calculate unit distances of the main tandem repeat in each Kgroup
---
Make multiple sequence alignment (MSA) of each Kgroup \
Filter tandem repeat results: Keep longer (unit length > 10 aa) and higher confidence (score > 100) repeats, with conserved unit length \
Find the main tandem repeat that present in > 50% of the genes in each Kgroup \

Notes:
 - Analysis_TR_unit_distance.txt

Script:
 - TR_dist.15.r


Step 3: Plot
---
Plot heatmap and eCDF

Notes:
 - Analysis_TR_unit_distance.txt