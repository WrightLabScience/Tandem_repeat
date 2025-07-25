# To generate figure 3C
## Ortholog benchmark

Step 0:
---
Source data:
 - KEGG orthology group (Kgroup)
 - aa and nt fasta files
 - 21211 Kgroups (K?????)
 - Random subsample to include at most 100 sequences per Kgroup

File paths: (sequence files are too large, can be provided upon request)
 - ./nt/K?????.fasta
 - ./aa/K?????.fasta 
 - Kgroup_21211.txt

Step 1:
---
Multiple sequence alignment
 - DECIPHER::AlignSeqs() for aa
 - DECIPHER::AlignTranslate() for nt 

Store sequence pairs in 50% ~ 60% PID range

Select Kgroup with >= 500 n pairs

File paths:
 - Kgroup_ge500.txt
 - pairs/ (can be provided upon request)

Notes: 
 - Analysis_orthology_benchmark.txt
   - input: Kgroup aa and nt fasta
   - output: msa, pairs

Scripts:
 - DECIPHER_MSA.5.r
 - consist_storeData.4.r
	 
Step 2:
---
Run TR detecting tools and parsing the results

Notes: 
 - Analysis_orthology_benchmark.txt

Scripts:
 - DECIPHER_DetectRepeats.20.r
 - parse_HHrepID.2.r
 - parse_RADAR.2.r
 - parse_trust.6.r
 - parse_treks.4.r
 - parse_TRF.1.py
 - parse_MREPS.1.py


Step 3:
---
Compute consistency and plot

Notes: 
 - Analysis_orthology_benchmark.txt

Scripts:
 - consist_all.3.r


Venn Diagrams
---
Compute consistency of tandem repeat identification across DNA and AA sequences

Notes: 
 - Analysis_orthology_benchmark.txt