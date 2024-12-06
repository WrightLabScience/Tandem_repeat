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
Multiple sequence alignment \
DECIPHER::AlignSeqs() for aa \
DECIPHER::AlignTranslate() for nt 

File paths: \
Notes: 
 - prepare_storeData.txt (TR02.94.txt)
   - input: Kgroup aa and nt fasta
   - output: t_aa.msa, mt_nt.msa, pairs

Scripts:
 - consist_storeData.4.r
	 
Step 2:
---
Run TR detecting tools and parsing the results


File paths:


Notes: 
 - run_parse_DetectRepeats.txt (TR02.99.txt)
 - run_parse_HHrepID.txt (TR02.95.txt)
 - run_parse_RADAR.txt (TR02.95.txt)
 - run_parse_TRUST.txt (TR02.95.txt)
 - run_parse_TREKS.txt (TR02.95.txt)
 - run_parse_XSTREAM.txt (TR02.95.txt)
 - run_parse_TRF.txt (TR02.95.txt)
 - run_parse_mrepstxt (TR02.95.txt)

Scripts:
 - parse_HHrepID.2.r


Step 3:
---
Compute consistency


Scripts:
 - consist_all.3.r
