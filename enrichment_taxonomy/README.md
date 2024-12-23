# Enrichment analysis for taxonomy
To generate figure 4A

Notes:
 - Analysis_enrichment_taxonomy.txt


Step 0: Prepare source data
---
Download aa fasta files from NCBI on 2022-11-18

File paths: (sequence files are too large, can be provided upon request)
 - 4028 prokaryote: ./faa_20221118/GCF_?????????.faa.gz
 - 590 eukaryote: ./faa_20221118/GCA_?????????.faa.gz 

	 
Step 1: Make species tree for genomes
---
Collect SSU RNA sequences from each genome using DECIPHER::FindNonCoding() \
Select unit with most SSU patterns matched \
Align SSU RNA sequences \
Reconstruct species tree using DECIPHER::TreeLine()

File paths:
 - pro_4077.rds
 - euk_1775.rds
 - species_4618.tree


Step 2: Run DetectRepeats
---
Run DetectRepeats on 4618 genomes' proteome

Scripts:
 - DECIPHER_DetectRepeats.18.r


Step 3: Compute phylogenetic signal Blomberg's K
---
Compute Blomberg's K for proteome TR fraction using phytools::phylosig() \

File paths:
 - all4618_ML_K_07012024.rds


Step 4: Plot species tree with TR fraction
---
Assign taxonomy using DECIPHER::IdTaxa() \
Plot tree with bar