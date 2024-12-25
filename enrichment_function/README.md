# Enrichment analysis for functional groups
To generate figure 4B

Notes:
 - Analysis_enrichment_function.txt


Step 0: Prepare source data
---
Protein sequences used in figure 4A (aa fasta files downloaded from NCBI on 2022-11-18)

File paths: (sequence files are too large, can be provided upon request)
 - 4028 prokaryote: ./faa_20221118/GCF_?????????.faa.gz
 - 590 eukaryote: ./faa_20221118/GCA_?????????.faa.gz 

	 
Step 1: Classify protein function
---
Identify KEGG orthologous group (KO group) for each protein using DECIPHER::IdTaxa()

File paths: (sequence files are too large, can be provided upon request)
 - 10,227 prokaryotic KO groups: ./prok_Kgroup_faa/K?????.faa.gz
 - 13,671 eukaryotic KO groups: ./euk_Kgroup_faa/K?????.faa.gz 

Scripts:
 - DECIPHER_IdTaxa.7.r
 - collect_Kgroup_seq.1.r
 - write_Kgroup_faa.r


Step 2: Run DetectRepeats
---
Run DetectRepeats on protein sequences of each KO group

Scripts:
 - DECIPHER_DetectRepeats.19.r


Step 3: Calculate values and plot
---
x-axis: Fraction of genomes have the Kgroup (perc_genomes) \
y-axis: Fraction of Kgroup have tandem repeat (perc_TR_genes) \
dot size: DetectRepeats average score (Scores_highest_ave) \ 
shade: Average sequence identity of genes in Kgroup (PID)

Scripts:
 - calculate_sample_PID.1.r

File paths:
 - table_function_prok_TR.tsv
 - table_function_euk_TR.tsv