# Enrichment analysis for genomic region
To generate figure 5

Notes:
 - Analysis_enrichment_genomic_region.txt


Step 0: Prepare source data
---
Whole genome sequence (*.fna.gz) and gene annotation (*.gff.gz) files downloaded from NCBI on 2024-07-02 \
121 eukaryotic genomes \
5374 prokaryotic genomes

File paths: (files are too large, can be provided upon request)
 - list of genomes to download: 
   - assembly_summary_prok_1.rds
   - assembly_summary_euk_1.rds
 - whole genome sequence: ./source_data_07022024/genomic.fna/*.fna.gz
 - gene annotation: ./source_data_07022024/genomic.gff/*.gff.gz

Scripts:
 - fetch_source_data.1.r


Step 1: Run DetectRepeats
---
Run DetectRepeats on whole genome sequence

Scripts:
 - DetectRepeats_gfna.3.r



Step 2: Calculate TR frequency
---
Calculate TR frequency for each genome \
Combine TR frequencies of all prokaryotic or eukaryotic genomes

Scripts:
 - CDSregion_TR.9.r

File paths: 
 - TR_freq.RData



Step 3: TR enrichment analysis
---
Plot TR frequency, median, genome-wide average \
Plot delta of TR frequency median between genic and intergenic \
region around start codon: distance to M \
region around stop codon: distance to * \
Find position(distance) with TR enrichment \
Adjust P-values for multiple comparisons

Scripts:
 - bgAvgTR.2.r

File paths: 
 - KStest_TR.RData
 - bgAvgTR_prok.tsv
 - bgAvgTR_euk.tsv