# Prokaryotic gene K00689 tandem repeat analysis
To generate figure 6C, 6D, S4

Step 0: Prepare source data
---
Prokaryotic K00689 gene sequences

File paths:
 - 24 gene sequences: prok_K00689_seqs.fasta
   - 23 K00689 genes with tandem repeat (TR) + 1 outgroup

Notes:
 - Analysis_prok_K00689_TR.txt


Step 1: Run DetectRepeats
---
Find tandem repeat in K00689 amino acid (aa) sequences using DetectRepeats()

Notes:
 - Analysis_prok_K00689_TR.txt

	 
Step 2: Make gene tree for K00689 genes
---
Make multiple sequence alignment (MSA) of K00689 amino acid sequences \
Remove tandem repeat region from MSA \
Make gene tree in Figure 6D

File paths:
 - MSA of K00689 aa seqs: prokK00689_faa_msa.3.fasta
 - MSA of K00689 aa seqs after removing TR region: prokK00689_faa_msa_rmTR.fasta
 - K00689 gene tree: prokK00689_faa_msa_rmTR_ML.tree

Notes:
 - Analysis_prok_K00689_TR.txt


Step 3: Make species tree for genomes included in this analysis
---
Collect 16S, rpoB, rpoC gene sequences \
Make multiple sequence alignment \
Make species tree in Figure S3B

File paths:
 - input 16S DNA sequences: proK00689_16S.fna
 - input rpoB DNA sequences: proK00689_rpoB.fna
 - input rpoC DNA sequences: proK00689_rpoC.fna
 - concat MSA of 16S, rpoB, rpoC: proK00689_16SrpoBCnt_MSA.3.fna
 - species tree: proK00689_16SrpoBCnt_ML.3.tree

Notes: 
 - Analysis_prok_K00689_species_tree.txt


Step 4: Make TR unit tree and cluster units
---
Define reference units and generate a fasta file that contains non-partial TR units \
Make MSA of TR units \
Make TR unit tree and distance matrix in Figure S2 \
Cluster TR units (colored clades in Figure S2) \
An example of K00689 TR unit MSA and distance matrix for Figure 6C

File paths:
 - TR unit sequences: prok_K00689_units_rmp.fasta
 - TR unit tree: prok_K00689_TR_units_ML.tree
 - TR unit MSA: prok_K00689_TR_units_MSA.fasta

Notes: 
 - Analysis_prok_K00689_units.txt