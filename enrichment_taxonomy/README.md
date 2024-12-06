# Enrichment analysis for taxonomy
To generate figure 4A

Step 0: Prepare source data
---
4618 aa fasta files downloaded from NCBI on 2022-11-18

File paths:
 - 4028 prokaryote: ./faa_20221118/GCF_?????????.faa.gz
 - 590 eukaryote: ./faa_20221118/GCA_?????????.faa.gz 

Notes:


Step 1: Run DetectRepeats
---
Run DetectRepeats on 4801 genomes' proteome

File paths:

Notes: 
 - DetectRepeats_proteome.txt (TR02.96.txt)
   - input: $genomeID.faa.gz
   - output: $genomeID.tr

Scripts:
 - DECIPHER_DetectRepeats.18.r

	 
Step 2: Make species tree for 4801 genomes
---
Collect SSU RNA sequences from each genome using DECIPHER::FinNonCoding() \
Select unit with most SSU patterns matched \
Reconstruct species tree using DECIPHER::TreeLine()



File paths:


Notes: 
 - make_4801_sp_tree.txt (TR02.99.txt)


Scripts:
 - XXX


Step 3: Compute phylogenetic signal Blomberg's K
---
Compute Blomberg's K for proteome TR fraction using phytools::phylosig() \
Plot species tree with TR fraction

Notes: 
 - phylogsig_plot.txt

Scripts:
 - XXX
