# Microbial-SNP-Annotation-Python

## Goal
To estimate the number of non-synonymous, synonymous and intergenic SNPs present in a sampled population.
  
    
## Data 
  - **SNP file**: 
    - A list of single nucleotide polymorphisms (SNPs) from a bacterial population. 
    - The information in this file was obtained by mapping the raw sequencing reads of individual bacterial isolates to a complete reference genome of the same species. The **contig of origin, position in the contig** and **nucleotide change** are displayed for every SNP.
  - **Gene file**: 
    - Contains the **geneID, gene description, contig of origin** and **genome coordinates (gene start and gene end)** for every gene encoded in the reference genome.
  - **Genome file**: 
    - a .fasta file containing the **complete genome** of the reference sample.

<img src=".\data\git_snp.png" alt="Drawing" style="width: 60px;"/>

## Key steps
1. Read in :
    1. .txt files into a pandas dataframe. 
    2. genome.fasta into a dictionary where 'Contig_#' is key & sequence is corresponding value. 
2. Compare SNP 'Position' with gene 'Start' & 'End'. We will thus identify **Intragenic & Intergenic SNPs**. For Intragenic SNPs - assign a foreign key that corresponds to row index of gene it lies in. 
3. Get **row indexes of Intragenic SNPs**. Loop through these to assign: 
    1. **Positional index of snp** in corresponding gene. 
    2. **Original amino acid** codon where the SNP lies. 
    3. **Modified codon**. 
    4. Original & Modified **Amino Acids** (by referencing BioPython module's 'Standard Codon Table').
    5. Finally, compare original & modified Amino Acids to annotate **Synonymous & Non-synomymous mutations**. 
4. Label the rest as Intergenic mutations. 
5. Generate **bar chart** to compare the number of non-synonymous, synonymous and intergenic SNPs present in a sampled population.

## Result
 - All three contigs have Synonymous mutations. 
 - **Future work** - Explore gene functions associated with the Synonymous SNPs & also explore if any Intergenic SNPs lie in genomic regions important for microbial regulatory functions. 
<img src=".\data\outputplot.png" alt="Drawing" style="width: 60px;"/>
