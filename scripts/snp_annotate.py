# Imports
from Bio.Data import CodonTable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Define some functions
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
def get_aa(codon):
    return standard_table.forward_table[codon]
def get_complement(codon):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement[codon[0]]+complement[codon[1]]+complement[codon[2]]


# Read in files
genes = pd.read_csv('genes.txt', sep='\t')
snps = pd.read_csv('SNPs.txt', sep='\t')
with open('genome.fasta') as genome:
    genome_line = genome.read().split('\n')
    contig = {'Contig_1': genome_line[1], 'Contig_2': genome_line[3], 'Contig_3': genome_line[5]}

    
# Add a column that represents directionality of genes
# Then correct 'Start' & 'End' to represent lower & higher index.
genes['if_forward'] = genes['Start'] < genes['End']
genes['start_adj'] = [min(i[1]['Start'], i[1]['End']) for i in genes.iterrows()]
genes['end_adj'] = [max(i[1]['Start'], i[1]['End']) for i in genes.iterrows()]


# Add some empty columns
snps = snps.reindex(columns=['Contig', 'Position', 'SNP', 'gene_index', 'snp_in_gene',
                             'ref_codon', 'mod_codon', 'ref_aa', 'mod_aa', 'kind'])


# Add a foreign key that references index of gene the SNP lies in
for contig_no in snps.Contig.unique():
    snps_sub = snps[snps.Contig == contig_no]
    genes_sub = genes[genes.Contig == contig_no]
    for i in snps_sub.index:
        for j in genes_sub.index:
            snp_position = snps.loc[i, 'Position']
            condition = snp_position in range(
                (genes.loc[j, 'start_adj']-1), (genes.loc[j, 'end_adj']))
            if condition:
                snps.loc[i, 'gene_index'] = j

                
# Get row indexes of Intragenic SNPs             
sub_index = snps[pd.notnull(snps.gene_index)].index


# Populate blank columns with additional information for Intragenic SNPs
for i in sub_index:
    gene_index = snps.loc[i, 'gene_index']
    gene_start = genes.loc[gene_index, 'Start']
    gene_forward = genes.loc[gene_index, 'if_forward']
    snp_position = snps.loc[i, 'Position']
    snp_in_gene = abs(gene_start-snp_position)
    contig_no = snps.loc[i, 'Contig']
    contig_seq = contig[contig_no]
    snp = snps.loc[i, 'SNP']

    
    # Add positional index of snp in corresponding gene
    snps.loc[i, 'snp_in_gene'] = snp_in_gene

    
    # Add original amino acid codon where the SNP lies
    if gene_forward:
        codon_start = (snp_in_gene//3*3)+gene_start
        snps.loc[i, 'ref_codon'] = contig_seq[codon_start-1:codon_start+2]
    else:
        codon_start = gene_start-(snp_in_gene//3*3)
        snps.loc[i, 'ref_codon'] = contig_seq[codon_start-1] + \
            contig_seq[codon_start-2] + contig_seq[codon_start-3]

        
    # Add modified codon
    ref_codon = snps.loc[i, 'ref_codon']
    if snp_in_gene % 3 == 0:
        snps.loc[i, 'mod_codon'] = snp+ref_codon[1:]
    elif snp_in_gene % 3 == 1:
        snps.loc[i, 'mod_codon'] = ref_codon[0] + snp + ref_codon[2]
    elif snp_in_gene % 3 == 2:
        snps.loc[i, 'mod_codon'] = ref_codon[:2] + snp

        
    # Add original & modified Amino Acids
    mod_codon = snps.loc[i, 'mod_codon']
    if gene_forward:
        snps.loc[i, 'ref_aa'] = get_aa(ref_codon)
        snps.loc[i, 'mod_aa'] = get_aa(mod_codon)
    else:
        snps.loc[i, 'ref_aa'] = get_aa(get_complement(ref_codon))
        snps.loc[i, 'mod_aa'] = get_aa(get_complement(ref_codon))

        
    # Annotate Synonymous & Non-synomymous mutations
    if snps.loc[i, 'ref_aa'] == snps.loc[i, 'mod_aa']:
        snps.loc[i, 'kind'] = 'Synonymous'
    else:
        snps.loc[i, 'kind'] = 'Non-synonymous'

        
# Populate kind column with information for Intergenic SNPs
snps['kind'].fillna('Intergenic', inplace=True)


#Plot bar chart. 
x=[1,2,3,4,5,6,7,8,9]
y=snps[['Contig','kind','Position']].groupby(['Contig','kind']).count().values
one = mpatches.Patch(color = "green", label = "1")
two = mpatches.Patch(color = "red", label = "2")
three = mpatches.Patch(color = "blue", label = "3")

plt.bar(x, y, color = ['green','red','blue','green','red','blue','green','red','blue'], edgecolor="black")
plt.xticks([2,5,8], ['Non-synonymous', 'Synonymous', 'Intergenic'], fontsize = 7)
plt.yticks(fontsize = 7)
plt.ylabel("Number of SNPs", fontname = "Arial", fontsize = 7)
plt.title("Mutation profile", fontname = "Arial", fontsize = 7)
plt.legend(handles = [one, two, three], title = "Contig", fontsize='xx-small')
plt.savefig("outputplot", format ='pdf')
