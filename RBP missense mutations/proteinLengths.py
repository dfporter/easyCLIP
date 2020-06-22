import pandas, os, sys, re, time, collections, Bio, pprint, random, glob

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO

class proteinLengths():

    def load(
        self, uniprot_fname='uniprot_id_to_gene_name.tab',
        aa_length_fname='./uniprot-proteome%3AUP000005640.fasta'):

        # Load IDs.
        df = pandas.read_csv(uniprot_fname, sep='\t', index_col=False)

        df.columns = ['Uniprot', 'Gene']
        uniprot_to_gene_name = dict(zip([x.upper() for x in df.Uniprot], [x.upper() for x in df.Gene]))

        # Get protein lengths.
        protein_length = collections.defaultdict(int)
        li = ''

        for seq_record in SeqIO.parse(aa_length_fname, "fasta"):
            m = re.search(' GN=(\w+) ', seq_record.description)
            gene = seq_record.id.split('|')[-1].split('_')[0].upper()

            if m is not None:
                gene = m.group(1)

            if 'Fragment' in seq_record.id:
                continue
            if gene in uniprot_to_gene_name:
                protein_length[uniprot_to_gene_name[gene]] = max([
                    protein_length[uniprot_to_gene_name[gene]], len(seq_record.seq)])
                protein_length[gene] = max([protein_length[uniprot_to_gene_name[gene]], len(seq_record.seq)])
            else:
                protein_length[gene] = max([protein_length[gene], len(seq_record.seq)])

            li += gene + '\n'
            
        with open('uniprot_gene_ids.txt', 'w') as f:
            f.write(li)

        return protein_length