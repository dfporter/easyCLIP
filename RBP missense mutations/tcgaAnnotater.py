import pandas, os, sys, re, time, collections, Bio, pprint, random, glob

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO

import proteinLengths

from typing import List, Mapping, Union

class tcgaAnnotater():

    def __init__(self, tcgaData, verbose=True):
        self.verbose = verbose
        self.load(tcgaData)

    def load(self, tcgaData):

        print("By specific mutations:")
        self.sort_by_mutations_and_print(tcgaData.n_by_mutation, verbose=self.verbose)
        print("By all mutations per gene:")
        self.sort_by_mutations_and_print(tcgaData.n_by_gene, verbose=self.verbose)
        
        print("Creating dataframes.")
        self.df_by_mutation = self.dataframe_sorted_by_mutations_per_length(
            tcgaData.n_by_mutation)
        self.df_by_gene = self.dataframe_sorted_by_mutations_per_length(
            tcgaData.n_by_gene)

        return (self.df_by_mutation, self.df_by_gene, 
            tcgaData.n_by_mutation, tcgaData.n_by_gene, 
            tcgaData.by_mutation, tcgaData.by_gene)

    @staticmethod
    def sort_by_mutations_and_print(by_mutation, verbose=True):
        """Just prints. No state changes."""
        mutations = sorted(by_mutation.keys(), key=lambda x: by_mutation[x])

        if verbose:
            print([f"{mut}, num patients={by_mutation[mut]}" for mut in mutations[-10:]])

        snp = [mut for mut in mutations if '*' not in mut[1]]
        print("Top mutations by frequency, no nonsense mutations:")
        print([f"{mut}, num patients={by_mutation[mut]}" for mut in snp[-10:]])
            #print(mut, by_mutation[mut])

        if verbose:
            print("Gene list for copy-pasting:")
            if type(snp[0]) != type(''):
                print('\n'.join([x[0] for x in snp[-10:]]))

    def get_protein_length(self, gene):
        if self.protein_lengths[gene] == 0:
            print('No length for ', gene, end=', ')
            length = 1E6
        else:
            length = self.protein_lengths[gene]    
        return length

    @staticmethod
    def top_snp_df(by_mutation):
        mutations = sorted(by_mutation.keys(), key=lambda x: by_mutation[x])    
        rows = []
        for (gene, snp), count in by_mutation.items():
            if '*' in snp:
                continue
            if count < 2:
                continue
            rows.append({
                'Gene': gene,
                'Variant': snp,
                'Count': count,
            })
        df = pandas.DataFrame(rows)
        df.sort_values(by='Count', ascending=False, inplace=True)
        return df

    def dataframe_sorted_by_mutations_per_length(self, by_mutation):

        if not hasattr(self, 'protein_length'):
            pl = proteinLengths.proteinLengths()
            self.protein_lengths = pl.load()

        by_gene = collections.defaultdict(int)
        rows = []
        
        try:
            for (gene, alt), muts in by_mutation.items():
                by_gene[gene] += muts
                
            sorted_by_n_per_gene = sorted(by_gene.keys(), key=lambda x: by_gene[x])

            for gene in by_gene.keys():
                length = self.get_protein_length(gene)
                rows.append({
                    'Gene': gene, 'N mutations': by_gene[gene], '# AA': length,
                    'Mut/length': by_gene[gene]/length, 
                })        
                
        except:  # This is by_gene already.
            for gene in by_mutation.keys():
                length = self.get_protein_length(gene)
                rows.append({
                    'Gene': gene, 'N mutations': by_mutation[gene], '# AA': length,
                    'Mut/length': by_mutation[gene]/length, 
                })     

        print('\nOutput from tcgaAnnotater:')
        df = pandas.DataFrame(rows)
        df.sort_values(by='Mut/length', inplace=True, ascending=False)
        print(df.head(5))
        print('\n'.join([x for x in df['Gene']][:20]))
        return df



