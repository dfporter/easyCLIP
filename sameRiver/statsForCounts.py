import pandas, collections, importlib, random, re, bisect, os, sys, dill

import sklearn
from sklearn.neighbors.kde import KernelDensity
import numpy as np
import scipy
import scipy.stats as scs

import statsmodels.api as sm

import sameRiver
import sameRiver.countsO, sameRiver.positiveCounts, sameRiver.negativeCounts

from typing import Union, List, Mapping

importlib.reload(sameRiver.countsO)


class statsForCounts():
    
    def __init__(
        self,
        negatives: Union[None, object, sameRiver.negativeCounts.negativeCounts]=None,
        positives: Union[None, object, sameRiver.positiveCounts.positiveCounts]=None,
        data_dir: str='./data/'):

        self.negatives = negatives
        self.positives = positives
        self.info = ''
        self.pvals = {}
        self.data_dir = data_dir

    def save(self, fname=''):

        if fname == '':
            fname = f"{self.positives.metadata.data_folder}/stats.dill"

        os.makedirs(os.path.dirname(fname), exist_ok=True)
        
        print(f"Saving to {fname}...")
        with open(fname, 'wb') as f:
            dill.dump(self, f)
        print("...Saved.")

    @staticmethod
    def load(fname=''):

        if fname == '':
            fname = f"{self.positives.metadata.data_folder}/stats.dill"

        print(f"Loading {fname}...")
        with open(fname, 'rb') as f:
            counts = dill.load(f)
        print("...Loaded.")
        return counts
    
    def counts_for_a_protein(
            self, protein='PCBP1', counts_by_protein=None, get_exon=True, get_intron=False,
            as_dict=True, as_list=False):

        with_protein = collections.defaultdict(dict)
        suffixes = set()
        if get_exon:
            suffixes.add('exon')
        if get_intron:
            suffixes.add('intron')
            
        for gene_type, dict_of_proteins in counts_by_protein.items():
            if gene_type.split('::')[-1] not in suffixes:
                continue
            
            with_protein[gene_type] = dict_of_proteins.get(protein, [np.nan])

        if as_list:
            return list(with_protein.values())
        
        if as_dict:
            return with_protein
       
    def example_random(self, _counts_by_protein=None, n_examples=1, verbose=True):
        if _counts_by_protein is None:
            _counts_by_protein = self.counts_by_protein
            
        n_examples = int(n_examples)
        genes = list(_counts_by_protein.keys())

        randoms = self.negatives.random_proteins

        for n in range(0, n_examples):
            gene = ''

            not_intron = False
            while not not_intron:
                gene = random.choice(genes)
                if not re.search('intron', gene):
                    not_intron = True

            if verbose:
                print("Gene: ", gene)
                print("Values: ", [(_, _counts_by_protein[gene][_]) for _ in randoms if _ in _counts_by_protein[gene]])

        return [_counts_by_protein[gene].get(_, 0) for _ in randoms]

    @staticmethod
    def scale_all_numbers(_counts_by_protein, scale_all_numbers_by):
        print("Scaling all count data points by {} before calculating p values...".format(
            scale_all_numbers_by))
        counts_by_protein = collections.defaultdict(dict)

        numbers_used = set()
        for rna in _counts_by_protein:
            if random.randint(0, 1E4) == 1:
                print(rna)
            for protein, val in _counts_by_protein[rna].items():
                ival = val[:]
                counts_by_protein[rna][protein] = [x * scale_all_numbers_by for x in val]
                #self.scaled_counts_by_protein = counts_by_protein
                if np.sum(ival) == 0 and np.sum([x * scale_all_numbers_by for x in val]) != 0:
                    print("Zero converted to nonzero somehow.")
        return counts_by_protein

    @staticmethod
    def flip_dict_order(dict_of_dicts):
        by_prot = collections.defaultdict(dict)

        for gene_type, dict_by_prot in dict_of_dicts.items():
            
            if dict_by_prot is None:
                continue
                
            for prot, pval in dict_by_prot.items():
                by_prot[prot][gene_type] = pval

        return by_prot

    def summarize_pvalues(self, pvals):
        by_prot = self.flip_dict_order(pvals)

        for prot in by_prot:
            print(prot)
            plt.clf(); plt.close()
            plt.hist(list(by_prot[prot].values()), bins=100)
            plt.show()

    def targets(self, protein, cutoff=0.01, which='per_read'):
        pvals = self.pvals_dfs[which]

        targets = pvals[pvals[protein]<cutoff].copy()
        targets = set(targets.index)
        return targets

    def write_targets(
        self, which='per_read', outfname='default',
        p_value_cutoff=None, apply_bh_adjust=False):

        self.write_pvals_single_file(which=which, outfname=outfname)

        if outfname == 'default':
            outfname = f'{self.positives.metadata.top_dir}/tables/pvals_{which}.xlsx'

        os.makedirs(os.path.dirname(outfname), exist_ok=True)

        writer = pandas.ExcelWriter(
            outfname + '_separate_sheets.xlsx', engine='xlsxwriter')
        workbook = writer.book

        # Add some cell formats.
        #format1 = workbook.add_format({'num_format': '#,##0.00'})
        #format2 = workbook.add_format({'num_format': '0%'})

        # Make by_prot = {protein => {RNA => p value}}
        pvals = self.pvals[which]
        by_prot = self.flip_dict_order(pvals)
        
        out_dir = os.path.dirname(outfname)
        if not os.path.exists(out_dir + '/sub/'):
            os.system('mkdir ' + out_dir + '/sub/')

        to_biotype = dict(zip(
            self.negatives.raw_counts_df.gene_name,
            self.negatives.raw_counts_df['Gene type']))

        for prot, rnas in by_prot.items():
            # Get the xlsxwriter workbook and worksheet objects.
            
            sheet_name = re.sub(':', ' ', prot)

            print('****', prot)
            df = pandas.DataFrame.from_dict(rnas, orient='index')
            df.columns = ['P value']
            df['gene_name'] = df.index

            # FYI:
            df_sig = df[df['P value']<0.01].copy()
            gene_names = set([x.split('::')[0] for x in df_sig.index])
            print("RNAs with P<{}: {}".format(0.01, len(gene_names)))
            # End FYI.

            neg = self.negatives.raw_counts_df.copy()
            # Remove non-numeric.
            neg = neg.loc[:, [col.dtype.kind not in 'bifc' for col in neg.columns]]
            neg.loc['Protein', :] = [scheme_.gene_from_fname(x) for x in neg.columns]

            if prot in self.negatives.metadata.random_proteins:
                df['Raw counts'] = [
                    np.mean(self.negatives.raw_counts_by_protein.get(gene, [0])[prot]) for gene in df.index]
                df['Reads per million'] = [
                    np.mean(self.negatives.per_million_by_protein.get(gene, [0])[prot]) for gene in df.index]
                df['Gene type'] = [
                    to_biotype[gene] for gene in df.index]
            else:
                df['Raw counts'] = [
                    np.mean(self.positives.raw_counts_by_protein.get(gene, [0])[prot]) for gene in df.index]
                df['Reads per million'] = [
                    np.mean(self.positives.per_million_by_protein.get(gene, [0])[prot]) for gene in df.index]
                df['Gene type'] = [
                    to_biotype[gene] for gene in df.index]

            df.sort_values(by=['P value'], inplace=True, ascending=True)

            df.to_csv(out_dir + '/sub/{}.txt'.format(sheet_name), sep='\t')

            if p_value_cutoff is not None:
                df = df[df['P value']<=p_value_cutoff]
            
            df.to_excel(writer, sheet_name=sheet_name)
            writer.sheets[sheet_name].set_column('A:A', 20)
            
        writer.close()

    def write_pvals_single_file(
        self, which='per_read', outfname='default'):

        if outfname == 'default':
            outfname = f'{self.positives.metadata.top_dir}/tables/pvals_{which}.xlsx'

        # Make by_prot = {protein => {RNA => p value}}
        pvals = self.pvals[which]
        
        os.makedirs(os.path.dirname(outfname), exist_ok=True)

        out_dir = os.path.dirname(outfname)

        try:
            if outfname.split('.')[-1] in ['csv', 'txt']:
                self.pvals_dfs[which].to_csv(outfname, sep='\t')
            else:
                self.pvals_dfs[which].to_excel(outfname)

        except:
            print("Tried to write a table of pvalues to {} but failed. Trying {}.".format(
                outfname, out_dir + '/pvalues.csv'))

            self.pvals_dfs[which].to_csv(
                        out_dir + f'/pvalues_{which}.csv', sep='\t')

