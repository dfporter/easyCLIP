import HTSeq, collections, pandas, os, re, dill, importlib
import numpy as np
import sameRiver
import sameRiver.countsO
import sameRiver.annCounts
import sameRiver.scheme
import sameRiver.biotypeUtils
import sameRiver.readsPerGene
importlib.reload(sameRiver.readsPerGene)
import sameRiver.readsPerGene as rpg
importlib.reload(sameRiver.annCounts)


class positiveCounts(sameRiver.annCounts.annCounts):

    def __init__(self, metadata, xl_rate_fname: str=''):
        self.metadata = metadata
        
        self.metadata.data_folder = metadata.top_dir + '/data/'

        os.makedirs(self.metadata.data_folder, exist_ok=True)

        self.scheme = sameRiver.scheme.scheme(metadata.scheme_file)

        self.raw_reads_per_gene = rpg.rawReadsPerGene(
            metadata.ann_counts_file, scheme_filename=metadata.scheme_file)

        print('positiveCounts loaded rpg. proteins=', self.raw_reads_per_gene .proteins())
        # Throw out columns that aren't for positive proteins.
        self.raw_reads_per_gene.df = self.raw_reads_per_gene.df.loc[:,
            [x for x in self.raw_reads_per_gene.df.columns if \
                (self.raw_reads_per_gene.df[x].dtype.kind not in 'bifc') \
                or (self.scheme.gene_from_fname(x) in self.metadata.positive_proteins)
            ]]

        self.reads_per_million = rpg.readsPerMillion(
            self.raw_reads_per_gene,
            load_total_read_numbers=f"{self.metadata.data_folder}/total_read_numbers.txt")

        if xl_rate_fname:
            self.reads_per_protein = rpg.xlsPerProtein(
                self.reads_per_million, xl_rate_fname=xl_rate_fname)

    @staticmethod
    def load(fname='data/positives_countsO.dill'):
        print(f"Loading {fname}...")
        with open(fname, 'rb') as f:
            counts = dill.load(f)
        print("...Loaded.")
        return counts

    def save(
        self, fname=None, top_dir=None, write_txt=False, write_object=True):

        if fname is None:
            fname = f'{self.metadata.data_folder}/positives_countsO.dill'
        if top_dir is None:
            top_dir = self.metadata.data_folder

        os.makedirs(top_dir, exist_ok=True)

        if write_object:
            print(f"Saving to {fname}...")
            with open(fname, 'wb') as f:
                dill.dump(self, f)
            print("...Saved.")

        if write_txt and (top_dir is not None):
            print("Writing txt of positive counts...")
            raw = top_dir + '/positive_counts_raw.txt'
            to_per_read = top_dir + '/positive_counts_per_read.txt'
            to_per_protein = top_dir + '/positive_counts_per_protein.txt'

            self.raw_reads_per_gene.df.to_csv(raw, sep='\t')
            self.reads_per_million.df.to_csv(to_per_read, sep='\t')
            if hasattr(self, 'reads_per_protein'):
                self.reads_per_protein.df.to_csv(to_per_protein, sep='\t')
            print(f"...Wrote txt files (e.g. {to_per_read}).")

        if not os.path.exists(f'{self.metadata.data_folder}/logs/'):
            os.system(f'mkdir {self.metadata.data_folder}/logs')

        with open(f'{self.metadata.data_folder}/logs/datasets_included_in_positiveCounts_object.py', 'w'
            ) as f:
            f.write(
                'datasets_in_positiveCounts = ' + str(list(self.reads_per_protein.df.columns)) + '\n')

