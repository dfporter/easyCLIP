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

class negativeCounts():

    def __init__(self, metadata: object, xl_rate_fname: str=''):

        self.info = ''
        self.metadata = metadata

        self.scheme = sameRiver.scheme.scheme(metadata.scheme_file_with_random_proteins)

        self.metadata.data_folder = metadata.top_dir + '/data/'

        os.makedirs(self.metadata.data_folder, exist_ok=True)

        self.raw_reads_per_gene = rpg.rawReadsPerGene(
            metadata.ann_counts_file, scheme_filename=metadata.scheme_file_with_random_proteins)

        # Throw out columns that aren't for negative proteins.
        self.raw_reads_per_gene.df = self.raw_reads_per_gene.df.loc[:,
            [x for x in self.raw_reads_per_gene.df.columns if \
                (self.raw_reads_per_gene.df[x].dtype.kind not in 'bifc') \
                or (self.scheme.gene_from_fname(x) in self.metadata.random_proteins)
            ]]

        self.reads_per_million = rpg.readsPerMillion(
            self.raw_reads_per_gene,
            load_total_read_numbers=f"{self.metadata.data_folder}/total_read_numbers.txt")

        if xl_rate_fname:
            self.reads_per_protein = rpg.xlsPerProtein(
                self.reads_per_million, xl_rate_fname=xl_rate_fname)

        self.get_lowest_values_for_all_proteins()

    @staticmethod
    def load(fname='data/negatives_countsO.dill'):

        print(f"Loading {fname}...")
        with open(fname, 'rb') as f:
            counts = dill.load(f)
        print("...Loaded.")

        return counts

    def save(
        self, fname=None, top_dir=None,
        write_txt=False, write_object=True):

        if fname is None:
            fname = f'{self.metadata.data_folder}/negatives_countsO.dill'
        if top_dir is None:
            top_dir = self.metadata.data_folder

        os.makedirs(top_dir, exist_ok=True)
        
        if write_object:
            print("Saving to {}...".format(fname))
            with open(fname, 'wb') as f:
                dill.dump(self, f)
            print("...Saved.")

        if write_txt and (top_dir is not None):
            print("Writing txt of negative counts...")
            raw = top_dir + '/negative_counts_raw.txt'
            to_per_read = top_dir + '/negative_counts_per_read.txt'
            to_per_protein = top_dir + '/negative_counts_per_protein.txt'
                
            self.raw_reads_per_gene.df.to_csv(raw, sep='\t')
            self.reads_per_million.df.to_csv(to_per_read, sep='\t')
            if hasattr(self, 'reads_per_protein'):
                self.reads_per_protein.df.to_csv(to_per_protein, sep='\t')
            print("...Wrote txt files (e.g. {}).".format(to_per_read))

        if not os.path.exists(f'{self.metadata.data_folder}/logs/'):
            os.system(f'mkdir {self.metadata.data_folder}/logs')

        with open(f'{self.metadata.data_folder}/logs/datasets_included_in_negativeCounts_object.py', 'w'
            ) as f:
            f.write(
                'datasets_in_negativeCounts = ' + str(list(self.raw_reads_per_gene.df.columns)) + '\n')

    def generate_average_of_RNAs_with_low_counts_among_randoms(
        self, _counts_by_protein):

        negatives = []
        for gene_type, dict_of_proteins in _counts_by_protein.items():

            randoms = [
                np.mean(dict_of_proteins[prot]) for prot in self.metadata.random_proteins]
            if any(randoms) > 0:
                negatives.append(randoms)

        negatives = sorted(negatives, key=sum)

        print("{0} RNAs with nonzero values for some random control.".format(len(negatives)))

        if len(negatives) > 1000:
            frac = 0.1
        elif len(negatives) > 100:
            frac = 0.05
        else:
            frac = 0.1

        self.bottom = negatives[:int(len(negatives) * frac)]

        self.ave_bottom = np.mean(self.bottom, axis=0)

        print("{0} RNAs in the bottom {1} fraction of that. The average is: {2}".format(
            len(self.bottom), frac, self.ave_bottom))

        return self.ave_bottom


    def randoms_at_gene(
        self, gene, counts_by_protein=None, verbose=False, get_exon=True, 
        get_intron=False):
        
        randoms = self.metadata.random_proteins
        
        if get_exon and get_intron:
            raise ValueError("randoms_at_gene() cannot only return an exon or intron, not both.")
        
        if len(gene.split('::')) > 1:
            gene_type = gene
        elif get_exon:
            gene_type = gene + '::exon'
        elif get_intron:
            gene_type = gene + '::intron'

        if gene_type not in counts_by_protein:
            return [[np.nan] * len(randoms)]

        return [counts_by_protein[gene_type][prot] for prot in randoms \
                    if prot in counts_by_protein[gene_type]]

    def randoms_at_genes(
        self, genes='all', counts_by_protein=None):
        
        if genes == 'all':
            genes = counts_by_protein.keys()
            genes = set([x.split('::')[0] for x in genes])
            genes -= set(['_ambigious', '_no_feature'])
            genes = list(genes)
        
        assert type(genes) == type([])
        
        counts_from_randoms = []
        for gene in genes:
            counts_from_randoms.append(
                self.randoms_at_gene(gene, counts_by_protein=counts_by_protein))
        
        return counts_from_randoms

    def get_lowest_values_for_all_proteins(self, maximum_per_read_floor=1E9):

#        print("""get_lowest_values_for_all_proteins():
#            for negatives, if 1E6/total_read_number > {}, the lowest possible positive value
#            for per (million) reads (replaces zero) is set to {}.""".format(
#                maximum_per_read_floor, maximum_per_read_floor))

        self.lowest_positive_vals = {
            'raw': collections.defaultdict(int),
            'per_read': collections.defaultdict(int),
            'per_protein': collections.defaultdict(int),      
        }

        def set_lowest(rpg_object, dict_to_set):

            # Use averages for a protein.
            for protein in rpg_object.proteins():
                for col in rpg_object.columns_for_a_protein(protein):
                    arr = rpg_object.df[col].to_numpy()
                    a = arr[np.where(arr > 0)]
                    if len(a) and (protein not in dict_to_set):
                        dict_to_set[protein] = np.min(a)
                    elif len(a):
                        dict_to_set[protein] = np.mean([np.min(a), dict_to_set[protein]])
                if (protein not in dict_to_set) or dict_to_set[protein]==0:
                    dict_to_set[protein] = 1E-9

            # Get values for individual columns as well.
            for col in rpg_object.numeric_columns(rpg_object.df):
                arr = rpg_object.df[col].to_numpy()
                a = arr[np.where(arr > 0)]
                if len(a):
                    dict_to_set[col] = np.min(a)
                else:
                    dict_to_set[col] = 1E-9

        set_lowest(self.raw_reads_per_gene, self.lowest_positive_vals['raw'])
        set_lowest(self.reads_per_million, self.lowest_positive_vals['per_read'])
        if hasattr(self, 'reads_per_protein'):
            set_lowest(self.reads_per_protein, self.lowest_positive_vals['per_protein'])


