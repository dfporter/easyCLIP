import HTSeq, collections, pandas, os, re, dill, importlib, random, glob
import numpy as np

import sameRiver
import sameRiver.scheme
import sameRiver.biotypeUtils
import sameRiver.countsO
importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.biotypeUtils)
importlib.reload(sameRiver.countsO)

class annCounts(sameRiver.countsO.countsO):
    """For loading an ann_counts.txt file of counts data output by
    sameRiver.countO.countsO.edit(). These files have biotype data
    and have been purged and filled in.

    Normalization to per read or per protein are handled by this
    object.
    """
    @staticmethod
    def numeric_columns(df):
        return [x for x in df.columns if (
            df[x].dtype.kind in 'bifc')]

    def normalize_ann_counts(
        self, xl_rate_fname='percentCrosslinked.xlsx',
        include_percent_xl=True,
        use_bed_dir=None, load_total_read_numbers=False):
        """ Set self.counts_per_million_df and self.counts_per_protein_df from
        self.raw_counts_df and a folder of bed files, the latter used to determine total
        read counts. After counts_per_million_df is set by normalizing raw_counts_df by
        total read number, counts_per_protein is set from counts_per_million_df using the
        xl rate filename (an excel file).
        """

        print('annCounts.normalize_ann_counts():')
        print(self.raw_counts_df.iloc[:4, :4])
        print(self.raw_counts_df.columns)

        # Set self.counts_per_million_df dataframe.
        self.normalize_counts_to_per_read(
            use_bed_dir=use_bed_dir, load_total_read_numbers=load_total_read_numbers)
        if include_percent_xl:
            # Sets counts_per_protein_df from counts_per_million_df and xl_rate_fname.
            self.normalize_counts_to_per_protein(xl_rate_fname=xl_rate_fname)

        #self.set_counts_by_protein()

        print('After all normalizations, annCounts.normalize_ann_counts():')
        print(self.raw_counts_df.iloc[:4, :4])
        print("And raw_counts_df.shape = ", self.raw_counts_df.shape)

    def get_total_reads(self, folder, load=False):

        if load and os.path.exists('./data/total_read_numbers.dill'):
            with open('./data/total_read_numbers.dill', 'rb') as f:
                self.total_counts = dill.load(f)
            return

        print("Calculating total read numbers from bed files...")
        self.total_counts = collections.defaultdict(int)
        for bedname in glob.glob(folder + '/*bed'):
            self.total_counts[bedname] = sum(1 for i in open(bedname, 'rb'))
            self.total_counts[os.path.basename(bedname).split('.bed')[0]] = self.total_counts[bedname]
        
        with open('./data/total_read_numbers.dill', 'wb') as f:
            dill.dump(self.total_counts, f)
        print("Finished calculating total read numbers.")

    def keep_only_these_proteins(self, protein_list):

        columns = [x for x in self.numeric_columns(self.raw_counts_df) \
            if (
                self.scheme.gene_from_fname(x) in protein_list
                ) or (
                self.scheme.gene_from_fname(x) in [
                    p.replace('-', ':') for p in protein_list]
                )]

        print('cols now ', columns)
        non_numeric = [x for x in self.raw_counts_df.columns if (
            self.raw_counts_df[x].dtype.kind not in 'bifc')]

        self.raw_counts_df = self.raw_counts_df[non_numeric+columns].copy()

        self.raw_counts = {}
        for row in self.raw_counts_df.to_dict('records'):
            self.raw_counts[row[self.indexing_col_name]] = row

    def bedfiles_for_a_protein(self, protein, df=None):
        if df is None:
            df = self.raw_counts_df
        return [x for x in self.numeric_columns(df) \
            if (self.scheme.gene_from_fname(x)==protein)]

    def normalize_counts_to_per_read(self, use_bed_dir=None, load_total_read_numbers=False):
        """Set self.counts_per_million_df from self.raw_counts_df and
        a folder of bed files (for total read count).
        """

        counts_cols = set([
            x for x in self.raw_counts_df.columns if self.scheme.gene_from_fname(x)])

        # Determine total read count for each replicate bed file.
        if use_bed_dir:
            self.get_total_reads(use_bed_dir, load=load_total_read_numbers)
        elif load_total_read_numbers:
            self.get_total_reads('', load=load_total_read_numbers)
        else:
            _ = self.raw_counts_df.loc[:, self.numeric_columns(self.raw_counts_df)].sum(axis=0).to_dict()
            self.total_counts = collections.defaultdict(int)
            self.total_counts.update(_)

        for bed_fname in [x for x in self.total_counts if x != 'gene_name']:
            if self.total_counts[bed_fname] == 0:
                print(f"No read counts in {bed_fname}. Blacklisting.")
                self.apply_blacklist(black_list=bed_fname)

        df = self.raw_counts_df.copy()

        for protein in self.proteins():

            for bed_fname in self.bedfiles_for_a_protein(protein):
                if self.total_counts[bed_fname] > 0:
                    df[bed_fname] = (1E6/self.total_counts[bed_fname]) * df[bed_fname]
                else:
                    df[bed_fname] = 0

        self.counts_per_million_df = df

        # Make a dict copy.
        self.counts_per_million_dict = {}
        for row in self.counts_per_million_df.to_dict('records'):
            self.counts_per_million_dict[row[self.indexing_col_name]] = row

    def normalize_counts_to_per_protein(
        self, xl_rate_fname='percentCrosslinked.xlsx'):
        """Sets counts_per_protein_df from counts_per_million_df and xl_rate_fname.
        """
        
        if not hasattr(self, 'counts_per_million_df'):
            self.normalize_counts_to_per_read()
        
        print("Multiplying by XL rates and outputing rates as PER 1E6 reads, PER 10000 proteins...")
        print("...(so if XL rate=0.01%, 3 reads per million stays 3. XL=1% -> 3 reads per million becomes 300.)")
        print('---.', self.proteins())
        df = self.counts_per_million_df.copy()

        for protein in self.proteins():
            for bed_fname in self.bedfiles_for_a_protein(protein):
                if self.total_counts[bed_fname] > 0:
                    xl_rate = self.scheme.percent_xl_of_protein_from_file(
                            protein, xl_rate_fname=xl_rate_fname)
                    print(f'Multiplying by {xl_rate} for {protein}.')
                    df[bed_fname] = 100. * float(xl_rate) * df[bed_fname].astype(float) 

        self.counts_per_protein_df = df

        # Make a dict copy.
        self.counts_per_protein_dict = {}
        for row in self.counts_per_protein_df.to_dict('records'):
            self.counts_per_protein_dict[row[self.indexing_col_name]] = row

        print("Finished multiplying by XL rates.")

    def set_counts_by_protein(self):
        """Make dicts ordered by protein.
        """
        
        print("Organizing protein replicates together...")
    
        self.raw_counts_by_protein = collections.defaultdict(dict)
        
        for gene_type, dict_of_beds in self.raw_counts.items():
            for bed_fname, value in dict_of_beds.items():
                
                protein = self.scheme.gene_from_fname(bed_fname)

                if protein:            
                    self.raw_counts_by_protein[gene_type].setdefault(protein, [])
                    self.raw_counts_by_protein[gene_type][protein].append(value)

        rows = []
        for rna, dict_by_protein in self.raw_counts_by_protein.items():
            rows.append({prot: np.mean(arr) for prot, arr in dict_by_protein.items()})
            rows[-1]['gene_name'] = rna

        self.raw_counts_by_protein_average_df = pandas.DataFrame(rows)

        if hasattr(self, 'counts_per_million_dict'):

            self.per_million_by_protein = collections.defaultdict(dict)

            for gene_type, dict_of_beds in self.counts_per_million_dict.items():
                for bed_fname, value in dict_of_beds.items():
                    
                    protein = self.scheme.gene_from_fname(bed_fname)
                    
                    if protein:
                        self.per_million_by_protein[gene_type].setdefault(protein, [])
                        self.per_million_by_protein[gene_type][protein].append(value)

            rows = []
            for rna, dict_by_protein in self.per_million_by_protein.items():
                rows.append({prot: np.mean(arr) for prot, arr in dict_by_protein.items()})
                rows[-1]['gene_name'] = rna
            self.per_million_by_protein_average_df = pandas.DataFrame(rows)

        if hasattr(self, 'counts_per_protein_dict'):
            
            self.per_protein_by_protein = collections.defaultdict(dict)

            for gene_type, dict_of_beds in self.counts_per_protein_dict.items():
                for bed_fname, value in dict_of_beds.items():

                    protein = self.scheme.gene_from_fname(bed_fname)

                    if protein:
                        self.per_protein_by_protein[gene_type].setdefault(protein, [])
                        self.per_protein_by_protein[gene_type][protein].append(value)                
                
            rows = []
            for rna, dict_by_protein in self.per_protein_by_protein.items():
                rows.append({prot: np.mean(arr) for prot, arr in dict_by_protein.items()})
                rows[-1]['gene_name'] = rna
            self.per_protein_by_protein_average_df = pandas.DataFrame(rows)

        print("Finished organizing protein replicates together.")