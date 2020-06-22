import pandas, re, collections, os
import sameRiver
import sameRiver.biotypeUtils

class biotyper():
    """Biotyping methods for sameRiver.countsO objects.
    Initiate with a countsO object.s
    """

    def __init__(self, countsO=None):
#        if countsO is not None:
#            assert(type(countsO) == type(sameRiver.countsO.countsO()))
        self.countsO = countsO
        self.counts_df = countsO.counts_df

    def cols_of_protein(self, protein):
        """Function identical in countsO."""
        return [col for col in self.counts_df.columns if re.search(protein, col)]

    def total_reads_per_protein(self, df):

        totals = collections.defaultdict(int)
        for protein in self.countsO.proteins():
            for col in self.cols_of_protein(protein):
                if col not in df:
                    totals[protein] += 0
                else:
                    totals[protein] += df[col].sum()
                
        return totals

    @staticmethod
    def total_reads_in_dataframe(df):
        hits_columns = find_hits_columns(df, verbose=False)
        total = 0
        for col in hits_columns:
            total += df[col].sum()
        return total

    @staticmethod
    def _genes(df):
        hits_columns = find_hits_columns(df)
        return set([x.split('_')[0] for x in hits_columns])

    @staticmethod
    def _cols_of_gene(df, gene):
        cols_to_keep = set()
        hits_columns = find_hits_columns(df, verbose=False)
        for col in hits_columns:
            if col.split('_')[0] == gene:
                cols_to_keep.add(col)
        return list(cols_to_keep)

    @staticmethod
    def biotypes_with_enough_reads(df, cutoff=50000):
        biotypes_to_keep = set()
        hits_columns = find_hits_columns(df)
        for biotype in set(df['Gene type'].tolist()):
            
            if biotype == 'processed_transcript':
                continue
            if biotype == 'retained_intron':
                continue
            
            all_for_type = 0
            sub = df[df['Gene type']==biotype]
            
            for col in hits_columns:
                all_for_type += sub[col].sum()
            if all_for_type >= cutoff:
                biotypes_to_keep.add(biotype)
        
        biotypes_to_keep -= set([
            'antisense', 'nonsense_mediated_decay', 
        'transcribed_processed_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene', 'sense_overlapping',
            'transcribed_unprocessed_pseudogene', 'TEC'
        ])
        print("Biotypes with at least {} reads across all datasets (Ignoring 'Unknown'): {}.".format(
            cutoff, len(biotypes_to_keep - set(['Unknown']))))
        
        return biotypes_to_keep - set(['Unknown'])

    def fraction_by_biotype(self, df=None):

        if df is None:
            df = self.counts_df

        for_df = []
        usable_types = self.biotypes_with_enough_reads(df)
        
        sub = df[[(x in usable_types) for x in df['Gene type']]].copy()

        totals = self.total_reads_per_protein(sub)
        
        for biotype in self.biotypes_with_enough_reads(df):
            
            sub = df[df['Gene type']==biotype]

            for gene in self.countsO.proteins():

                if totals[gene] == 0:
                    continue

                try:
                    reads_in_biotype = self.total_reads_in_dataframe(
                        sub[self.cols_of_protein(gene)])
                    for_df.append(
                        {'Gene type': biotype, 'Protein': gene,
                         'Reads': reads_in_biotype,
                         'Reads (% total)': 100 * reads_in_biotype/totals[gene]})
                except:
                    pass
                
        return pandas.DataFrame(for_df)
