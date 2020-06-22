import os, sys, re, pandas

class intronExonSplit():
    """Holds methods to split a pandas DataFrame by intron/exon.
    """

    @classmethod
    def split_exon_introns_mrna(cls, df):
        _df = cls.mRNA_only(df)
        (exons, introns) = cls.split_types(df)
        return exons, introns
    
    @staticmethod
    def mRNA_only(df):
        _df = df[df['Gene type']=='protein_coding'].copy()
        return _df
    
    @staticmethod
    def split_types(df, col=None):

        if col is None:
            names = df.index
        else:
            names = df[col].tolist()
        
        def later_part(_str):
            try:
                return _str.split('::')[1]
            except:
                return ''
            
        exons = [later_part(x) == 'exon' for x in names]
        introns = [later_part(x) == 'intron' for x in names]
        
        return df[exons].copy(), df[introns].copy()

    @classmethod
    def percent_introns(cls, counts, gene_list=None):
        #if col_list is None:
        #    col_list = df.columns
        df = counts.counts_df
        
        (exons, introns) = cls.split_exon_introns_mrna(df)
        
        exonic_back = total_reads_in_dataframe(exons[counts.cols_of_protein('No vector')])
        intronic_back = total_reads_in_dataframe(introns[counts.cols_of_protein('No vector')])        
        
        background_total = exonic_back + intronic_back
        
        results = []
        for gene in gene_list:
            exonic = total_reads_in_dataframe(exons[counts.cols_of_protein(gene)])
            intronic = total_reads_in_dataframe(introns[counts.cols_of_protein(gene)])
            total = exonic + intronic
            
            if (total < 1E3):
                continue
            norm_exon = exonic - exonic_back
            norm_intron = intronic - intronic_back
            norm_total = max([1, total - background_total])
            
            results.append({
                '% Exonic': 100 * exonic/total,
                '% Intronic': 100 *intronic/total,
                'Total': exonic + intronic,
                'Total (log10)': np.log10(exonic + intronic),
                'Intronic / total': (intronic/total) / np.log10(exonic + intronic),
                'Protein': gene,
                
                '% Exonic (- Back.)': 100 * norm_exon/total,
                '% Intronic (- Back.)': 100 * norm_intron/total,
                'Total (- Back.)': norm_total,
                'Total (log10) (- Back.)': np.log10(norm_total),               
                
            })
    
        return pandas.DataFrame(results)




def find_hits_columns(df, verbose=False):

    hits_columns = [col for col in df.columns if ((len(col.split('_'))>2) and (col != 'gene_name'))]
            
    if verbose:
        print("Found {0} HITS columns out of {1} total columns (had '_' pattern).".format(
            len(hits_columns), len(df.columns)))
        
    return hits_columns


def total_reads_in_dataframe(df):
    hits_columns = find_hits_columns(df, verbose=False)
    total = 0
    for col in hits_columns:
        total += df[col].sum()
    return total
    