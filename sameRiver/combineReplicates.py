import collections, pandas
import numpy as np

def combine_protein_replicates(df, nb, scheme):
    """ Averages columns of the same protein.
    df: pandas.DataFrame.
    nb: statsForCountsNB.
    scheme: schemeO.
    """
    # Combine protein replicates:
    by_prot = collections.defaultdict(list)
    for col in nb.positives.numeric_columns(df):
        prot = scheme.gene_from_fname(col)
        by_prot[prot].append(df[col].values)

    for prot in by_prot:
        arrs = by_prot[prot]
        n = len(by_prot[prot])
        by_prot[prot] = np.sum(arrs, axis=0)/n

    if 'gene_name' in df.columns:
        by_prot['gene_name'] = df['gene_name']
    if 'Gene type' in df.columns:
        by_prot['Gene type'] = df['Gene type']
    
    df = pandas.DataFrame(by_prot)

    return df