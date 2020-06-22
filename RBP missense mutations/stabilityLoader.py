import pandas, os, sys, re, collections, importlib, glob, itertools
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import numpy as np
from typing import List, Tuple, Union, Mapping
import nameUtils
importlib.reload(nameUtils)

def is_mut(_str: str) -> str:
    if len(_str.split(' ')) > 1:
        return 'Mut'
    return 'WT'

def load_data(input_filename: str, input_sheet: str) -> pandas.DataFrame:
    df = pandas.read_excel(input_filename, sheet_name=input_sheet)
    df['sig'] = df['Norm signal']

    df = df.loc[[bool(x=='Test') for x in df['Test or control']]]

    df = df.loc[:, [_ for _ in df.columns if _ not in ['Mut/WT', 'T Test P value', 'Averages', 'Notes']]]

    df['WT or mut'] = [is_mut(_) for _ in df.Protein]
    df = df.loc[[pandas.isna(x) for x in df.Discard], :]
    
    df['Protein group'] = [x.split(' ')[0] for x in df.Protein]

    wt_g = df.loc[[x=='WT' for x in df['WT or mut']], :]
    wt_g = wt_g.groupby(by=['Image', 'Protein group'])['sig'].mean()

    df['Abundance/WT normalized by image'] = [
        sig/wt_g.loc[im, prot] for sig, im, prot in zip(df.sig, df.Image, df['Protein group'])]
    
    return df

def cohend(x,y):
    return (np.mean(x) - np.mean(y)) / np.sqrt((np.std(x, ddof=1) ** 2 + np.std(y, ddof=1) ** 2) / 2.0)

def get_stats(
    df: pandas.DataFrame) -> Mapping[str, dict]:
    """For each set of proteins, get statistics for their comparison.
    Returns:
    {type of statistic: {mutant protein name: value}}
    """

    # Mapping[str, float]:
    signal_by_prot = df.groupby(by=['Protein'])['Abundance/WT normalized by image'].apply(list).to_dict()
    #print('>>>>', signal_by_prot)
    _stats = collections.defaultdict(dict)

    for wt in set(df['Protein group'].to_list()):

        wt_sig = signal_by_prot[wt]
        muts = set(df.loc[[bool(x==wt) for x in df.Reference], 'Protein'].to_list()) - set([wt])

        for mut in muts:

            mut_sig = signal_by_prot[mut]
            _stats['# stability replicates'][mut] = len(mut_sig)

            cohen = cohend(mut_sig, wt_sig)
            t, pval = sp.stats.ttest_ind(wt_sig, mut_sig)

            _stats['cohend'][mut] = cohen
            _stats['pval'][mut] = pval
            _stats['Protein group'][mut] = wt
            if pval < 0.01:
                _stats['ttest sig'][mut] = '<0.01'
            elif pval < 0.1:
                _stats['ttest sig'][mut] = '<0.1'
            elif pval > 0.1:
                _stats['ttest sig'][mut] = '>0.1'
                
    return _stats