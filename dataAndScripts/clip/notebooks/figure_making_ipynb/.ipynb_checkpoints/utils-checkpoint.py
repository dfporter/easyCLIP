import os, sys, re, pandas, collections
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as scs
import scipy as sp

obj_color = {
    'αL3': '#DF587A',
    'αL5': '#74C6A0',
    700: '#DF587A',
    800: '#74C6A0',
    'L5': '#DF587A',
    'L3': '#74C6A0',
    "5'": '#DF587A',
    "3'": '#74C6A0',
    'Hairpin': '#74C6A0',
    'Staple': 'k',
    '10 mM MgCl2': '#8DBFC6',
    '50 mM NaCl': '#5B8C88',
}

def load_sheet(
    fname='',#'../Shift based quantification protocol and results.xlsx', 
    sheet_name='',#'180423_1'
    ):
    
    df = pandas.read_excel(fname, sheet_name=sheet_name)

    if ('Skip' in df.columns) and any([type(x) == type('') for x in df['Skip'].tolist()]):
        df = df[df['Skip']!='Yes'].copy()
        df = df[df['Skip']!='Skip'].copy()
        
    df['FU'] = df['Signal']
    df['Signal/fmol'] = [x/y for x,y in zip(df.Signal, df.fmols)]
    df['FU/fmol'] = [x/y for x,y in zip(df.Signal, df.fmols)]
    df['log_Signal'] = [np.log2(np.max([x, 0.1])) for x in df.Signal]
    df['log_fmols'] = [np.log2(x) for x in df.fmols]
    return df

def get_Ka(df, v=False):

    staples = df[df['Complex']=='Staple']

    means = staples.groupby(by=['Object', 'fmols']).mean()

    if v:
        print(means)
        
    Ka5 = means.loc[('αL5', 50)]['Signal']/50

    # 4/3 correction for incomplete label.
    Ka3 = (4/3) * means.loc[('αL3', 50)]['Signal']/50

    K = {'αL3': Ka3, 'αL5': Ka5}
    
    def to_est_fmols(fu, obj):
        if obj in K:
            return fu/K[obj]
        else:
            return 0
    
    df['Est. fmols'] = [to_est_fmols(fu, obj) for fu,obj,fmols in zip(
        df['FU'], df['Object'], df['fmols'])]
    
    return df, K

def draw(fig, filename):
    fig.savefig(filename)
    plt.show()
    plt.clf()
    plt.close()
    