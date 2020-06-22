import os, sys, re, pandas, collections
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as scs
import scipy as sp

import sameRiver
import sameRiver.adapter_fluorescence_utils

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

def draw(fig, filename):
    fig.savefig(filename)
    plt.show()
    plt.clf()
    plt.close()

def hyperbolic(x, a, b):
    return (x*a)/(x+b)

def inverse_hyperbolic(x, a, b):
    return (-x*b)/(x-a)

def load_sheet(
    fname='',#'../Shift based quantification protocol and results.xlsx', 
    sheet_name='',#'180423_1'
    ):
    
    df = pandas.read_excel(fname, sheet_name=sheet_name)

    if ('Skip' in df.columns) and any([type(x) == type('') for x in df['Skip'].tolist()]):
        df = df[df['Skip']!='Yes'].copy()
        df = df[df['Skip']!='Skip'].copy()

    return df

def get_parameters(fname='180523.xlsx', sheet_name='180523_1'):

    df = load_sheet(
        fname=fname, sheet_name=sheet_name)
    
    #df, K = sameRiver.adapter_fluorescence_utils.get_Ka(df)
    df = df[df['fmols']<100]
    df = df[df['fmols']>3]

    _df = df.loc[[(x in ['αL3/L3', 'αL5/L5']) for x in df['Complex']], :]
    l5 = _df[_df['Object']=='αL5']
    l3 = _df[_df['Object']=='αL3']

    m5, b5 = slope_intercept(*fmols_sig(l5))
    l5['Est. fmols'] = [(x - b5)/m5 for x in l5['Signal']]
    m3, b3 = slope_intercept(*fmols_sig(l3))
    l3['Est. fmols'] = [(x - b3)/m3 for x in l3['Signal']]

    return [[m5, b5], [m3, b3]]

def get_parameters_and_draw_staple(fname='180523.xlsx', sheet_name='180523_1'):
    
    df = sameRiver.adapter_fluorescence_utils.load_sheet(fname=fname, sheet_name=sheet_name)
    
    plot_ip_control_signal_vs_fmols(df)
        
    staple = df[df['Complex']=='Staple']
    
    df = df[df['fmols']<300]
    
    l5 = df[df['Object']=='αL5']
    l3 = df[df['Object']=='αL3']

    # Linear fit.
    m5, b5 = slope_intercept(*fmols_sig(l5))
    l5['Est. fmols'] = [(x - b5)/m5 for x in l5['Signal']]
    m3, b3 = slope_intercept(*fmols_sig(l3))
    l3['Est. fmols'] = [(x - b3)/m3 for x in l3['Signal']]

    # Draw FU vs fmols for L5.
    fig = plt.figure()
    xmax = 250
    xarr = np.arange(start=0, stop=xmax)
    plt.plot(
        xarr, [b5 + m5*x for x in xarr], 'k-', linewidth=1)
    plt.plot(
        l5['fmols'], l5['Signal'], obj_color['L5'],
        markersize=14, marker='.', linewidth=0, alpha=0.5)

    #plt.plot(range(0, 250), [hyperbolic(x, *params5) for x in range(0, 250)], 'r:')
    plt.xlabel('Approximate fmols input')
    plt.ylabel('Fluorescence')
    plt.title('αL5 fluorescence vs input fmols')
    fig.set_figwidth(5)
    fig.set_figheight(5)
    draw(fig, './figs/staple_aL5_signal_vs_fmols.pdf')

    # Draw FU vs fmols for L3.
    fig = plt.figure()
    plt.plot(
        xarr, [b3 + m3*x for x in xarr], 'k-', linewidth=1)
    plt.plot(
        l3['fmols'], l3['Signal'], obj_color['L3'], 
        markersize=14, marker='.', linewidth=0, alpha=0.5)

    plt.xlabel('Approximate fmols input')
    plt.ylabel('Fluorescence')
    plt.title('αL3 fluorescence vs input fmols')
    fig.set_figwidth(5)
    fig.set_figheight(5)
    draw(fig, './figs/staple_aL3_signal_vs_fmols.pdf')

    # Draw linear fit vs fmols for L5/L3.
    
    fig = plt.figure()
    #plt.plot(range(0, xmax), range(0, xmax), 'k-', alpha=0.5, linewidth=1)
    plt.scatter(
        l5['fmols'], l5['Est. fmols'], c=obj_color['L5'],
        linewidth=0,
        s=100,
        #markersize=14, marker='.',
        alpha=0.3
        )
    plt.scatter(
        l3['fmols'], l3['Est. fmols'], c=obj_color['L3'],
        s=200, marker='X',
        linewidth=0,
        #markersize=14, marker='.',
        alpha=0.3
        )

    
    plt.xlim(0, xmax)
    plt.ylim(0, xmax)
    plt.xlabel('Approximate fmols input')
    plt.ylabel('Estimated fmols (linear fit)')
    fig.set_figwidth(5)
    fig.set_figheight(5)
    draw(fig, './figs/linear_estimated_fmols_vs_input_fmols.pdf')


def plot_ip_control_signal_vs_fmols(df):
    
    ip_control = df[df['Object']=='IP control'].copy()
    
    plt.clf()
    fig = plt.figure()
    sns.set_style('ticks')

    sns.lmplot(data=ip_control, y='Signal', x='fmols')
    
    plt.title("IP control signal vs fmols")
    plt.show()
    plt.clf()

def fmols_sig(_x):
    return _x['fmols'].tolist(), _x['Signal'].tolist()

def fmols_and_est(_x):
    return _x['fmols'].tolist(), _x['Est. fmols'].tolist()

def y_f(m, b, arr):
    return [val*m + b for val in arr]

def slope_intercept(x, y):
    a = scs.linregress(x, y)
    print("Linear fit: ", a)
    return a.slope, a.intercept

def est_fmols_linear_using_params(df, params, staples=None):

    [[m5, b5], [m3, b3]] = params

    l5 = df[df['Object']=='αL5']
    l3 = df[df['Object']=='αL3']

    if staples is not None:
        l5_pred_at_50 = b5 + 50 * m5
        l3_pred_at_50 = b3 + 50 * m3

        print("Predicted L3 signal at 50 fmols: {0}".format(l3_pred_at_50))

        l5_obs = l5[l5['fmols']==50]['Signal'].mean()
        l3_obs = l3[l3['fmols']==50]['Signal'].mean()
        scale_l5 = l5_pred_at_50/l5_obs
        scale_l3 = l3_pred_at_50/l3_obs

        print("-> actual L3 {} scale obs/pred {}.".format(l3_obs, scale_l3))

    else:
        scale_l5 = 1
        scale_l3 = 1

    def to_est(signal, object):
        if object=='αL5':
            return scale_l5 * (signal-b5)/m5
        if object=='αL3':
            return scale_l3 * (signal-b3)/m3
        return np.nan

    df['Est. fmols'] = [to_est(x, object) for x, object in zip(
        df['Signal'], df['Object'])]    

    return df

def est_fmols_linear(df, staples=None):

    _df = df.loc[[(x in ['αL3/L3', 'αL5/L5']) for x in df['Complex']], :].copy()

    l5 = _df[_df['Object']=='αL5'].copy()
    l3 = _df[_df['Object']=='αL3'].copy()

    l5 = l5[l5['fmols']>0].copy()
    l3 = l3[l3['fmols']>0].copy()

    # Linear fit.
    m5, b5 = slope_intercept(*fmols_sig(l5))
    m3, b3 = slope_intercept(*fmols_sig(l3))

    def to_est(signal, object):
        if object=='αL5':
            return (signal-b5)/m5
        if object=='αL3':
            return (signal-b3)/m3
        return np.nan

    df['Est. fmols'] = [to_est(x, object) for x, object in zip(
        df['Signal'], df['Object'])]

    return df

def est_fmols_hyperbolic(_x, params, channel=800, staples=None):
    
    scale = 1

    if staples is not None:
        pred_at_50 = hyperbolic(50, *params)
        print("Predicted signal at 50 fmols: {0}".format(pred_at_50))

        try:
            on_channel = staples[staples['Channel']==channel].copy()
            
            obs = on_channel[on_channel['fmols']==50]['Signal'].mean()
            dev = obs - pred_at_50
            scale = obs/pred_at_50
            print("-> actual {0} deviation {1} scale obs/pred {2}.".format(obs, dev, scale))
        except:
            print("No fmols at 50")
            scale = 1
        
    _x['Est. fmols'] = [
        1/scale * inverse_hyperbolic(signal, *params) for obj, signal in zip(
            _x['Object'], _x['Signal'])]

    return _x

