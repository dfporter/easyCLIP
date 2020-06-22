import sys, dill, pickle, pandas, re, HTSeq, os, importlib, operator, functools, random
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from matplotlib import patches

import sklearn
from sklearn.neighbors.kde import KernelDensity


def draw_ivs(ax, ivs, **kwargs):
    n = 0
    drawn = set()
    for (a, b) in ivs:
        n += 1
        if (a, b) in drawn:
            continue
        draw = True
        #for (_a, _b) in list(drawn):
        #    if (_a <= a <= _b) or (_a <= b <= _b):
        #        draw = False
        if draw:
            print(a, b)
            
            exon_y_len = 0.5
            center_y = 0.5*exon_y_len
            intron_y_len = 0.05
            intron = bool('intron' in kwargs)
            
            ax.add_patch(
                patches.Rectangle(
                (a, {True: center_y-intron_y_len/2, False:center_y-exon_y_len/2}[intron]),# Lower left corner x, y
                    b-a,  # width in x axis.
                    {True: intron_y_len, False:exon_y_len}[intron],  # length in y axis
                    #fill=intron,      # remove background
                    facecolor={True:'k', False:'#142C37'}[intron],
                 ))#np.linspace(a, b-1, b-a), [1] * (b-a), **kwargs)
            
            drawn.add((a,b))
    return ax
    
def draw_gene_body(rnas_object, name, intervals=None):
    #if name not in scheme_signal_rnas_obj.RNAs.mRNAs:
    #    return None

    exon_ivs = rnas_object.mRNAs[name].elements['exon']
    intron_ivs = rnas_object.mRNAs[name].elements['intron']
    all_ivs = exon_ivs + intron_ivs

    min_left = min([x.start for x in all_ivs])
    max_right = max([x.end for x in all_ivs])
    
    print(exon_ivs)
    print('{}:{}-{}'.format(exon_ivs[0].chrom, min_left, max_right))
    exon_ivs = [(_iv.start - min_left, _iv.end - min_left) for _iv in exon_ivs]
    intron_ivs = [(_iv.start - min_left, _iv.end - min_left) for _iv in intron_ivs]
    
    plt.clf()
    fig, ax = plt.subplots()
    ax = draw_ivs(ax, exon_ivs, linewidth=1)
    ax = draw_ivs(ax, intron_ivs, linewidth=1, c='k', intron=1)
    plt.ylim(-1, 4)
    plt.xlim(0-100, max_right-min_left+100)
    fig.savefig('/Users/dfporter/pma/dataAndScripts/clip/figs/gene_body_{}.svg'.format(name))
    plt.show()


def genomic_iv_of_gene(scheme_signal_rnas_obj, name):

    if name not in scheme_signal_rnas_obj.RNAs.mRNAs:
        return None

    ivs = scheme_signal_rnas_obj.RNAs.mRNAs[name].elements['exon'][:]
    ivs = ivs + scheme_signal_rnas_obj.RNAs.mRNAs[name].elements['intron']

    min_left = min([x.start for x in ivs])
    max_right = max([x.end for x in ivs])

    return HTSeq.GenomicInterval(ivs[0].chrom, min_left, max_right, ivs[0].strand)


def array_from_htseq(ga, iv):

    arr = []
    try:
        for n, (step_iv, value) in enumerate(ga[iv].steps()):
            if (step_iv.end - step_iv.start) < 1E6:
                arr += [value] * (step_iv.end - step_iv.start)
    except:
        print(iv.chrom)
        iv.chrom = re.sub('chr', '', iv.chrom)
        for n, (step_iv, value) in enumerate(ga[iv].steps()):
            if (step_iv.end - step_iv.start) < 1E6:
                arr += [value] * (step_iv.end - step_iv.start)
        
    return np.array(arr)


def signal_for_protein_in_region(
    protein, ssr, rna, bedfname_to_total_reads, iv=None):

    if iv is not None:
        assert(type(iv) == type(HTSeq.GenomicInterval('a',1,20,'-')))
        _iv = iv
    elif type(rna) == type(HTSeq.GenomicInterval('a',1,20,'-')):
        _iv = rna
    else:
        _iv = genomic_iv_of_gene(ssr, rna)

    _iv.chrom = _iv.chrom.replace('chr', '')
    arrs = []
    bed_fnames = ssr.gene_name_to_bed_filenames(protein)

    for bedfname in bed_fnames:
        bedO = ssr.beds.bedgraphs[bedfname]
        arr = array_from_htseq(bedO.ga, _iv)

        if _iv.chrom == 'repeats':

            sum_this_strand = np.sum(arr)
            alt_strand = dict([('+', '-'), ('-', '+')])[str(_iv.strand)]
            alt_iv = HTSeq.GenomicInterval(_iv.chrom, _iv.start, _iv.end, alt_strand)
            alt_arr = array_from_htseq(bedO.ga, alt_iv)
            if np.sum(alt_arr) > sum_this_strand:
                _iv = alt_iv
                arr = alt_arr
            print('{}: repeats chrom, using {} strand ({} cf {})'.format(
                protein, _iv.strand, sum_this_strand, np.sum(alt_arr)))

        arr = 1E6 * arr/bedfname_to_total_reads[bedfname]
        arrs.append(arr)

    zer = np.zeros(len(arrs[0]))
    _y = functools.reduce(operator.add, arrs, zer)
    print("{} Dividing by {}".format(protein, len(bed_fnames)))
    _y = _y/len(bed_fnames)
    _x = np.linspace(0, len(zer), len(zer))

    return _y, _x

    
def smooth(y, box_pts):
    #https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth   

def smooth_y_values(arr):
    #return np.linspace(0, len(arr)-1, len(arr)), smooth(arr, 30)
    pos_hits = []
    print('max value {}'.format(np.max(arr)))
    pos = [_ for _ in arr if _>0]
    min_pos = np.min(pos)
    
    if len(pos) == 0:
        print("No counts passed. Empty array given to smooth_y_values.")
        return 
    
    if min_pos >= 1:
        ratio_factor = 1
    else:
        ratio_factor = int(1/min_pos)
        
    print('ratio_factor = {}'.format(ratio_factor))
    
    #ratio_factor = np.min([ratio_factor, 1E])
    
    for pos, val in enumerate(arr):
        if val > 0:
            pos_hits.extend([pos] * int(ratio_factor*val))

    if len(pos_hits) > 10000:
        print('randoming: {}'.format(arr[:20]))
        pos_hits = random.sample(pos_hits, 10000)

    if len(arr) < 200:
        bw = 5
    elif len(arr) < 500:
        bw = 10
    else:
        bw = 100

    kde = KernelDensity(
        kernel='gaussian', bandwidth=bw).fit(np.array(pos_hits)[:, np.newaxis])

    x = np.linspace(0, len(arr)-1, len(arr))
    y = kde.score_samples(x[:, np.newaxis])
    y = [10**_ for _ in y]
    print("Sum value in arr: {}. Sum value in y: {}. Ratio: {}".format(
        np.sum(arr), np.sum(y), np.sum(arr)/np.sum(y)))
    ratio = np.sum(arr)/np.sum(y)
    
    y = [ratio * _ for _ in y]
    print("New sum of smoothed array: ", np.sum(y))
    return x, y

def draw3lines(_x, y_wt, y_100p, y_100q):
    linewidth=1#0.5
    a = 1
    plt.plot(_x, y_wt, label='WT', c='#142C37', linewidth=linewidth, alpha=a)
    plt.plot(_x, y_100p, label='L100P', c='#A5689A', linewidth=linewidth, alpha=a)
    plt.plot(_x, y_100q, label='L100Q', c='#899FC9', linewidth=linewidth, alpha=a)

def plot3(_x, y_wt, y_100p, y_100q, title='', rna='rna', smooth=False):
    fig = plt.figure()
    plt.clf()
    
    if smooth:
        (_x, sy_wt) = smooth_y_values(y_wt)
        (_x, sy_100p) = smooth_y_values(y_100p)
        (_x, sy_100q) = smooth_y_values(y_100q)
        draw3lines(_x, sy_wt, sy_100p, sy_100q)
    else:
        draw3lines(_x, y_wt, y_100p, y_100q)
    plt.legend()
    plt.title(title)
    fig.set_figheight(2)
    fig.set_figwidth(6)
    fig.savefig(
        '/Users/dfporter/pma/dataAndScripts/clip/figs/signal_at_{}_{}.pdf'.format(rna, title))
    plt.show()
    plt.clf()

#_x = [x for x in range(0, len(_y))]

