import pandas, re, random, collections, HTSeq, os 
import operator, importlib, glob, functools, csv, random

import numpy as np
#import sklearn
#from sklearn.neighbors.kde import KernelDensity
#import matplotlib.pyplot as plt
#import seaborn as sns

import sameRiver
import sameRiver.artist
import sameRiver.scheme
import sameRiver.set_of_named_mRNAs
importlib.reload(sameRiver.artist)


def for_split_lines(fname, skip_first_line=False, only_chrom=None):

    with open(fname) as f:
        if skip_first_line:
            next(f)

        if only_chrom:
            for li in f:
                s = li.rstrip('\n').split('\t')
                if re.match('chr.*', s[0]):
                    s[0] = s[0][3:]
                if s[0] != only_chrom:
                    continue
                yield s                

        else:
            for li in f:
                s = li.rstrip('\n').split('\t')
                if re.match('chr.*', s[0]):
                    s[0] = s[0][3:]
                yield s


def read_bedgraph(fname, only_chrom=None):
    ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
    
    for s in for_split_lines(fname, skip_first_line=True, only_chrom=only_chrom):
        ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] += int(float(s[3]))

    minus_fname = fname.split('+.wig')[0] + '-.wig'
    
    for s in for_split_lines(minus_fname, skip_first_line=True, only_chrom=only_chrom):
        ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] += int(float(s[3]))

    return ga


class peakArtist():
    
    def __init__(
        self, scheme_filename: str, bedgraph_folder: str,
        rnas_obj: sameRiver.set_of_named_mRNAs.set_of_named_mRNAs=None,
        file_with_total_read_counts: str = '',
        ):

        self.scheme = sameRiver.scheme.scheme(scheme_filename)
        self.get_total_reads(file_with_total_read_counts)
        self.rnas_obj = rnas_obj
        self.rnas_obj.define_a_genomic_array_of_sets()
        self.priority = ['rRNA', 'snRNA', 'scaRNA', 'snoRNA', 'tRNA', 'protein_coding']
        self.bedgraph_folder = bedgraph_folder
        self.gas = {}

    def gene_name_to_bedgraph_filenames(self, gene_name, bedgraph_filenames):
        return list(set([
            bedfname for bedfname in bedgraph_filenames \
            if (self.scheme.gene_from_fname(os.path.basename(bedfname)) == gene_name)
            ]))

    def get_total_reads(self, fname=False):
        
        print(f"Loading total read counts from {fname}.")
        with open(fname) as fh:
            next(fh)  # Skip header.
            self.total_counts = dict((r[0],float(r[1])) for r in csv.reader(fh, delimiter='\t'))

    def signal_for_protein_in_region(
        self, protein_name, rna, iv=None):

        if iv is not None:  # Passed a genomic interval to plot.
            assert(type(iv) == type(HTSeq.GenomicInterval('a',1,20,'-')))
            _iv = iv
        # Passed a genomic interval as an RNA?
        elif type(rna) == type(HTSeq.GenomicInterval('a',1,20,'-')):
            _iv = rna
        else:  # Passed a gene name. Determine its genomic interval first.
            _iv = self.genomic_iv_of_gene(rna)

        bedfnames = self.gene_name_to_bedgraph_filenames(
                protein_name, glob.glob(self.bedgraph_folder + '/*_+.wig'))

        # Load any bedgraphs that haven't been loaded already.
        for bedgraph_fname in filter(lambda x: x not in self.gas, bedfnames):
            self.gas[bedgraph_fname] = read_bedgraph(bedgraph_fname)
        
        bn = lambda x: os.path.basename(x).rstrip('_+.wig')

        # Read the signal across bedgraphs.
        arrs = []
        for bedgraph_fname in bedfnames:
            print(f"Using {bedgraph_fname} as a {protein_name} dataset.")
            arr = self.array_from_htseq(self.gas[bedgraph_fname], _iv)
            print(f'total->{np.sum(arr)}')
            arr = 1E6 * arr/self.total_counts[bn(bedgraph_fname)]
            print(f'total->{np.sum(arr)}')
            arrs.append(arr)

        zer = np.zeros(len(arrs[0]))  # Array of zeroes the length of the genomic interval.
        # Add the signal across bedgraphs.
        _y = functools.reduce(operator.add, arrs, zer)
        # Make _y the averaged signal across bedgraphs.
        _y = _y/len(arrs)
        # Make an 'x' axis the length of the genomic interval.
        _x = np.linspace(0, len(zer), len(zer))

        return _x, _y

    def array_from_htseq(self, ga, iv):
        arr = []
        
        # If we are removing reads mapping to other genomic elements, iv is a list.
        # In between the intervals in iv are the other genomic elements, and we zero out those
        # areas.
        if type(iv) == type([]):
            for n, _iv in enumerate(iv):

                # We're assuming ivs are start inclusive and not end-inclusive: [start, end).
                if (n > 0) and (iv[n].start > iv[n-1].end):
                    arr += [0] * (iv[n].start - iv[n-1].end)

                for step_iv, value in ga[_iv].steps():
                    if (step_iv.end - step_iv.start) < 1E6:
                        arr += [value] * (step_iv.end - step_iv.start)

            return np.array(arr)
        
        # If passed a single HTSeq.GenomicInterval, not a list:
        for n, (step_iv, value) in enumerate(ga[iv].steps()):
            if (step_iv.end - step_iv.start) < 1E6:
                arr += [value] * (step_iv.end - step_iv.start)

        return np.array(arr)

    def genomic_iv_of_gene(self, name: str):
        
        if name not in self.rnas_obj.mRNAs:
            return None

        ivs = self.rnas_obj.mRNAs[name].elements['exon'][:]
        ivs = ivs + self.rnas_obj.mRNAs[name].elements['intron'][:]

        min_left = min([x.start for x in ivs])
        max_right = max([x.end for x in ivs])

        return HTSeq.GenomicInterval(re.sub('\Achr', '', ivs[0].chrom), min_left, max_right, ivs[0].strand)
