import pandas, re, random, collections, HTSeq, os, time, importlib

import numpy as np

import sklearn
from sklearn.neighbors.kde import KernelDensity
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
from numpy.random import choice


class peakSeqFinder():
    
    def __init__(self, fasta_dict=None):
        self.fasta_dict = fasta_dict
    
    def get_seqs(self, peak_locations_df):
            
        self.seqs = collections.defaultdict(dict)

        for protein in peak_locations_df.columns:

            self.seqs[protein] = {
                target:self.seq_around_point_of_highest_coverage(point) for target, point in zip(
                    peak_locations_df.index, peak_locations_df[protein])}

        self.nt_fractions()
        return self.seqs

    @staticmethod
    def _mk(name):
        if not os.path.exists(name):
            os.system('mkdir ' + name)
    
    @staticmethod
    def make_random_seqs(seqs, minimum=20000):
        
        outlines = ''
        alphabet = ['A', 'T', 'C', 'G']
        
        if len(seqs) > minimum:
            n_seqs = len(seqs)
        else:
            n_seqs = minimum

        if len(seqs) == 0:
            return ''
        
        seq_len = len(list(seqs.values())[0])

        for n in range(0, n_seqs):
            outlines += ">{}\n".format(n)
            outlines += "".join(random.choices(alphabet, k=seq_len)) + '\n'
            
        return outlines

    def make_scrambled_seqs(self, seqs, minimum=1000):

        outlines = ''
        if len(seqs) > minimum:
            _seqs = seqs.values()
        else:
            _seqs = np.random.choice(list(seqs.values()), size=minimum)

        for n, seq in enumerate(_seqs):
            seq = seq.upper()
            scrambled = ''.join(np.random.choice(list(seq), size=len(seq), replace=False))
            outlines += '>{}\n{}\n'.format(n, scrambled)

        return outlines

    def write_seqs_for_motif_searches(self, outfolder='tables/seqs/', write_controls=True):
        os.makedirs('tables/', exist_ok=True)
        os.makedirs(outfolder, exist_ok=True)
        
        for prot, target_seqs in self.seqs.items():
            print('Writing {}'.format(prot))
            target_seqs = dict([(k, v) for k, v in target_seqs.items() if len(v)])
            
            if len(target_seqs) == 0:
                print("No sequences for {}.".format(prot))
                continue
            
            if write_controls:
                randoms = self.make_random_seqs(target_seqs)
                with open(outfolder + '/randoms_for_{0}.fa'.format(prot), 'w') as f:
                    f.write(randoms)
                
            with open(outfolder + '/{0}.fa'.format(prot), 'w') as f:
                f.write('\n'.join(['>{0}\n'.format(k) + str(v) for k, v in target_seqs.items()]))
                
    def nt_fractions(self):
        for prot, target_seqs in self.seqs.items():

            all_seqs = ''.join(target_seqs.values())
            print(prot, self.frac_nt(all_seqs))
            
    def seq_around_point_of_highest_coverage(self, point):
        if point:

            try:
                point.strand
            except:
                
                try:
                    _str = point.split(':')
                    (pos, strand) = _str[-1].split('/')
                    point = HTSeq.GenomicPosition(_str[0], int(pos), strand)
                except:
                    return ''

            if point.strand == '+':
                iv_for_seq = HTSeq.GenomicInterval(
                    point.chrom, point.start - 15, point.end + 30, point.strand)
#                iv_for_seq.start -= 15  # 5
#                iv_for_seq.end += 30  # 20
            else:
                iv_for_seq = HTSeq.GenomicInterval(
                    point.chrom, point.start - 30, point.end + 15, point.strand)
#                iv_for_seq.start -= 30  # 20
#                iv_for_seq.end += 15  # 5

            return self.grab_sequence_from_iv_with_offset(iv_for_seq)
        else:
            return ''
        
    @staticmethod
    def rc(_seq):
        _dict = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        return ''.join([_dict[x.upper()] for x in _seq[::-1]])

    @staticmethod
    def frac_nt(_seq):
        import collections
        count = collections.defaultdict(int)
#        count = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
        
        for base in list(_seq):
            count[base] += 1
            
        total = sum(count.values())
        
        if total == 0:
            return count
            
        count = dict([(k, int(100 * v/total)) for k, v in count.items()])
        return count
        
    def grab_sequence_from_iv_with_offset(self, iv, offset=0):        
        
        seq = ''
        _chrom = iv.chrom
        if iv.chrom in self.fasta_dict:
            _chrom = iv.chrom
        elif 'chr' + str(iv.chrom) in self.fasta_dict:
            _chrom = 'chr' + str(iv.chrom)
        else:
            print("Chrom name {0} not found".format(iv.chrom))
            return ''

        seq = self.fasta_dict[_chrom][iv.start - 1 + offset:iv.end + offset].upper()

        if iv.strand == '-':
            seq = self.rc(seq)

        return seq