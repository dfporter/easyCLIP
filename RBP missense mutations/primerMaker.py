import pandas, os, sys, re, collections
import Bio
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
        
class primerMaker():
    """
    Generates InFusion stitching primers for a set of inserts.
    
    Arguments:
    1. Left and right overlap regions with vector
    2. A list of sequences to be stitched together.
    
    Outputs:
    1. A set of primers
    """
    
    def __init__(self, left='', right='', blocks=[], fname=None, overlap_len=9):
        self.left, self.right = (left.upper(), right.upper())
        self.blocks = [x.upper() for x in blocks]
        self.overlap_len = overlap_len
        self.fname = fname
        self.seq = ''.join(blocks)
        #self.pick_primers()

        
    def pick_primers(self):
        self.primers = {'No overlaps': {}, 'stitch': {}}
        self.p = self.primers
        for n, seq in enumerate(self.blocks):
            
            is_first, is_last = (0, 0)
            if n == 0:
                is_first = 1
            if n + 1 == len(self.blocks):
                is_last = 1
                
            left_name = str(n) + '_left'
            right_name = str(n) + '_right'
            
            self.p['No overlaps'][left_name], self.p['No overlaps'][right_name] = (
                      self.find_left_primer(seq), self.find_right_primer(seq))
            
            if is_first and is_last:
                self.p['stitch'][left_name], self.p['stitch'][right_name] = (
                      self.left + self.find_left_primer(seq), self._rc(self.right) + self.find_right_primer(seq))          
                
            if is_first and not is_last:
                self.p['stitch'][left_name], self.p['stitch'][right_name] = (
                    self.left + self.find_left_primer(seq), 
                    self.reverse_left_overlap(self.blocks[n+1]) + self.find_right_primer(seq))

            if (not is_first) and (not is_last):
                self.p['stitch'][left_name], self.p['stitch'][right_name] = (
                    self.for_right_overlap(self.blocks[n-1]) + self.find_left_primer(seq), 
                    self.reverse_left_overlap(self.blocks[n+1]) + self.find_right_primer(seq))
            
            if (not is_first) and is_last:
                self.p['stitch'][left_name], self.p['stitch'][right_name] = (
                    self.for_right_overlap(self.blocks[n-1]) + self.find_left_primer(seq), 
                    self._rc(self.right) + self.find_right_primer(seq))
        #pp(self.p)
    
    def _rc(self, seq):
        rc = Seq(seq).reverse_complement()
        return rc._data
    
    def for_right_overlap(self, seq):
        return seq[-self.overlap_len:]
    
    def reverse_left_overlap(self, seq):
        rc = Seq(seq).reverse_complement()
        return rc._data[-self.overlap_len:]
    
    def find_left_primer(self, seq, optimal_tm=54):
        seqO = Seq(seq)
        seqs = []
        for _len in range(10, 60):
            if _len <= len(seq):
                seqO = Seq(seq[:_len])
                seqs.append([seq[:_len], abs(optimal_tm - mt.Tm_NN(seqO))])
            else:
                seqO = Seq(seq)
                seqs.append([seq, abs(optimal_tm - mt.Tm_NN(seqO))])
                break
        seqs = sorted(seqs, key=lambda x: x[1])
        best = seqs[0]
        return best[0]
    
    def find_right_primer(self, seq):
        rc = Seq(seq).reverse_complement()
        return self.find_left_primer(rc._data)
    
    def primers_fasta(self, prefix=''):
        out = ''
        name = ''
        for top_level_name, vals in self.p.items():
            for second_level_name, seq in vals.items():
                name = prefix + '_' + second_level_name + top_level_name
                out += '>{0}\n{1}\n'.format(name, seq)
        return out
    
    def fasta(self, prefix=''):
        out = self.primers_fasta(prefix=prefix)
        out += '>' + prefix + '_stitched_seq\n'
        out += ''.join([self.left, ''.join(self.blocks), self.right])
        #print(out)
        return(out)

    def get_seq_as_fasta(self, prefix=''):
        return '>' + prefix + '_stitched_seq\n' + self.seq
    
    def one_per_line(self, prefix=''):
        out = ''
        name = ''
        for top_level_name, vals in self.p.items():
            for second_level_name, seq in vals.items():
                name = prefix + '_' + second_level_name + top_level_name
                out += '{0}\t{1}\n'.format(name, seq)
        #print(out)
        return out        
            