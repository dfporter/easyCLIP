
import os, sys
from pathlib import Path

#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

__all__ = [
    'RNA',
    'exp',
    'scheme_signal_RNAs',
    '__init__',
    'gapped_area',
    'set_of_named_mRNAs',
    '__pycache__',
    'gapped_area_mRNA',
    'signal_RNAs',
    'area',
    'makeBedDataFile',
    'signal_normalizer',
    'bedgraphs',
    'metaExp',
    'snoUtils',
    'biotypeUtils',
    'peakFinder',
    'stacked_bargraph',
    'refseqUtils',
    'statsForCounts',
    'collapse_duplicates',
    'rnaDataFileMaker',
    'countsO',
    'scheme',]

def rc(seq):
    _dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(_dict[x.upper()] for x in seq[::-1])


def add_paths(_paths):
    _paths.append('os.path.dirname(os.path.realpath(__file__))')
    return _paths
