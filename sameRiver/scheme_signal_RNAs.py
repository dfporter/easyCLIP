import importlib

import sameRiver
import sameRiver.signal_RNAs
import sameRiver.scheme

importlib.reload(sameRiver.signal_RNAs)

class scheme_signal_RNAs(sameRiver.signal_RNAs.signal_RNAs):
    """Useful to output counts.txt files from RNA and bed data."""
    
    def __init__(self, set_of_bedgraphs=None, set_of_RNAs=None):
        if (set_of_bedgraphs is not None) and (set_of_RNAs is not None):
            super().__init__(set_of_bedgraphs=set_of_bedgraphs, set_of_RNAs=set_of_RNAs)
    
    def add_signal_RNAs(self, signal_RNAs_object=None, set_of_bedgraphs=None, set_of_RNAs=None):
        
        if signal_RNAs_object is not None:
            self.__dict__.update(signal_RNAs_object.__dict__)
        else:
            super().__init__(set_of_bedgraphs=set_of_bedgraphs, set_of_RNAs=set_of_RNAs)
    
    def add_scheme(self, scheme_fname=''):
        self.scheme = sameRiver.scheme.scheme(scheme_fname=scheme_fname)
    
    def gene_name_to_bed_filenames(self, gene_name):
        return list(set([
            bedfname for bedfname in self.beds.bedgraphs if (self.scheme.gene_from_fname(bedfname) == gene_name)
                    ]))
        
