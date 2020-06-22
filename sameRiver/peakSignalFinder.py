import pandas, re, dill, HTSeq, collections, os, glob
import sameRiver
from sameRiver.biotypeAdder import biotypeAdder
from sameRiver.peakFinder import read_bedgraph, for_split_lines
from typing import List, Mapping, Union
import numpy as np

class heights():

    def __init__(self, peak_locations=Mapping[str, HTSeq.GenomicInterval], name=''):
        """Hold heights for various datasets at a given set of positions,
        the positions being given in a dict with names: {rna -> genomic location}.
        """

        # The protein that defines this peak set.
        self.name = name

        # peak_locations = {RNA name -> genomic interval}
        self.peak_locations = peak_locations

        # {protein -> {replicate -> height}}
        self.heights = {}

    def add_heights(self, ga: HTSeq.GenomicInterval, protein_name='', replicate_name=''):
        self.heights.setdefault(protein_name, {})
        self.heights[protein_name][replicate_name] = collections.defaultdict(float)

        for target, iv in self.peak_locations.items():
            if iv is None:
                self.heights[protein_name][replicate_name][target] = np.nan 
                continue
            for s, val in ga[iv].steps():
                self.heights[protein_name][replicate_name][target] += (s.end - s.start) * int(val)

    def average_bedgraphs(self):

        self.averaged = {}
        for protein_name in set(self.heights.keys()):
            df = pandas.DataFrame.from_dict(self.heights[protein_name])
            df.fillna(0., inplace=True)
            self.averaged[protein_name] = df.sum(axis=1).to_dict()

class peakSignalFinder():

    def __init__(
        self, scheme_filename: str,
        rnas_obj: sameRiver.set_of_named_mRNAs.set_of_named_mRNAs=None):
        self.scheme = sameRiver.scheme.scheme(scheme_filename)
        self.rnas_obj = rnas_obj
        # Dict {Protein => {RNA => position}}
        #self.peak_positions = locations
        # Dict {Protein => {RNA => height in raw reads}}
        #self.peak_heights = heights

    def gene_name_to_bedgraph_filenames(self, gene_name, bedgraph_filenames):
        return list(set([
            bedfname for bedfname in bedgraph_filenames \
            if (self.scheme.gene_from_fname(os.path.basename(bedfname)) == gene_name)
            ]))

    @staticmethod
    def str_to_iv(position):
        if type(position) != type(''):
            return   # Nan value.

        _ = re.split('[:/]', position)
        _[1] = int(_[1])
        bw = 50
        iv = HTSeq.GenomicInterval(_[0], max([_[1]-bw, 0]), _[1] + bw, _[2])
        return iv


    def find_signal_at_peaks(
        self, excel_of_peak_locations='data/peak_locations.xlsx',
        bedgraph_folder: str='', only_do: list=[],
        output_folder='./data/'):

        os.makedirs(os.path.dirname(output_folder), exist_ok=True)

        df = pandas.read_excel(excel_of_peak_locations, index_col=0)
        #df = df.loc[:, [(x<3) for x in range(len(df.columns))]]

        # Get HTSeq.GenomicIntervals for each peak location.
        set_of_all_positions = set()
        for prot in df.columns:
            set_of_all_positions |= set(df[prot].apply(self.str_to_iv).tolist()) - set([None])

        # Dict of {protein X -> {gene -> peak in gene for protein X}}
        all_positions = {}
        for prot in df.columns:
            all_positions[prot] = dict(zip(df.index, df[prot].apply(self.str_to_iv).tolist()))

        # If only loading signal/looking at peak locations for certain proteins.
        if only_do:
            all_positions = {k:v for k,v in all_positions.items() if k in only_do}

        # {proteinX -> heights object}, where the heights obj. holds locations for protein X's peaks.
        self._heights = {}

        # Organize the sets of peak locations by their defining protein.
        for protein, peak_locations in all_positions.items():
            self._heights[protein] = heights(peak_locations=peak_locations, name=protein)

        # For each protein we will load bedgraph data for:
        for n, (protein1, peak_locations1) in enumerate(all_positions.items(), start=1):
            print(f"Finding signal for {protein1} at all peak locations... {n}/{len(all_positions)}")

            # Load genomic bedgraph files for this protein.
            bedfnames1 = self.gene_name_to_bedgraph_filenames(
                protein1, glob.glob(bedgraph_folder + '/*_+.wig'))
            gas1 = {bedgraph_fname: read_bedgraph(bedgraph_fname) for bedgraph_fname in bedfnames1}

            if (not bedfnames1) or (len(bedfnames1) == 0):
                continue

            # For each set of peak locations:
            for protein2, peak_locations2 in all_positions.items():
                for bedfname1, ga1 in gas1.items():
                    # Match ann_counts name format (ensure this isn't a file path.)
                    col_name = bedfname1.split('/')[-1] 
                    col_name = re.sub('_[+-]\.wig', '', col_name)

                    self._heights[protein2].add_heights(
                        ga=ga1, protein_name=protein1, replicate_name=col_name)

        # For each protein we will output a file for:
        for protein, _heights in self._heights.items():

            # _heights.heights was set in this form:
            # heights[protein_name][replicate_name][target] = float

            _heights.average_bedgraphs()

            _log = collections.defaultdict(int)
            _log['Rows'] = len(df.index)
            #_log['RNAs'] = len(set(targets))

            output_df = pandas.DataFrame.from_dict(_heights.averaged)
            output_df['gene_name'] = output_df.index

            #if os.path.exists('./data/biotypeAdder.dill'):
            #    ba = biotypeAdder.load()
            #else:
            #    ba = biotypeAdder()
            #    ba.create()
            #    ba.save()

            #output_df = ba.add_biotypes_column_from_gene_name(output_df)

            out_fname = f'{output_folder}/signal_at_{protein}_peak_locations.xlsx'
            print(f"Writing to {out_fname}...")
            output_df.to_excel(out_fname)

            # End column block (one per protein).

    # End method function.