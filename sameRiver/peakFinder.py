import pickle, pandas, re, random, collections, HTSeq, os, time, importlib, glob

import numpy as np
import sklearn
from sklearn.neighbors.kde import KernelDensity
import matplotlib.pyplot as plt
import seaborn as sns

import sameRiver
import sameRiver.artist
import sameRiver.scheme
import sameRiver.set_of_named_mRNAs
importlib.reload(sameRiver.artist)

def for_split_lines(fname, skip_first_line=False):
    with open(fname) as f:
        if skip_first_line:
            next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            if re.match('chr.*', s[0]):
                s[0] = s[0][3:]
            yield s


def read_bedgraph(fname):
    ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
    
    for s in for_split_lines(fname, skip_first_line=True):
        ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] += int(float(s[3]))

    minus_fname = fname.split('+.wig')[0] + '-.wig'
    
    for s in for_split_lines(minus_fname, skip_first_line=True):
        ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] += int(float(s[3]))

    return ga


class peakFinder():
    
    def __init__(
        self, scheme_filename: str,
        rnas_obj: sameRiver.set_of_named_mRNAs.set_of_named_mRNAs=None):

        self.scheme = sameRiver.scheme.scheme(scheme_filename)
        self.rnas_obj = rnas_obj
        self.rnas_obj.define_a_genomic_array_of_sets()
        self.priority = ['rRNA', 'snRNA', 'scaRNA', 'snoRNA', 'tRNA', 'protein_coding']

    def resolve_ambiguity(self, name, possible_features=[]):
        """Determine if reads in this region should be assigned to the gene "name" out of
        the overlapping genomic elements given in possible_features. {name}::exon or {name}::intron
        is expected to be in possible_features.
        name = e.g. "SNORD10", not "SNORD10::exon".
        possible_features = e.g. ['SNORD10::exon', 'RACK1::intron']
        """
        pri = 100

        _type ,= self.rnas_obj.mRNAs[name].transcript_biotypes  # Unpack a singleton set.
        (_type in self.priority) and (pri := self.priority.index(_type))
                    
        alt_pris = []

        # Does this region overlap with the exon of the GOI?
        cf_with_exon = bool(f"{name}::exon"  in possible_features)

        for alt_name in possible_features:

            (alt_gene, alt_ex_or_intron) = alt_name.split('::')

            if alt_gene == name:
                continue

            alt_pris.append(100)
            
            alt_type ,= self.rnas_obj.mRNAs[alt_gene].transcript_biotypes  # Unpack a singleton set.
            if (alt_type in self.priority):
                alt_pris[-1] = self.priority.index(alt_type)
        
            # If this is cf. with the GOI exon, de-priotitize introns.
            if cf_with_exon and (alt_ex_or_intron == 'intron'):
                alt_pris[-1] += 50
        
        if pri < min(alt_pris):
            return True

        return False        

    def genomic_iv_of_gene(self, name: str):
        
        if name not in self.rnas_obj.mRNAs:
            return None

        ivs = self.rnas_obj.mRNAs[name].elements['exon'][:]
        ivs = ivs + self.rnas_obj.mRNAs[name].elements['intron'][:]

        min_left = min([x.start for x in ivs])
        max_right = max([x.end for x in ivs])

        gene_iv = HTSeq.GenomicInterval(ivs[0].chrom, min_left, max_right, ivs[0].strand)

        list_of_ivs = []
        for iv, genes in self.rnas_obj.gaos[gene_iv].steps():

            if len(genes) == 1:
                list_of_ivs.append(iv)

            elif self.resolve_ambiguity(name, genes):
                list_of_ivs.append(iv)

        return list_of_ivs 
    
    @staticmethod
    def array_from_htseq(ga, iv):
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
    
    @classmethod
    def highest_point(cls, ga, iv, arr=None, gene_name='', **kwargs):

        if arr is None:
            arr = cls.array_from_htseq(ga, iv)  # Replace this with sensible code.
        
        if len(arr) > 0:
            _max = np.max(arr)
        else:
            _max = 0
            return False

        indexes_of_maxes = []
        for n, v in enumerate(arr):
            if v == _max:
                indexes_of_maxes.append(n)

        if len(indexes_of_maxes) > 1:
            consecutive = []
            for n, index in enumerate(indexes_of_maxes):
                if n == 0:
                    consecutive.append(index)
                    continue
                if index <= (consecutive[n-1]+3):  # Allow 2 nt gaps.
                    consecutive.append(index)
                else:  # Reject if there are nonadjacent maxima
                    return None
            # Take the peak as the middle point.
            max_index = int((consecutive[0] + consecutive[-1])/2)

#            max_index = random.choice(indexes_of_maxes)
        else:
            max_index = indexes_of_maxes[0]

        if type(iv) == type([]):
            genomic_point_of_highest_coverage = HTSeq.GenomicPosition(
                iv[0].chrom, iv[0].start + max_index, iv[0].strand)
        else:
            genomic_point_of_highest_coverage = HTSeq.GenomicPosition(
                iv.chrom, iv.start + max_index, iv.strand)

        return genomic_point_of_highest_coverage

    @classmethod
    def is_artifact(cls, arr: np.array, pos: int, window=50):
        _arr = arr[pos-window:pos+window+1]
        s = np.sum(_arr)
        fraction = [(_arr[n-1]+_arr[n])/s for n in range(1, len(_arr))]
        if any([f>0.8 for f in fraction]):  # Used 0.7 for the first set of motif identifications.
            return True
        return False

    @staticmethod
    def smooth(y, box_pts):
        #https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth   

    @classmethod
    def convolve(
        cls, ga, iv_list, make_fig=False, gene_name='',
        verbose=False):
            
        arr = cls.array_from_htseq(ga, iv_list)
        
        #print(f'pizza: arr={arr}')

        if len(arr) == 0:
            return None

        if len(arr) < 200:
            bw = 10
        elif len(arr) < 2000:
            bw = 20
        else:
            bw = 50

        y = cls.smooth(arr, bw)
 
        point = cls.highest_point(ga, iv_list, arr=y)  # Returns an HTSeq.GenomicPosition or None.

        if point is None:
            return None
        
        if cls.is_artifact(arr, int(point.pos-iv_list[0].start)):
            return None

        if make_fig:
            print(f'****\n{gene_name}')
            print('>point ', point)
            plt.plot(np.arange(len(arr)), y, 'k-')#, c='k-')
            # Red dot at the position of the called peak.
            plt.scatter(point.pos-iv_list[0].start, np.max(y), c='r')
            plt.show()
            plt.clf()
            plt.close()
        
        return point

    def gene_name_to_bedgraph_filenames(self, gene_name, bedgraph_filenames):
        return list(set([
            bedfname for bedfname in bedgraph_filenames \
            if (self.scheme.gene_from_fname(os.path.basename(bedfname)) == gene_name)
            ]))

    def find_peaks(
        self, excel_of_target_RNAs='tables/pvals.xlsx',
        bedgraph_folder: str='',
        datasets_included=None,
        method='highest_nt', verbose=False, only_do=None,
        top_n_targets=None, p_cutoff=0.001, outfile='./data/peak_locations.xlsx'):
        
        # Whether or not averaging is used:
        peak_method = {'highest_nt': self.highest_point,
                      'convolve': self.convolve}[method]

        # Read target excel.
        df = pandas.read_excel(excel_of_target_RNAs, index_col=0)

        # Pizza: Remobe non-snoRNA.
        #df = df.loc[['RACK1' in x for x in df.index], :]

        print("P values for determing targets are e.g.:")
        print(df.head(2))
        print(f"Shape of p value dataframe {df.shape}")

        overlaps_snoRNA = self.rnas_obj.RNAs_overlapping_snoRNA()
        print(f"{len(overlaps_snoRNA)} RNAs overlapping snoRNA (exons/introns counted separately).")
        #df = df.loc[[x.split('::')[0] not in overlaps_snoRNA for x in df.index], :]
        #print(f"Shape of p value dataframe after removing RNAs overlapping snoRNA {df.shape}")

        # Dict of {protein: {target_RNA: HTSeq.GenomicPosition, ...}, ...}
        self.peaks = collections.defaultdict(dict)
        #Dict of {protein: {target_RNA: float, ...}, ...}
        self.peak_heights = collections.defaultdict(dict)
        
        # For each list of RNA targets, corresponding to a CLIP'd protein:
        for protein_name in df.columns:

            # Get the bedgraph filenames corresponding to this protein.
            bedfnames = self.gene_name_to_bedgraph_filenames(
                protein_name, glob.glob(bedgraph_folder + '/*_+.wig'))

            # If only looking at a given subset of proteins...
            if (only_do is not None) and (protein_name not in only_do):
                continue
            
            print("=======" * 4, '\n', f"Finding peaks for {protein_name}...")

            # If only using a specified subset of bedgraph files...
            if datasets_included is not None:
                bedfnames = [x for x in bedfnames if \
                    os.path.basename(x).split('_+.wig')[0] in datasets_included]

            print(f"Bed files: {bedfnames}")

            # If no bedgraph files were found, move on.
            if (not bedfnames) or (len(bedfnames) == 0):
                continue

            gas = [read_bedgraph(bedgraph_fname) for bedgraph_fname in bedfnames]

            # Apply a p_cutoff if given, to determine target RNAs.
            if p_cutoff < 1:
                _df = df[df[protein_name]<=p_cutoff].copy()
            else:
                _df = df.copy()
                
            # Rank target RNAs.
            _df.sort_values(by=protein_name, ascending=True, inplace=True)

            # Reduce 'GAPDH::exon' to 'GAPDH' in targets.
            targets = set([x.split('::')[0] for x in list(_df.index)])

            print("# RNAs {0}".format(len(set(targets))))
            print('RNAs of starting df:', len(_df.index))

            _log = collections.defaultdict(int)

            _log['Rows below P value'] = len(_df.index)
            _log['RNAs'] = len(set(targets))

            rnas_processed = set()
                
            # For each RNA target of protein_name:
            for target in targets:

                # Intron and exon processing is identical, so don't do twice.  
                if target.split('::')[0] in rnas_processed:
                    continue
                
                if top_n_targets is not None:
                    if len(self.peaks[protein_name]) >= top_n_targets:
                        break
                    
                if verbose:
                    print('peakFinder.find_peaks(): {0}'.format(target))

                # iv actually a list of ivs - we subtract regions that overlap other genes.
                iv = self.genomic_iv_of_gene(target)
                #print(f'pizza: {target}: ivs = {iv}')

                if (iv is None) or len(iv) == 0:
                    continue
                
                if iv[-1].end - iv[0].start > 1E5:
                    _log['Gene too long'] += 1
                    continue
                
                # ga is the sum of the individual bedgraphs in the given interval.
                ga = HTSeq.GenomicArray('auto', stranded=True)

                # ga objects always have chr cut off.
                for _iv in iv:
                    _iv.chrom = re.sub('\Achr', '', _iv.chrom) 

                # Read the bedgraphs in the given interval.
                auc = 0
                for _ga in gas:
                    try:
                        for _iv in iv:
                            for step, val in _ga[_iv].steps():
                                ga[step] += int(val)
                                auc += int(val)
                    except:
                        _log['No reads found.'] += 1
                
                if auc < 10:
                    #print("pizza: had <10 reads.")
                    _log['Less than 10 reads'] += 1
                    continue

                #print("pizza: had >=10 reads.")

                # This is an HTSeq.GenomicPosition    
                self.peaks[protein_name][target] = peak_method(
                    ga, iv, make_fig=(random.randint(0,2000)==1), gene_name=target)
                
                if self.peaks[protein_name][target] is None:
                    _log['KDE method returned None'] += 1
                    del self.peaks[protein_name][target]
                else:
                    _log['Found peak'] += 1
                    pos = self.peaks[protein_name][target]
                    self.peak_heights[protein_name][target] = ga[pos]
                    
                if verbose:
                    print('peakFinder.find_peaks(): finished {0}'.format(target))
            
                rnas_processed.add(target.split('::')[0])

                if random.randint(0, 10000) == 1:
                    print('In progress, current log={}'.format(_log))
            print(_log)

            # End target block (one per target RNA).

        # End column block (one per protein).
        df = pandas.DataFrame.from_dict(self.peaks)
        df = df.astype(str)
        df.to_excel(outfile)

    # End method function.

    def save(self, fname='data/peaksFinder.dill'):
        if not os.path.exists('./data/'):
            os.system('mkdir ./data/')
        
        print("Saving to {}...".format(fname))
        with open(fname, 'wb') as f:
            dill.dump(self, f)
        print("...Saved.")

    @staticmethod
    def load(fname='data/peaksFinder.dill'):

        print("Loading {}...".format(fname))
        with open(fname, 'rb') as f:
            counts = dill.load(f)
        print("...Loaded.")
        return counts