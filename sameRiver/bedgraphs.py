import HTSeq, glob, re, os, random

"""
Example usage:
import importlib
import sameRiver.makeBedDataFile

bedMaker = sameRiver.makeBedDataFile.bedDataFileMaker()
beds = bedMaker.make_bed_data_file(
    bed_folder='/Users/dfporter/pma/miseq/Runs/170924_hiseq/beds/',
    output_data_filename='data/bed_170924.data',
    use='read start')

"""

class bedgraph():
    
    def __init__(self, fname='', from_list=False, **kwargs):
        
        if from_list:
            self.ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
            for t in from_list:
                self.ga[HTSeq.GenomicInterval(t[0], t[1], t[2], t[3])] += t[4]

        elif fname.split('.')[-1] == 'bed':
            self.ga = read_bed(fname=fname, **kwargs)
        else:
            self.ga = read_bedgraphs(fname=fname, **kwargs)
            
        self.fname = fname
        self.name = os.path.basename(fname).split('.wig')[0]
    
    def total_area_under_curve(self, iv='all'):
        total_auc = 0

        # If using all regions/chromosomes:
        if (type(iv) == type('')) and iv == 'all':
            for iv, value in self.ga.steps():
                total_auc += value * (iv.end - iv.start)

        # If using a specific region:
        else:
            for iv, value in self.ga[iv].steps():
                total_auc += value * (iv.end - iv.start)
                
        return total_auc            
        

class set_of_bedgraphs():
    
    def make_serializable(self):
        """Get rid of HTSeq objects, which have no __dict__ and can cause errors when pickling.
        Replace with a dict of {fname1: [[*iv, value], [*iv, value], ...], ..}
        """

        self.can_pickle = {}
        for fname in self.bedgraphs:
            self.can_pickle[fname] = []
            for iv, value in self.bedgraphs[fname].ga.steps():
                self.can_pickle[fname].append([iv.chrom, iv.start, iv.end, iv.strand, value])

        self.bedgraphs = self.can_pickle

    def load_serializable(self, _dict):
        for fname in _dict:
            self.bedgraphs[fname] = bedgraph(fname=fname, from_list=_dict[fname])

    def recover_bedgraph_objects(self):
        for fname in self.can_pickle:
            self.bedgraphs[fname] = bedgraph(fname=fname, from_list=self.can_pickle[fname])

    def __init__(self, bedgraphs_folder=None, bedgraphs_list=None, bed_folder=None,
                **kwargs):
        
        self.bedgraphs = {}
        
        print(f"Creating bedgrapsh.set_of_bedgraphs (**kwargs = {kwargs})...")
        if bedgraphs_list is not None:
            for fname in bedgraphs_list:
                if fname[-5:] == '-.wig':
                    continue
                self.bedgraphs[fname] = bedgraph(fname=fname, **kwargs)
        
        if bed_folder is not None:
            for fname in glob.glob(bed_folder + '/*.bed'):
                self.bedgraphs[fname] = bedgraph(fname, **kwargs)
                
            if len(glob.glob(bed_folder + '/*.bed')) < 1:
                print("WARNING: did not find any bed files in bed folder {}".format(bed_folder))
                
            self.bed_folder = bed_folder
            
        if bedgraphs_folder is not None:
            for fname in glob.glob(bedgraphs_folder + '/*_+.wig'):
                self.bedgraphs[fname] = bedgraph(fname=fname, **kwargs)
                
            if len(glob.glob(bedgraphs_folder + '/*_+.wig')) < 1:
                print("WARNING: did not find any wig/bedgraph files in bedgraphs folder {}".format(bedgraphs_folder))
                
    def total_area_under_curve(self, verbose=False):
        
        self.aucs = {}
        
        for fname, bed in self.bedgraphs.items():
            self.aucs[fname] = bed.total_area_under_curve()
            
            if verbose:
                print(f"{fname}: {self.aucs[fname]}")
            
    def collapse_duplicates_in_bed_dir(self, bed_folder=None):
        if bed_folder is None:
            bed_folder = self.bed_folder
        
        sameRiver.collapse_duplicates.collapse_bed_dir(bed_folder, 'temp_bed/')
        os.system('mv temp_bed/*bed ' + bed_folder + '/')
        #os.system('mv temp_bed/*/*bed ' + bed_folder + '/'
        
    def write_bedgraphs(self, top_dir='./data/wigs', verbose=False):
        
        if not os.path.exists(top_dir):
            os.system('mkdir ' + top_dir)
            
        for bedfname, bedgraph in self.bedgraphs.items():
            
            out_fname = top_dir + '/' + os.path.basename(bedfname).split('.')[0]
            out_fname = out_fname.replace('_+', '')
            out_fname = out_fname.replace('_-', '')
            bedgraph.ga.write_bedgraph_file(out_fname + '_+.wig', strand='+')
            bedgraph.ga.write_bedgraph_file(out_fname + '_-.wig', strand='-')
            

def for_split_lines(fname, skip_first_line=False):
    with open(fname) as f:
        if skip_first_line:
            next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            if re.match('chr.*', s[0]):
                s[0] = s[0][3:]
            yield s

def chromosomes(fname='/opt/genomes/gencode.v29/GRCh38.primary_assembly.genome.fa.gz.chrom_lengths'):
    """Not used."""
    chrom = {}
    with open(fname) as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            chrom[str(s[0])] = int(s[1])
            chrom[str(s[0]).replace('chr', '')] = int(s[1])
            chrom[str(s[0]).split('.')[0]] = int(s[1])
            chrom['repeats'] = 1E7  # Incorrect.
    return chrom

def read_bedgraphs(
        fname='/Users/dfporter/pma/miseq/Runs/170830/sams/consensus/GTCGTC_TCA_deletions_+.wig',
        **kwargs):
    
    print(f'read_bedgraphs({fname}), **kwargs={kwargs}')
    if 'verbose' in kwargs and kwargs['verbose']:
        print(f"Loading bedgraph file {fname}")
        
    if 'use_first_n_lines' in kwargs:
        use_first_n_lines = kwargs['use_first_n_lines']
        if use_first_n_lines is not None:
            print(f"bedgraphs.read_bedgraphs(): Using only the first {use_first_n_lines} lines of bedgraph file.")
    else:
        use_first_n_lines = False
        
    if use_first_n_lines is None:
        use_first_n_lines = False
        
    ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
    
    if (use_first_n_lines is False) or (use_first_n_lines is None):
        for s in for_split_lines(fname, skip_first_line=True):
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] += int(float(s[3]))
    else:
        for line_n, s in enumerate(for_split_lines(fname, skip_first_line=True)):
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] += int(float(s[3]))
            
            if line_n >= use_first_n_lines:
                break  # If you found this comment, email me and I will mail you $20.
                
    minus_fname = fname.split('+.wig')[0] + '-.wig'
    
    if use_first_n_lines is False  or (use_first_n_lines is None):
        for s in for_split_lines(minus_fname, skip_first_line=True):
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] += int(float(s[3]))

    else:
        for line_n, s in enumerate(for_split_lines(minus_fname, skip_first_line=True)):
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] += int(float(s[3]))
            
            if line_n >= use_first_n_lines:
                break
                
    return ga


def for_split_bed_lines(fname):
    with open(fname) as f:
        for li in f:
            sp = li.rstrip('\n').split('\t')
            yield (sp[0], int(sp[1]) , int(sp[2]), sp[5])
            
            
def read_bed(fname='', use_first_n_lines=False, **kwargs):
    """Bed format is 0-based and [a,b).
    So a read on the + strand has it's 5' end at [a], and a
    read on the - strand has its read at [b-1].
    """
    ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
    
    if ('verbose' in kwargs) and (kwargs['verbose']):
        print('Loading bed file {0}...'.format(fname))
    
    if 'use' not in kwargs:
        use = 'read start'
    else:
        use = kwargs['use']
    
    # If using the first n lines, read the bed:
    if (use_first_n_lines) and (use_first_n_lines>=1):
        if use == 'bed coord':
            for n, s in enumerate(for_split_bed_lines(fname)):
                ga[HTSeq.GenomicInterval(s[0], s[1], s[2], s[3])] += 1
                
                if n >= use_first_n_lines:
                    break

        if use == 'read start':
            for n, s in enumerate(for_split_bed_lines(fname)):
                if s[3] == '+':
                    ga[HTSeq.GenomicPosition(s[0], s[1], s[3])] += 1
                else:
                    ga[HTSeq.GenomicPosition(s[0], s[2]-1, s[3])] += 1
                    
                if n >= use_first_n_lines:
                    break
                    
    # If using the entire bed file, read the bed:
    else:
        if use == 'bed coord':
            for s in for_split_bed_lines(fname):
                ga[HTSeq.GenomicInterval(s[0], s[1], s[2], s[3])] += 1

        if use == 'read start':
            for s in for_split_bed_lines(fname):
                if s[3] == '+':
                    ga[HTSeq.GenomicPosition(s[0], s[1], s[3])] += 1
                else:
                    ga[HTSeq.GenomicPosition(s[0], s[2]-1, s[3])] += 1
                        
    return ga
