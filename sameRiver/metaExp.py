import collections, pandas, os, re, glob, sys, importlib, random, pickle
from pathlib import Path
from typing import List, Tuple, Union, Mapping

import sameRiver
import sameRiver.scheme
import sameRiver.exp

import sameRiver.makeBedDataFile
import sameRiver.statsForCounts

importlib.reload(sameRiver.exp)
importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.makeBedDataFile)
importlib.reload(sameRiver.statsForCounts)

class metaExp():
    """Combines sameRiver.exp objects, which each represent one sequencing run,
    into one meta experiment covering multiple sequencing runs.
    
    This is the master class for sequencing analysis after individual sequencing runs are
    split by barcode and mapped to the genome using STAR.
    
    The exp objects are made first using individually passed file paths for their
    respective sequencing data and scheme.xlsx info, and then added using metaExp.add_exp().
    metaExp is also passed folders paths for sequencing data and a scheme file path, like exp objects,
    but these are created by metaExp from the input exp objects.
    metaExp can then combine the exp data (metaExp.combine_exps()), outputing to the
    locations specified by file_paths.
    
    metaExp can combine exp files, using combine_exps, without any processing of the added
    exp objects, just concateninating information. In this case, it would be used:
    m = metaExp(file_paths={...})
    m.add_exp(...)
    m.add_exp(...)
    m.combine_exps()
    
    However, metaExp can also call the processing methods of the individual exp objects.
    In this case, it would be used:
    m = metaExp(file_paths={...})
    m.add_exp(...)
    m.add_exp(...)
    # Then any of the following:
    m.make_signal_data_file()
    # Or:
    maker = sameRiver.rnaDataFileMaker.rnaDataFileMaker()
    RNAs = maker.make_from_gtf_file(...)
    m.make_scheme_signal_RNA_data_files(rna_data_object=RNAs)
    # Or (this gets mapping/fastq numbers):
    m.determine_statistics()
    
    When combine_exps() is called, it looks for counts files to combine. The counts
    files can only be output after signal_RNAs data files have been made, as counts are
    determined using the signal and RNA location data held in signal_RNAs objects.
    As a result, combine_exps() assumes this processing has been done. So a complete
    analysis pipeline loads the individual exps into metaExp, calls make_scheme_signal_RNA_data_files()
    and then combine_exps().
    """
    
    def __init__(self, name='meta', file_paths: Mapping[str, str] = {}):
        self.exps = {}
        self.name = name
        self.file_paths = file_paths

        # Set default file paths.
        self.file_paths.setdefault('top_dir', os.path.abspath(os.path.dirname('.')))
        self.file_paths.setdefault('scheme', self.file_paths['top_dir'] + '/scheme.xlsx')
        self.file_paths.setdefault('counts', self.file_paths['top_dir'] + '/counts.txt')
        self.file_paths.setdefault('beds', self.file_paths['top_dir'] + '/beds/')
        self.file_paths.setdefault('wigs', self.file_paths['top_dir'] + '/beds/read_start/')

        # Make any needed directories.
        for _path in self.file_paths.values():
            if os.path.dirname(Path(_path)) != '':
                os.makedirs(os.path.dirname(Path(_path)), exist_ok=True)

        print('meta object created with: ', locals())

    def add_exp(self, exp_obj: sameRiver.exp.exp):
        
        if exp_obj.name in self.exps:
            print(f"Warning! metaExp already has an experiment named {exp_obj.name} - overwriting...")
            
        self.exps[exp_obj.name] = exp_obj
        
    def combine_exps(self, exp_list: Union[str, List[str]] = 'all'):
        
        if exp_list == 'all':
            exp_list = list(self.exps.keys())

        self.combine_scheme_files()
        self.combine_bed_files()
        self.combine_wig_files(exp_list)
        self.combine_counts_files(exp_list)
    
    def check_what_paths_exist(self):
        results = f"\n-**** {self.name} ***-\n"
        for name, exp in self.exps.items():
            results += exp.check_what_exists()

        with open('what_file_paths_exist.txt', 'w') as f:
            f.write(results)

    def make_bedgraphs_from_bed(
        self, use_first_n_lines=False, verbose=False, clobber=True, use='read start'):
        
        for name, exp in self.exps.items():
            exp.bed_to_wig(
                use_first_n_lines=use_first_n_lines, verbose=verbose, clobber=clobber, use=use)
            
    def determine_statistics(self, **kwargs):
        """This function modifies the exp.scheme dataframes with statistics
        on the fastq, bed, and bedgraph files. Pass a use_top_n keyword argument 
        to truncate processing of the fastq/bed/begraph files."""
        
        for name, exp in self.exps.items():
            if 'verbose' in kwargs and (kwargs['verbose']):
                print(f'Determining statistics for {name}')
            self.exps[name].determine_statistics(**kwargs)
            if 'verbose' in kwargs and (kwargs['verbose']):
                print(f'...Determined statistics for {name}')
                print('exp.scheme_df:', self.exps[name].scheme.scheme_df)
    
    def combine_scheme_files(self):
        schemes = {}
        scheme_list_of_rows = []
        
        for name, exp in self.exps.items():
            
            print(name, exp.file_paths['scheme'])
            
            # Load a scheme if needed.
            if not hasattr(exp, 'scheme'):
                schemes[name] = sameRiver.scheme.scheme(exp.file_paths['scheme'])
            else:
                schemes[name] = exp.scheme
            
            # Make sure P6_BC is a string.
            schemes[name].scheme_df['P6_BC'] = [
                str(x) for x in schemes[name].scheme_df['P6_BC'].tolist()]
            
            # Convert to dict to simplify combining them later.
            _d = schemes[name].scheme_df.to_dict('records')

            # Add the run code (e.g. 'rb' or 'm0')
            [x.update({'Seq run code': name}) for x in _d]

            scheme_list_of_rows.extend(_d)

        # Combine the list of dicts into a single dataframe.
        df = pandas.DataFrame(scheme_list_of_rows)

        # Write.
        os.makedirs(os.path.dirname(self.file_paths['scheme']), exist_ok=True)
        df.to_excel(self.file_paths['scheme'])
    
    def combine_scheme_stat_files(self):
        schemes = {}
        scheme_list_of_rows = []
        

        for name, exp in self.exps.items():

            scheme_with_stat = os.path.dirname(exp.file_paths['scheme']) + '/stats_in_scheme.xlsx'
            
            _d = sameRiver.scheme.scheme(scheme_with_stat).scheme_df.to_dict('records')

            [x.update({'Seq run code': name}) for x in _d]
            
            scheme_list_of_rows.extend(_d)

        df = pandas.DataFrame(scheme_list_of_rows)

        os.makedirs(os.path.dirname(self.file_paths['scheme']), exist_ok=True)
        df.to_excel(Path(os.path.dirname(self.file_paths['scheme']), 'stats_in_scheme.xlsx'))

    @staticmethod
    def proclaim(cmd: str, verbose=True):
        if verbose:
            print(cmd)
        os.system(cmd)

    @staticmethod
    def vprint(verbose, *args):
        if verbose:
            if len(args) == 1:
                print(args[0])
            else:
                print(args)
            
    def combine_bed_files(self, no_clobber=False, verbose=True):
        
        for name, exp in self.exps.items():
            
            self.vprint(verbose, name, exp.file_paths['beds'])
            os.makedirs(self.file_paths['beds'], exist_ok=True)

            # For each bed file in this exp's bed folder:
            for fname in glob.glob(exp.file_paths['beds'] + '/*.bed'):
                out_name = Path(self.file_paths['beds'], os.path.basename(fname))

                if not os.path.exists(self.file_paths['beds']):
                    os.system('mkdir ' + self.file_paths['beds'])
                
                if no_clobber and os.path.exists(out_name):
                    self.vprint(verbose, f"Already exists {out_name}, skipping...")
                    continue
                
                # Copy the bed file to the meta exp's bed folder.
                self.proclaim(f'rsync {fname} {out_name}', verbose=verbose)

    def combine_wig_files(self, no_clobber=False, verbose=False):
        
        for name, exp in self.exps.items():
            
            self.vprint(verbose, name, exp.file_paths['wigs'])
            os.makedirs(self.file_paths['wigs'], exist_ok=True)
            
            # For each bedgraph in this exp's begraph folder:
            for fname in glob.glob(exp.file_paths['wigs'] + '/*.wig'):
                out_name = Path(self.file_paths['wigs'], os.path.basename(fname))
                
                if no_clobber and os.path.exists(out_name):
                    self.vprint(verbose, f"Already exists: {out_name}; skipping...")
                    continue
                
                # Copy the wig file to the meta exp's bedgraph folder.
                self.proclaim(f'rsync {fname} {out_name}', verbose=verbose)         
    
    def combine_counts_files(self, verbose=False):

        # Create dict to hold the combined counts.
        all_counts = collections.defaultdict(dict)

        for name, exp in self.exps.items():

            self.vprint(verbose, name, exp.file_paths['counts'])

            # Load the scheme file.
            if not hasattr(exp, 'scheme'):
                exp.read_scheme()

            _counts = sameRiver.countsO.countsO(
                filename=exp.file_paths['counts'],
                scheme=exp.scheme,
                indexing_col_name='gene_name',
                index_col=0,
                exclude_unmapped=False,  # Removes _no_feature, _ambiguous if True.
                )
            
            _counts.simplify_column_names()

            # For each gene (GAPDH::exon, GAPDH::intron, ect.) and {bed_name -> read #}:
            for gene, dict_of_beds in _counts.raw_counts.items():
                
                if random.randint(0, 10000) == 1:
                    self.vprint(verbose, gene, dict_of_beds)

                # Add to the combined counts {gene -> {bed_name: read #}} dict.
                all_counts[gene].update(dict_of_beds)

        # Convert the {gene -> {bed_name: read #}} dict into a dataframe.
        df = pandas.DataFrame.from_dict(all_counts, orient='index')

        # Remove gene_name columns.
        df = df.loc[:, [x for x in df.columns if not (re.search('gene_name', x))]]
        
        # Write.
        df.to_csv(self.file_paths['counts'], sep='\t')

        return df

    def annotate_counts_file(self, annotate_individual_exps: bool = True,
        biotype_mapping_fname='enst_transcript_id_name_biotype_map.txt'):
        
        ann_counts_fname = os.path.dirname(self.file_paths['counts']) + \
            '/ann_' + os.path.basename(self.file_paths['counts'])
        self.file_paths['ann_counts'] = ann_counts_fname
        
        print("Annotating counts file and outputing to {}".format(ann_counts_fname))
                                                                                                
        self.counts = sameRiver.countsO.countsO(
            filename=self.file_paths['counts'],
            scheme=self.read_scheme(),
            index_col=0,
            indexing_col_name='gene_name')

        if 'xl_rate_fname' not in self.file_paths:
            print("No xl_rate_fname found in self.file_paths.")
            if os.path.exists('percentCrosslinked.xlsx'):
                print("Found ./percentCrosslinked.xlsx, using that.")
                self.file_paths['xl_rate_fname'] = 'percentCrosslinked.xlsx'
        
            else:
                print("Creating a dummy percentCrosslinked.xlsx file in the current directory.")
                rows = []
                for protein in list(self.scheme.proteins):
                    rows.append({'Exp': 'Exp', 'Protein': protein, 'Value': 100, 'Label': '% XL'})
                    rows.append({'Exp': 'Exp', 'Protein': protein, 'Value': 100, 'Label': '% XL (monomeric)'})
                pandas.DataFrame(rows).to_excel('percentCrosslinked.xlsx')
                self.file_paths['xl_rate_fname'] = 'percentCrosslinked.xlsx'

        self.counts.edit(
            verbose=True, xl_rate_fname=self.file_paths['xl_rate_fname'],
            mapping_fname=biotype_mapping_fname)

        self.counts.raw_counts_df.to_csv(ann_counts_fname, sep='\t')
#        pandas.DataFrame.from_dict(self.counts.raw_counts, orient='index').to_csv(
#            ann_counts_fname, sep='\t')
        
        # Annotate individual experiments:
        if annotate_individual_exps:
            for name, exp in self.exps.items():
                exp.annotate_counts_file()
        
    def make_signal_data_file(self, clobber=True, verbose=True):
        for exp_name in self.exps:
            self.vprint(verbose, "Making signal data file for {0}...".format(exp_name))
            self.exps[exp_name].make_signal_data_file(clobber=clobber)
            
    def make_scheme_signal_RNA_data_files(self, rna_data_object=None, no_clobber=False, verbose=True):
        
        if rna_data_object is None:
            self.vprint(verbose, "Loading ./data/rna.data...")
            rna_data_object = pickle.load(open('./data/rna.data', 'rb'))
            rna_data_object.find_introns()
        
        elif type(rna_data_object) == type(''):
            self.vprint(verbose, "Loading {}...".format(rna_data_object))
            rna_data_object = pickle.load(open(rna_data_object, 'rb'))
            rna_data_object.find_introns()

        for exp_name in self.exps:
            self.exps[exp_name].make_scheme_signal_RNA_data_files(
                rna_data_object=rna_data_object, no_clobber=no_clobber)
            
    def duplication_rate(self, outfname='tables/duplication_rates.txt'):
        
        reports = ""
        
        for exp_name in self.exps:
            reports += self.exps[exp_name].duplication_rate()
        
        if not os.path.exists(os.path.dirname(outfname)):
            os.system('mkdir {}'.format(os.path.dirname(outfname)))
        
        with open(outfname, 'w') as f:
            f.write(reports)

    def read_scheme(self, fname: Union[None, str] = None) -> sameRiver.scheme.scheme:
        
        if fname is not None:
            self.file_paths['scheme'] = fname
        self.scheme = sameRiver.scheme.scheme(self.file_paths['scheme'])
        
        return self.scheme
    
    def only_columns_with_known_proteins(self):
        
        print("Subsetting columns for annotated counts file based on whether they contain the CLIP'd protein.")
        
        if 'ann_counts' not in self.file_paths:
            ann_counts_fname = os.path.dirname(self.file_paths['counts']) + \
                '/ann_' + os.path.basename(self.file_paths['counts'])
            self.file_paths['ann_counts'] = ann_counts_fname
            
            print("The metaExp object had no known ann_counts file, so it was guessed to be {}.".format(ann_counts_fname))

            if not os.path.exists(self.file_paths['ann_counts']):
                print("...But this file did not exist.")
                return
            
        self.counts = sameRiver.countsO.countsO(
            filename=self.file_paths['ann_counts'],
            scheme=self.read_scheme(),
            index_col=0,
            indexing_col_name='gene_name')
        
        self.counts.only_columns_with_known_proteins(output_fname=self.file_paths['ann_counts'])
        
    def process(self, clobber=False, verbose=False):
        no_clobber = not(clobber)
        
        print(self.file_paths)
        self.combine_scheme_files()
        self.make_bedgraphs_from_bed(use_first_n_lines=False, verbose=verbose, clobber=clobber)
        self.make_signal_data_file(clobber=clobber, verbose=verbose)
        self.make_scheme_signal_RNA_data_files(rna_data_object=None, no_clobber=no_clobber, verbose=verbose)
        self.combine_bed_files(no_clobber=no_clobber, verbose=verbose)
        self.combine_wig_files(no_clobber=no_clobber, verbose=verbose)
        self.combine_counts_files(verbose=verbose)
        self.annotate_counts_file(annotate_individual_exps=False)
        #self.only_columns_with_known_proteins()
        self.save(clobber=clobber)
    
    def save(self, clobber=True):
        print("Saving to data/meta.data")
        if clobber or not(os.path.exists('data/meta.data')):
            with open('data/meta.data', 'wb') as f:
                pickle.dump(obj=self, file=f)

    def gather_repeats_files(self):

        for name, _exp in self.exps.items():
            name = name.replace('/', '-')

            # Inputs filenames.
            repEnrich_folder = Path(_exp.file_paths['sams'], 'RepEnrich_output/')
            family_fraction_fname = Path(repEnrich_folder, 'sample_name_prefix_family_fraction_counts.txt')
            fraction_fname = Path(repEnrich_folder, 'sample_name_prefix_fraction_counts.txt')
            
            # Output filenames.
            td = Path(self.file_paths['top_dir'], 'RepEnrich_output')
            os.makedirs(td, exist_ok=True)
            output_family_fname = Path(td, f'{name}_family_fraction_counts.txt')
            output_fname = Path(td, f'{name}_fraction_counts.txt')

            cmd = f"rsync {family_fraction_fname} {output_family_fname}"
            self.proclaim(cmd)
            cmd = f"rsync {fraction_fname} {output_fname}"
            self.proclaim(cmd)

    def split_repEnrich2_files_by_sample(self):
        """For each experiment, go through the lists of reads mapped to repeats and
        consolidate into a single file counting the number of reads in each sample
        mapped to each repeat.
        """

        counts_by_name = {}  # By sample name.
        counts_by_exp = {}  # To compare with the RepEnrich2 output.

        # Add the read to the correct sample based on the read name.
        def add_read(read, _gene):
            try:
                counts_by_name[_gene][read.split('rand=')[0]] += 1
            except:
                print(f"Failed in {_gene} on read {read} to parse out the read name.")

        # Initialize the dictionary counters with {repeat name: dict counter}.
        for name, _exp in self.exps.items():
            just_names = Path(_exp.file_paths['sams'], 'just_read_names')
            counts_by_name.update(
                {n: collections.defaultdict(int) for n in os.listdir(just_names)})
            counts_by_exp.update(
                {n: collections.defaultdict(int) for n in os.listdir(just_names)})

        # Add reads to their appropriate samples.
        for name, _exp in self.exps.items():

            print(f"\nReading {name}: ")
            dirname = Path(_exp.file_paths['sams'], 'just_read_names')

            for gene in os.listdir(dirname):
                print(f'{gene}... ', end='')
                fname = Path(dirname, gene)
                with open(fname) as f:
                    [add_read(read, gene) for read in f.readlines()]
                with open(fname) as f:
                    counts_by_exp[gene][_exp.name] = sum(1 for li in f)

        # Index is gene, column is sample.
        self.repeat_counts = pandas.DataFrame(counts_by_name).T
        # Index is gene, column is experiment name.
        self.repeat_counts_by_exp = pandas.DataFrame(counts_by_exp).T

        # Output to text files.
        outdir = Path(self.file_paths['top_dir'], 'data')
        os.makedirs(outdir, exist_ok=True)
        self.repeat_counts.to_csv(Path(outdir, 'repeat_counts.txt'), sep='\t')
        self.repeat_counts_by_exp.to_csv(Path(outdir, 'repeat_counts_by_exp.txt'), sep='\t')


"""
install.packages("readxl")
library("readxl")
library('edgeR')

# Get the repeats mapping data from file with columns as sample names,
# created with metaExp.split_repEnrich2_files_by_sample().
counts = read.csv(file='repeat_counts.txt', sep='\t', row.names='X')
counts[is.na(counts)] <- 0

# Look up library sizes:
s = read_excel('stats_in_scheme.xlsx')

# Get sample names.
n = unlist(s[,'long_fname_R1'])
n = unlist(strsplit(n, '.fastq'))

# Make a lookup structure.
x = structure(unlist(s[,'r1r2_split Reads']),names=n) #=unlist(s[,'long_fname_R1']))

# Put lookup structure into a function.
libsize = function(name) { return(x[name]) }

# Get a vector of library sizes in the order of colnames(counts).
n_reads = unlist(lapply(colnames(counts), libsize))

# Get rid of columns with NA read size.
counts = counts[!is.na(n_reads)]
n_reads = unlist(lapply(colnames(counts), libsize))

a = strsplit(colnames(counts), '_')  # colnames returns a vector.
get2 = function(x) { return(x[2]) }
proteins = lapply(a, get2)  # Returns a list.
proteins = unlist(proteins)  # Covert from list to vector.

meta = data.frame(row.names=colnames(counts), condition=proteins, libsize=n_reads)


# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)

# Build a DGE object for the GLM
y <- DGEList(counts=counts, lib.size=libsize)

"""

