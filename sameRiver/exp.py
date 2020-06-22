import collections, pandas, os, re, glob, sys, importlib, pickle, subprocess, time, dill
from typing import List, Tuple, Union
from pprint import pprint
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import sameRiver
import sameRiver.scheme
import sameRiver.makeBedDataFile
import sameRiver.scheme_signal_RNAs
import sameRiver.bedAndFastaStats
import sameRiver.bedgraphs
import sameRiver.countsO
import sameRiver.clip_adapters
import sameRiver.combine_fastqs_for_mapping
import sameRiver.sam_to_bed_and_wig
import sameRiver.remove_short_reads
import sameRiver.collapse_duplicates
import sameRiver.statsForCounts
import sameRiver.split_r1
import sameRiver.split_r2
import sameRiver.mapping

importlib.reload(sameRiver.scheme_signal_RNAs)
importlib.reload(sameRiver.makeBedDataFile)
importlib.reload(sameRiver.bedAndFastaStats)
importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.bedgraphs)
importlib.reload(sameRiver.countsO)
importlib.reload(sameRiver.clip_adapters)
importlib.reload(sameRiver.combine_fastqs_for_mapping)
importlib.reload(sameRiver.sam_to_bed_and_wig)
importlib.reload(sameRiver.remove_short_reads)
importlib.reload(sameRiver.collapse_duplicates)
importlib.reload(sameRiver.statsForCounts)
importlib.reload(sameRiver.split_r1)
importlib.reload(sameRiver.split_r2)
importlib.reload(sameRiver.mapping)

class exp(sameRiver.mapping.mappingMethods):
    """
    This is a holder of scheme, bed, count, and rna data locations.
    It differs from scheme_signal_RNAs objects, and the like, in that it
    only points to signal data locations, and never contains signal
    data, except during the creation of bed and scheme_signal_RNAs .data files, after
    which the large variables should go out of scope and be cleared from memory.
    
    This is intended to be used for dataset combinations using metaExp data objects.
    """
    
    def __init__(self, name: str = '0', file_paths: dict = {}):
        self.name = name

        assert(type(file_paths) == type({}))

        self.log = {'__init__ arguments': locals()}

        self.file_paths = file_paths

        # If 'scheme' and 'fastq' are given, guesses anything not already set.
        self.set_default_paths_if_needed()

        if not os.path.exists(self.file_paths['R1_fastq']):
            raise IOError(f"Need input fastq file (R1_fastq) in file paths. Guessed {self.file_paths['R1_fastq']} but no such file exists.")
        if not os.path.exists(self.file_paths['R2_fastq']):
            raise IOError(f"Need input fastq file (R1_fastq) in file paths. Guessed {self.file_paths['R2_fastq']} but no such file exists.")

        with open(self.file_paths['log'], 'w') as f:
            pprint(self.log, f)

    def write_to_log(self):
        with open(self.file_paths['log'], 'a') as f:
            f.write('\n--***--\n')
            pprint(self.log, f)

    def write_read_len_stats(self):
        stats_fname = os.path.dirname(file_paths['log']) + '/read_lens.txt'

    def clean_intermediates(self):
        for fname in glob.glob(self.file_paths['r1_split'] + '/*.fastq'):
            subprocess.check_output(f'gzip {fname}'.split(' '))

        for fname in glob.glob(self.file_paths['r1_clipped'] + '/*.fastq'):
            subprocess.check_output(f'gzip {fname}'.split(' '))

        # Should clean more...

    def check_what_exists(self) -> str:
        to_finder = {
            'Fastq': lambda x: all([os.path.exists(x.file_paths['R1_fastq']),
                os.path.exists(x.file_paths['R2_fastq'])]),
            'Scheme': lambda x: os.path.exists(x.file_paths['scheme']),
            'r1r2_split': lambda x: os.path.exists(x.file_paths['r1r2_split']),
            'r1r2_clipped': lambda x: os.path.exists(x.file_paths['r1r2_clipped']),
            'cutadapt split': lambda x: os.path.exists(
                x.file_paths['fastq'] + '/ready_to_map/cutadapt/split/'),
            'beds': lambda x: os.path.exists(x.file_paths['beds']),
            'bedgraphs': lambda x: os.path.exists(x.file_paths['wigs']),
            'scheme_with_stat': lambda x: os.path.exists(
                os.path.dirname(x.file_paths['scheme']) + '/stats_in_scheme.xlsx'),
        }

        results = f"\n-*-{self.name}\t{os.path.dirname(self.file_paths['scheme'])}-*-\n"
        for name, exists in to_finder.items():
            result = {True: '[x]', False: '[ ]'}[exists(self)]
            results += f'\t{result}\t{name} exists?\n'
        results += '-===-\n'

        print(results)

        return results

    @staticmethod
    def autogenerate_paths(top_dir: str) -> dict:
        file_paths = {
            'scheme': top_dir + '/scheme.xlsx',
            'beds': top_dir + '/beds/',
            'wigs': top_dir + '/beds/read_start/',
            'fastq': top_dir + '/fastq/',
            'R1_fastq': top_dir + '/fastq/raw/R1.fastq',
            'R2_fastq': top_dir + '/fastq/raw/R2.fastq',
            'sams': top_dir + '/sams/',
            'counts': top_dir + '/counts.txt',
            'log': top_dir + '/log.txt',
            'data': top_dir + '/data/',  # For statistics/read count data after mapping.
            'STAR index': '/oak/stanford/groups/khavari/users/dfporter/genome/GRCh38.gencode.29/star_index',
            'STAR': 'STAR',
            'bowtie2': 'bowtie2',
            'bowtie2 repeats index': '/oak/stanford/groups/khavari/users/dfporter/easyCLIP-dev/RepEnrich2/hg38re',
            'repeats annotation': '/oak/stanford/groups/khavari/users/dfporter/easyCLIP-dev/RepEnrich2/collapsed_by_family_hg38_repeatmasker_clean.txt',
            'STAR repeats index': '/oak/stanford/groups/khavari/users/dfporter/genome/repeats/star_repeats',
        }
        return file_paths

    def read_scheme(self, fname: Union[str, None] = None) -> sameRiver.scheme.scheme:
        
        if fname is not None:
            self.file_paths['scheme'] = fname
        self.scheme = sameRiver.scheme.scheme(self.file_paths['scheme'])
        
        return self.scheme

    def set_default_paths_if_needed(self):

        self.file_paths.setdefault('bowtie2', 'bowtie2')
        self.file_paths.setdefault(
            'bowtie2 repeats index', 
            '/oak/stanford/groups/khavari/users/dfporter/easyCLIP-dev/RepEnrich2/hg38re')
        self.file_paths.setdefault(
            'repeats annotation',
            '/oak/stanford/groups/khavari/users/dfporter/easyCLIP-dev/RepEnrich2/collapsed_by_family_hg38_repeatmasker_clean.txt')
        self.file_paths.setdefault('R1_fastq', self.file_paths['fastq'] + '/raw/R1.fastq')
        self.file_paths.setdefault('R2_fastq', self.file_paths['fastq'] + '/raw/R2.fastq')
        self.file_paths.setdefault('r1_split', self.file_paths['fastq'] + '/r1_split/')
        self.file_paths.setdefault('r1r2_split', self.file_paths['fastq'] + '/r1r2_split/')
        self.file_paths.setdefault('r2_clipped', self.file_paths['fastq'] + '/r2_clipped/')
        self.file_paths.setdefault('r1r2_clipped', self.file_paths['fastq'] + '/r1r2_clipped/')
        self.file_paths.setdefault('cutadapt', self.file_paths['fastq'] + '/ready_to_map/cutadapt/')
        
        self.file_paths['R1_fastq'] = self.use_gz_version_of_file_if_present(
            self.file_paths['R1_fastq'])
        self.file_paths['R2_fastq'] = self.use_gz_version_of_file_if_present(
            self.file_paths['R2_fastq'])
        
        auto_paths = self.autogenerate_paths(os.path.dirname(self.file_paths['scheme']))
        for path_name, path in auto_paths.items():
            if path_name not in self.file_paths:
                self.file_paths[path_name] = path

    def line_numbers(self) -> None:
        return self.read_stats()

    def read_stats(self) -> None:
        by_basename = collections.defaultdict(dict)

        print("Getting line numbers in {}...".format(self.file_paths['beds']))
        
        for fname in glob.glob('{}/*.bed'.format(self.file_paths['beds'])):
            cmd = f'wc -l {fname}'
            line_num = int(subprocess.check_output(cmd.split(' ')).split()[0])
#            print(out)
#            line_num = int(out.split(' ')[-2])
#            print(line_num)
            basename = os.path.splitext(os.path.basename(fname))[0]
            by_basename[basename]['bed line number'] = line_num

        r1r2_split = self.file_paths['r1r2_split']
        r1r2_clipped = self.file_paths['r1r2_clipped']
        cutadapt = self.file_paths['cutadapt']

        for (name, path) in [
            ('r1r2_split', r1r2_split), ('r1r2_clipped', r1r2_clipped),
            ('cutadapt', cutadapt)
            ]:
            print(f"Getting read numbers in {path}...")
            if os.path.exists(path):

                for fname in glob.glob(f'{path}/*.fastq'):
                    if '_R2_' in fname or '_R2.fastq' in fname:
                        continue

                    reads = int(subprocess.check_output(f'wc -l {fname}'.split(' ')).split()[0])/4
                    basename = os.path.splitext(os.path.basename(fname))[0]
                    by_basename[basename][name + ' reads'] = reads

        # Convert to dataframe.
        rows = []
        for basename, row in by_basename.items():
            rows.append(row)
            rows[-1]['basename'] = basename
        self.line_nums = pandas.DataFrame(rows)#, index='basename')

        # Save dataframe
        print(self.line_nums)
        self.line_nums.to_excel('line_numbers.xlsx')

    def determine_statistics(
        self, use_top_n: Union[int, None] = None,
        file_of_line_numbers: Union[str, None] = None, **kwargs):
        """Adds mapping statistics to the scheme dataframe from the fastq,
        bed, and bedgraph folders in self.file_paths. Can take a long time
        to calculate if not truncated with 'use_top_n'.
        """

        if 'use_top_n_reads' in kwargs and (use_top_n is None):
            use_top_n=kwargs['use_top_n_reads']

        if 'use_first_n_lines' in kwargs and (use_top_n is None):
            use_top_n=kwargs['use_first_n_lines']

        not hasattr(self, 'scheme') and self.read_scheme()
            
        stats = sameRiver.bedAndFastaStats.bedAndFastaStats(
            scheme_object=self.scheme,
            repeat_chr_names_file=Path(self.file_paths['sams'], "data", "names_of_repeats.txt"))

        skip = """
        if 'r1r2_split' in self.file_paths and (
            os.path.exists(self.file_paths['r1r2_split'])):
            
            stats.fastq_numbers(
                self.file_paths['r1r2_split'], use_top_n_reads=use_top_n,
                col_name='r1r2_split Reads')

        if 'r1r2_clipped_fastq' in self.file_paths and (
            os.path.exists(self.file_paths['r1r2_clipped_fastq'])):
            
            stats.fastq_numbers(
                self.file_paths['r1r2_clipped_fastq'], use_top_n_reads=use_top_n,
                col_name='r1r2_clipped_fastq Reads')
        """
        if 'fastq' in self.file_paths and (
            os.path.exists(self.file_paths['fastq'] + '/r1r2_split/')):
        
            stats.fastq_numbers(
                self.file_paths['fastq'] + '/r1r2_split/', use_top_n_reads=use_top_n,
                col_name='r1r2_split Reads')

        if 'fastq' in self.file_paths and (
            os.path.exists(self.file_paths['fastq'] + '/r1r2_clipped/')):
        
            stats.fastq_numbers(
                self.file_paths['fastq'] + '/r1r2_clipped/', use_top_n_reads=use_top_n,
                col_name='r1r2_clipped Reads')

        if 'fastq' in self.file_paths and (
            os.path.exists(self.file_paths['fastq'] + '/cutadapt/split/')):
        
            stats.fastq_numbers(
                self.file_paths['fastq'] + '/r1r2_clipped/', use_top_n_reads=use_top_n,
                col_name='r1r2_clipped Reads')

        print("Determining mapping numbers from bedgraph directory...")
        print(f"use_first_n_lines={use_top_n}")
        
        if 'wigs' not in self.file_paths:
            self.file_paths['wigs'] = self.file_paths['beds'] + '/read_start/'

        stats.mapping_numbers_from_bedgraph_dir(
                self.file_paths['wigs'], use_first_n_lines=use_top_n)

        print("Determining mapping numbers from bed directory...")
        stats.mapping_numbers_from_bed_dir_by_chrom(
                self.file_paths['beds'], use_top_n_reads=use_top_n, )

        self.scheme = stats.scheme

        outf = os.path.dirname(self.file_paths['scheme']) + '/stats_in_scheme.xlsx'
        self.scheme.scheme_df.to_excel(outf)

    def make_signal_data_file(self, clobber=True, serialize=True):
        """ Make data/bed_x.data object of sequencing coverage over chromosome locations.
        Will serialize regardless of serialize argument."""
        
        data = self.file_paths['data']
        os.makedirs(data, exist_ok=True)

        if serialize:
            output_data_filename = f"{data}/bed_{self.name}.data"
            if (not clobber) and os.path.exists(output_data_filename):
                return
        
        if 'wigs' not in self.file_paths:
            self.file_paths['wigs'] = self.file_paths['beds'] + '/read_start/'

        bedMaker = sameRiver.makeBedDataFile.bedDataFileMaker()
        
        # Originally, this function call was to dump to the .data file and returns a large object,
        # which we discarded. HTSeq.GenomicArray objects that are not empty will not serialize
        # currently on MacOS Mojave.
        """This fails with a seg fault on Mac in a Jupyter notebook (python3.8):
        _ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
        _ga[HTSeq.GenomicInterval('1', 0, 100, '+')] += 1
        with open('./data/bed_m0.data', 'wb') as f:
            dill.dump(_ga, f)
        """
        print("Making bedgraphs objects from the folder {}, which has {} plus strand files.".format(
            self.file_paths['wigs'], len(glob.glob(self.file_paths['wigs'] + '/*_+.wig')))
             )
        
        # Returns the serializable form.
        self.beds = bedMaker.make_bed_data_file(
#            bed_folder=self.file_paths['beds'],
            bedgraphs_folder=self.file_paths['wigs'],
            output_data_filename=f'{data}/bed_{self.name}.data',
            use='read start',
            collapse_reads=False)
        
        return self.beds

    def make_scheme_signal_RNA_data_files(
        self,
        rna_data_object: Union[None, sameRiver.set_of_named_mRNAs.set_of_named_mRNAs] = None,
        no_clobber: bool = False, serialize: bool = True):
        
        data = self.file_paths['data']

        if no_clobber and os.path.exists(f'{data}/signal_rnas_{self.name}.data'):
            print(f"scheme_signal_rna data file exists: {data}/signal_rnas_{self.name}.data." + \
                  " Not overwriting because no_clobber is on.")
            return
        
        if rna_data_object is None:
            print(f"Loading RNA data file {data}/rna.data'...", end='')
            rna_data_object = pickle.load(open(f'{data}/rna.data', 'rb'))
            #rna_data_object.find_introns()
            print("Finished loading RNA data file.")

        print("ENST00000383925 in rna_data_object?")
        if 'ENST00000383925' in rna_data_object.mRNAs:
            print('>yes')
        else:
            print('>no')

        print("RNU1-1 in rna_data_object?")
        if 'RNU1-1' in rna_data_object.mRNAs:
            print('>yes')
        else:
            print('>no')

        # Always returns the serializable form.
        self.beds = self.make_signal_data_file(clobber=not(no_clobber), serialize=serialize)

        # If serialize:
        print(f"Loading bed data file {data}/bed_{self.name}.data")
        self.beds = dill.load(open(f'{data}/bed_{self.name}.data', 'rb'))
        print("Finished loading bed data file.")
        
        self.beds.recover_bedgraph_objects()

        ssr = sameRiver.scheme_signal_RNAs.scheme_signal_RNAs(
            set_of_bedgraphs=self.beds, set_of_RNAs=rna_data_object)
        
        ssr.add_scheme(self.file_paths['scheme'])
        ssr.make_genomic_array(no_clobber=no_clobber)
        ssr.assign_to_genes(verbose=True, no_clobber=no_clobber)
        
        if 'counts' not in self.file_paths:
            self.file_paths['counts'] = 'counts.txt'

        ssr.write_counts_file(fname=self.file_paths['counts'], report_introns_separately=False)
        
        if serialize:
            with open(f'{data}/signal_rnas_{self.name}.data', 'wb') as f:
                pickle.dump(ssr, file=f)
    
    def duplication_rate(self) -> str:
        
        if not hasattr(self, 'scheme'):
            self.read_scheme()
        #beds = pickle.load(open('./data/bed_{0}.data'.format(self.name), 'rb'))
        
        #beds.duplication_rate()
        
        stats = sameRiver.bedAndFastaStats.bedAndFastaStats(
            scheme_object=self.scheme)

        reports = stats.duplication_rate(self.file_paths['beds'])
        
        if type(reports) != type(''):
            return 'No report\n {}'.format(reports)
        
        return reports
    
    def make_bedgraphs_from_bed(self,  use_first_n_lines=False, verbose=False, use='read start', **kwargs):
        self.bed_to_wig(use_first_n_lines=use_first_n_lines, verbose=verbose, use=use, **kwargs)
        
    def bed_to_wig(
        self, use_first_n_lines=False, verbose=False, use='read start',
        clobber=True, **kwargs):
        
        print(f"Making bedgraph files from bed files for {self.name}.")
        
        if verbose:
            print("beds: {}\nbedgraphs: {}".format(
                self.file_paths['beds'], self.file_paths['wigs']))
        
        if not clobber:
            have_all = True
            for bedfname in glob.glob(self.file_paths['beds'] + '/*.bed'):
                out_fname = self.file_paths['wigs'] + '/' + os.path.basename(bedfname).split('.')[0]
                out_fname = out_fname.replace('_+', '')
                out_fname = out_fname.replace('_-', '')
                if not(os.path.exists(out_fname + '_+.wig') and os.path.exists(
                    out_fname + '_-.wig')):
                    have_all = False
                    
            if have_all:
                return
            
        sbedgraphs = sameRiver.bedgraphs.set_of_bedgraphs(
            use_first_n_lines=use_first_n_lines, use=use,
            bed_folder=self.file_paths['beds'], clobber=clobber,
            **kwargs
            )

        if verbose:
            print("Writing bedgraphs...")
            
        sbedgraphs.write_bedgraphs(top_dir=self.file_paths['wigs'], verbose=verbose)

    def annotate_counts_file(self):
        
        ann_counts_fname = os.path.dirname(self.file_paths['counts']) + \
            '/ann_' + os.path.basename(self.file_paths['counts'])
        
        print(f"Annotating counts file and outputing to {ann_counts_fname}.")

        self.counts = sameRiver.countsO.countsO(
            filename=self.file_paths['counts'],
            scheme=self.read_scheme(),
            index_col=0,
            indexing_col_name='gene_name')

        self.counts.edit()

        self.counts.raw_counts_df.to_csv(ann_counts_fname, sep='\t')
        #pandas.DataFrame.from_dict(self.counts.raw_counts, orient='index').to_csv(
        #    ann_counts_fname, sep='\t')

#########################################################
# The methods below are all preprocessing scripts.
# Namely, going: raw fastq files -> split barcode fastq 
# -> clipped adapters (ours) -> clipped adapters (cutadapt)
# -> map to repeats -> map to genome -> split sam by barcode
# -> collapse duplicates in sam -> convert to bed and bedgraph
#########################################################
    def use_gz_version_of_file_if_present(self, fname):
        if not os.path.exists(fname) and os.path.exists(fname + '.gz'):
            print(f"Using {fname + '.gz'}...")
            return fname + '.gz'
        return fname

    def erase_fastq_files_in_directory(self, dirname):

        print(f"Deleting *.fastq files in {dirname}...")
        for fname in glob.glob(f'{dirname}/*.fastq'):
            os.remove(fname)


    def split_by_barcode(self):

        # Erase existing files in fastq intermediate directories.
        self.erase_fastq_files_in_directory(self.file_paths['r1_split'])
        self.erase_fastq_files_in_directory(self.file_paths['r1r2_split'])
        self.erase_fastq_files_in_directory(self.file_paths['r1r2_split'] + '/*/')
        self.erase_fastq_files_in_directory(self.file_paths['r2_clipped'])
        self.erase_fastq_files_in_directory(self.file_paths['r1r2_clipped'])
        self.erase_fastq_files_in_directory(self.file_paths['fastq'] + '/ready_to_map/')
        self.erase_fastq_files_in_directory(self.file_paths['fastq'] + '/ready_to_map/' + '/*/')

        if 'cutadapt' not in self.file_paths:
            self.file_paths['cutadapt'] =  self.file_paths['fastq'] + '/ready_to_map/cutadapt/'
        self.erase_fastq_files_in_directory(self.file_paths['cutadapt'])
        self.erase_fastq_files_in_directory(self.file_paths['cutadapt'] + '/*/')
        
        # No reads are discarded by either splitting function.
        sameRiver.split_r1.split_a_fastq_by_R1(
            self.file_paths['R1_fastq'],
            read2=self.file_paths['R2_fastq'],
            output_dir=self.file_paths['r1_split'],
            scheme_fname=self.file_paths['scheme'],
            )

        sameRiver.split_r2.split_a_dir_by_R2(
            self.file_paths['r1_split'],
            self.file_paths['r1r2_split'],
            self.file_paths['scheme'],
            )

        self.convert_split_r1r2_folder_to_long_filenames()

    def convert_split_r1r2_folder_to_long_filenames(self):

        cmds = []
        for dirname, subdirs, files in os.walk(self.file_paths['r1r2_split']):
            print(dirname)

            for fname in files:

                if len(fname.split('_')) < 2:
                    continue

                print(fname)
                p6 = fname.split('_')[0]
                p3 = fname.split('_')[-1].split('.fastq')[0]

                if (p6, p3) not in self.scheme.p6p3_to_long_filename_r1:
                    print(f"Skipping {fname} as it wasn't in the scheme. Looked for {p6}, {p3}.")
                    continue

                if '_R2_' in fname:
                    new_fname = self.scheme.p6p3_to_long_filename_r2[(p6, p3)]
                else:
                    new_fname = self.scheme.p6p3_to_long_filename_r1[(p6, p3)]

                cmds.append(
                    f"mv {dirname}/{fname} {self.file_paths['r1r2_split']}/{new_fname}")

        [os.system(cmd) for cmd in cmds]


    def preprocess_split_barcodes_folder(self, dirname: Union[str, None] = None):
        """
        Requires self.scheme, self.file_paths['fastq']
        """
        if dirname is None:
            dirname = self.file_paths['fastq'] + '/r1r2_split/'
        
        r2_clipped_dir = self.file_paths['r2_clipped']
        if not os.path.exists(self.file_paths['r2_clipped']):  
            os.system('mkdir ' + self.file_paths['r2_clipped'])
            
        r1r2_clipped_dir = self.file_paths['r1r2_clipped']
        
        if not os.path.exists(r1r2_clipped_dir):  
            os.system('mkdir ' + r1r2_clipped_dir)
            
        def get_p6p3_r1_fastq(_r2_fastq):

            basename = os.path.basename(_r2_fastq)
            if basename in self.scheme.long_fname_to_info:
                file_info = self.scheme.long_fname_to_info[basename]
            else:
                print(f"Not found in scheme: {_r2_fastq}")
                print(f"scheme values: {self.scheme.long_fname_to_info}")
            
            p6 = file_info['P6_BC']
            p3 = file_info['P3_BC']
            r1_fastq = os.path.dirname(_r2_fastq) + '/' + file_info['long_fname_R1']
            return p6, p3, r1_fastq
        
        for r2_fastq in glob.glob(dirname + '/*_R2_*.fastq'):

            (p6, p3, r1_fastq) = get_p6p3_r1_fastq(r2_fastq)
            
            self.log['clip_by_r2'] = sameRiver.clip_adapters.clip_by_r2(
                r1_fastq, r2_fastq, r2_clipped_dir)

        for r2_fastq in glob.glob(r2_clipped_dir + '*_R2_*.fastq'):
            (p6, p3, r1_fastq) = get_p6p3_r1_fastq(r2_fastq)
            
            self.log['clip_by_r1'] = sameRiver.clip_adapters.clip_by_r1(
                r1_fastq, r2_fastq, r1r2_clipped_dir, bc_l3=p3)

        self.write_to_log()

        # The above two clipping calls to clip_adapters only:
        # 1. Remove the L3 barcode and anything after from R1.
        # 2. Remove the L5 random hexamer and anything after from R2.
        # 3. Remove read pairs where one of the reads ends up too short.
        # At this point the L3 barcode is still in R2 and the L5 barcode/hexamer is still in R1.
        
        # Move barcodes to the read name.
        # No reads are discarded in make_mappable.
        (r1_fnames, r2_fnames) = sameRiver.combine_fastqs_for_mapping.make_mappable(
            r1r2_clipped_dir, self.file_paths['fastq'] + '/ready_to_map/', self.file_paths['scheme'])
        print("exp returned from make_mappable(): ", (r1_fnames, r2_fnames))

        if 'cutadapt' not in self.file_paths:
            self.file_paths['cutadapt'] =  self.file_paths['fastq'] + '/ready_to_map/cutadapt/'

        if not os.path.exists(self.file_paths['cutadapt']):
            os.system('mkdir ' + self.file_paths['cutadapt'])
        
        # Remove adapters using a system call to cutadapt.
        # Cutadapt removes reads < 17 len.
        cut_fnames = []
        for r1_fname, r2_fname in zip(r1_fnames, r2_fnames):
            cut_fnames.append(sameRiver.clip_adapters.clip_a_file_pair_using_cutadapt(
                r1_fname, r2_fname, self.file_paths['cutadapt'] + '/split/'
                #self.file_paths['fastq'] + '/ready_to_map/concatenated_R1_reads_for_mapping.fastq',
                #self.file_paths['fastq'] + '/ready_to_map/concatenated_R2_reads_for_mapping.fastq',
                ))
        r1_fnames = [x[0] for x in cut_fnames]
        r2_fnames = [x[1] for x in cut_fnames]
        print("clip_a_file_pair_using_cutadapt() returned: ", (r1_fnames, r2_fnames))

        sameRiver.combine_fastqs_for_mapping.concatenate_fastqs(
            r1_fnames=r1_fnames, r2_fnames=r2_fnames,
            concat_file_r1=self.file_paths['cutadapt'] + '/concatenated_R1_reads_for_mapping.fastq',
            concat_file_r2=self.file_paths['cutadapt'] + '/concatenated_R2_reads_for_mapping.fastq',
        )

        try:
            self.log['reads_pre_cutadapt'] = int(subprocess.check_output(
                f"wc -l {self.file_paths['fastq'] + '/ready_to_map/concatenated_R1_reads_for_mapping.fastq'}".split(' ')).split()[0])/4    
            self.log['reads_post_cutadapt'] = int(subprocess.check_output(
                f"wc -l {self.file_paths['cutadapt'] + '/concatenated_R1_reads_for_mapping.fastq'}".split(' ')).split()[0])/4
        except:  # Not really important.
            pass

        self.write_to_log()

