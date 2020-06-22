import os, HTSeq, importlib, collections, glob, re, subprocess
import sameRiver
import sameRiver.scheme
import sameRiver.bedgraphs
from typing import List, Mapping, Union
from pprint import pprint

importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.bedgraphs)

class bedStats:
    
    def __init__(self, scheme, bed_dir, bedgraph=False, repeat_chr_names=set()):
        self.scheme = scheme
        self.bed_dir = bed_dir
        self.is_bedgraph = bedgraph
        self.repeat_chr_names = repeat_chr_names
    
    def mapping_numbers_from_bedgraph_dir(
        self, bedgraphs_folder=None, bedgraphs_object=None, use_first_n_lines=None, **kwargs):
        
        print(f"Counting mapped reads from bedgraph dir. kwargs={kwargs}.")
        print(f" use_first_n_lines={use_first_n_lines}")

        kwargs.setdefault('verbose', False)
            
        if bedgraphs_folder is not None and bedgraphs_object is not None:
            raise ValueError(
                "Pass either a sameRiver.bedgraphs.set_of_bedgraphs object or a folder, but not both.")
        
        # Get a set_of_bedgraphs object from a passed folder.
        if bedgraphs_folder is not None:
            wig = sameRiver.bedgraphs.set_of_bedgraphs(
                bedgraphs_folder=bedgraphs_folder, use_first_n_lines=use_first_n_lines)

        # Use the given set_of_bedgraphs object is passed one.
        if bedgraphs_object is not None:
            wig = bedgraphs_object
        
        # Calculate the number of reads for each bedgraph.
        wig.total_area_under_curve(verbose=kwargs['verbose'])
        
        self.counts = {}
        for bedgraph_fname, auc in wig.aucs.items():
            p6p3 = os.path.basename(bedgraph_fname).split('.wig')[0].rstrip('_+').rstrip('_-')
            self.counts[p6p3] = auc
            
        total = sum(self.counts.values())
        self.frac_counts = dict([(k, 100 * v/total) for k, v in self.counts.items()])
        
        counts_str = {k:f"{float(v):,}" for k,v in self.counts.items()}
        frac_counts_str = {k:f"{float(v):.5}" for k,v in self.counts.items()}

        self.scheme.fill_in_info_column_by_filename(
            'Mapped reads (bedgraph)', counts_str)
        self.scheme.fill_in_info_column_by_filename(
            'Mapped reads (bedgraph, fraction total)', frac_counts_str)
        
        print(self.counts)
        
        return self.scheme

    def mapping_numbers_from_bed_dir(
        self, use_top_n_reads=None):
        """Fill in self.scheme with usable read counts."""

        self.counts = {}

        for bedfname in glob.glob(self.bed_dir + '/*.bed'):

            p6p3 = os.path.basename(bedfname).split('.bed')[0]
            self.counts[p6p3] = self.count_reads_bed(
                bedfname, use_top_n_reads=use_top_n_reads, repeat_chr_names=self.repeat_chr_names)

        total = sum(self.counts.values())
        self.frac_counts = dict([(k, 100 * v/total) for k, v in self.counts.items()])

        print(self.counts)
        print(self.frac_counts)

        counts_str = {k:f"{float(v):,}" for k,v in self.counts.items()}
        frac_counts_str = {k:f"{float(v):.5}" for k,v in self.counts.items()}
        
        self.scheme.fill_in_info_column_by_filename('Mapped reads', counts_str)
        self.scheme.fill_in_info_column_by_filename(
            'Mapped reads (fraction total)', frac_counts_strs)

        return self.scheme
        
    def _mapping_numbers_from_bed_dir_by_chrom(self, simplify=True):
        """Fill in self.scheme with usable read counts, keeping chromosome info."""

        new_info_cols = collections.defaultdict(dict)

        for bedfname in glob.glob(self.bed_dir + '/*.bed'):
            p6p3 = os.path.basename(bedfname).split('.bed')[0]

            # Get {chromosome_category: read_count}
            chrom_to_counts = self.count_reads_bed_by_chrom(
                bedfname, simplify=simplify, repeat_chr_names=self.repeat_chr_names)

            print(chrom_to_counts)

            total_for_this_bed = sum(chrom_to_counts.values())

            frac_each_chrom_for_this_bed = dict([
                (k, 100 * v/total_for_this_bed) for k, v in chrom_to_counts.items()])

            new_info_cols['Mapped reads'][p6p3] = total_for_this_bed
            for chr_name, count_at_chr in chrom_to_counts.items():
                new_info_cols[f'Mapped reads at {chr_name}'][p6p3] = \
                    f'{count_at_chr:,} ({frac_each_chrom_for_this_bed[chr_name]:.3}%)'

        self.counts = new_info_cols['Mapped reads']
        total = sum(new_info_cols['Mapped reads'].values())
        self.frac_counts = dict([(k, 100 * v/total) for k, v in self.counts.items()])
        self.scheme.fill_in_info_column_by_filename('Mapped reads', self.counts)
        self.scheme.fill_in_info_column_by_filename('Mapped reads (fraction total)', self.frac_counts)

        for col_name, data_dict in new_info_cols.items():
            print(f"\n\n###\nFilling in {col_name}: {data_dict}")
            self.scheme.fill_in_info_column_by_filename(col_name, data_dict)

        print(self.scheme.scheme_df)
        print(self.scheme.scheme_df.columns)
        print(self.scheme.scheme_df.iloc[0])
        return self.scheme

    def duplication_rate(
        self, bed_dir, use_top_n_reads=None):
        
        reports = "-" * 10 + ' \nDuplication rates from bed files:'
        
        for bedfname in glob.glob(bed_dir + '/*.bed'):
            
            p6p3 = os.path.basename(bedfname).split('.bed')[0]
            reports += f"\n***\nDuplication rate for {bedfname}:\n"
            report = self.count_unique_reads_bed(bedfname, use_top_n_reads=use_top_n_reads)
            reports += report + '\n'
            
        return reports

    @classmethod
    def count_reads_bed(
        cls, bedfname, use_top_n_reads=None, by_chrom=False, repeat_chr_names=set()) -> int:
        
        if by_chrom:
            return cls.count_reads_bed_by_chrom(bedfname, repeat_chr_names=repeat_chr_names)

        if use_top_n_reads is None:
            return int(subprocess.check_output(f'wc -l {bedfname}'.split(' ')).split()[0])

        n = 0
        with open(bedfname) as fh:
            for li in fh:
                n += 1
                if n >= use_top_n_reads:
                    break

        return n

    @staticmethod
    def count_reads_bed_by_chrom(
        bedfname: str, simplify: bool = False, repeat_chr_names=set()) -> Mapping[str, int]:

        # Count the instances of reads at each chromosome (column 1).
        chrs = subprocess.check_output(f'cut -f1 {bedfname}'.split(' '))
        chrs = collections.Counter(chrs.decode().split('\n'))

        rRNA = ['LSU-rRNA_Hsa', 'SSU-rRNA_Hsa', '5S']

        if '' in chrs:
            del chrs['']

        print(f"Before simpl: {chrs}")
        
        if simplify:
            simple_dict = collections.defaultdict(int)
            keep = ['chrY', 'chrX', 'repeats'] # 'chrM', 
            odds = set()
            for a_chr, count_at_chr in chrs.items():
                if a_chr in keep:
                    simple_dict[a_chr] += count_at_chr
                elif a_chr in rRNA:
                    simple_dict['rRNA'] += count_at_chr
                elif a_chr in repeat_chr_names:
                    simple_dict['repeats'] += count_at_chr
                elif a_chr[:3] == 'chr':
                    simple_dict['autosomes'] += count_at_chr
                else:
                    simple_dict['Other genomic sequence'] += count_at_chr
                    odds.add(a_chr)
            print(f"Odd chromosomes {odds}.")
            pprint(f"count_reads_bed_by_chrom(): {simple_dict}")

            return simple_dict

        return chrs
    
    @staticmethod
    def count_unique_reads_bed(bedfname, use_top_n_reads=None):
        """Format:
        ==> /Users/dfporter/pma/miseq/Runs/180815/beds/TGTTGG_NNN.bed <==
        chr12	20551477	20551540	TGCTCGG	1	+
        s[0] chr
        s[1], s[2] start, end
        s[3] 7N random barcode on L5
        s[4] count
        s[5] strand
        """
        n = 0

        locations = collections.defaultdict(list)
        with open(bedfname) as fh:

            if use_top_n_reads is None:
                for li in fh:
                    s = li.rstrip('\n').split('\t')
                    if s[5] == '+':
                        locations[(s[0], s[1], s[5])].append(s[3])
                    else:
                        locations[(s[0], s[2], s[5])].append(s[3])
                    n += 1
                    

            else:
                for li in fh:
                    s = li.rstrip('\n').split('\t')
                    if s[5] == '+':
                        locations[(s[0], s[1], s[5])].append(s[3])
                    else:
                        locations[(s[0], s[2], s[5])].append(s[3])
                    n += 1
                    if n >= use_top_n_reads:
                        break

        unique = {}
        n_unique = 0
        for iv in locations:
            unique[iv] = set(locations[iv])
            n_unique += len(unique[iv])
            
        histogram_raw = collections.defaultdict(int)
        histogram_unique = collections.defaultdict(int)
        for iv in locations:
            histogram_raw[len(locations[iv])] += 1
            histogram_unique[len(unique[iv])] += 1
            
        report = """Total reads {:,} Unique reads {:,} % unique {:,}\n""".format(
            n, n_unique, n_unique/max([n, 1]))
        report += """Histograms of raw read counts/location:\n"""
        for n_reads in sorted(histogram_raw.keys(), key=lambda x: x):
            report += "{}: {}\n".format(n_reads, histogram_raw[n_reads])  # n_reads, frequency
        report += """Histograms of unique read counts/location:\n"""
        for n_reads, frequency in histogram_unique.items():
            report += "{}: {}\n".format(n_reads, frequency)
        
        return report
    
    
class bedAndFastaStats:
    
    def __init__(self, scheme_fname=None, scheme_object=None, repeat_chr_names_file=None):
        
        if (scheme_fname is None) and (scheme_object is None):
            raise ValueError("Must give a scheme filename or object.")
        if (scheme_fname is not None) and (scheme_object is not None):
            raise ValueError("Must give either a scheme filename or object, and not both.")
        
        if scheme_fname is not None:
            self.scheme = sameRiver.scheme.scheme()
            self.scheme.read_scheme(scheme_fname)

        if scheme_object is not None:
            self.scheme = scheme_object

        # Read in the chromosome names from a text file if passed one.
        if (repeat_chr_names_file is not None) and os.path.exists(repeat_chr_names_file):
            with open(repeat_chr_names_file) as fh:
                self.repeat_chr_names = set([x.rstrip('\n') for x in fh.readlines() if x!=''])
        else:
            self.repeat_chr_names = {'repeats'}
            
    def mapping_numbers_from_bed_dir(
        self, bed_dir, use_top_n_reads=None):
        
        bs = bedStats(self.scheme, bed_dir)
        self.scheme = bs.mapping_numbers_from_bed_dir(use_top_n_reads=use_top_n_reads)
    
    def mapping_numbers_from_bed_dir_by_chrom(
        self, bed_dir, **kwargs):

        bs = bedStats(self.scheme, bed_dir, repeat_chr_names=self.repeat_chr_names)
        self.scheme = bs._mapping_numbers_from_bed_dir_by_chrom()

    def mapping_numbers_from_bedgraph_dir(
        self, bedgraphs_folder, use_first_n_lines=False):
        
        bs = bedStats(self.scheme, bedgraphs_folder)
        self.scheme = bs.mapping_numbers_from_bedgraph_dir(
            bedgraphs_folder=bedgraphs_folder, use_first_n_lines=use_first_n_lines)
    
    def fastq_numbers(
        self, initial_fastq, star_text_file=None, col_name='Reads', use_top_n_reads=None):

        fs = fastqStats('')
        self.scheme = fs.fastq_numbers(
            self.scheme, initial_fastq, star_text_file=star_text_file,
            col_name=col_name, use_top_n_reads=use_top_n_reads)
        
    def duplication_rate(self, bed_dir, use_top_n_reads=None):
        
        bs = bedStats(self.scheme, bed_dir)
        return bs.duplication_rate(bed_dir, use_top_n_reads=use_top_n_reads)
        

class fastqStats:
    
    def __init__(self, fastq):
        self.fastq = fastq
    
    @staticmethod
    def count_r1_reads(fastq, use_top_n_reads=None):
        reader = HTSeq.FastqReader(fastq)
        r1_counts = collections.defaultdict(int)

        if use_top_n_reads is None:
            for n, read in enumerate(reader):
                r1_counts[read.seq[:6]] += 1
        else:
            for n, read in enumerate(reader, start=1):
                r1_counts[read.seq[:6]] += 1
                if n >= use_top_n_reads:
                    break

        r1_counts_strs = {}
        for k in list(r1_counts.keys())[:]:
            r1_counts_strs[str(k, 'utf-8')] = r1_counts[k]

        return r1_counts_strs

    @staticmethod
    def count_reads_fastq(fastq, use_top_n_reads=None):

        if use_top_n_reads is None:
            return int(subprocess.check_output(f'wc -l {fastq}'.split(' ')).split()[0])/4

        else:

            lines = 0
            _use_top_n_reads = 4 * use_top_n_reads

            with open(fastq) as f:
                for n, li in enumerate(f):
                    lines += 1
                    if n >= _use_top_n_reads:
                        break

            return int(lines/4)
    
    @classmethod
    def count_r1_reads_with_scheme(cls, scheme, fastq, use_top_n_reads=None):

        r1_counts = cls.count_r1_reads(fastq, use_top_n_reads=use_top_n_reads)
        not_in_scheme = {}
        in_scheme = {}

        for bc, count in r1_counts.items():

            if bc in scheme.p6_to_gene:
                in_scheme[bc] = count            
            else:
                not_in_scheme[bc] = count

        scheme.fill_in_info_column_by_p6('R1 reads with P6 barcode', in_scheme)
        
        total_not_in_scheme = sum(not_in_scheme.values())
        total_in_scheme = sum(in_scheme.values())
        total = total_not_in_scheme + total_in_scheme

        fractions = dict([(k, 100 * v/total_in_scheme) for k, v in in_scheme.items()])
        scheme.fill_in_info_column_by_p6('R1 reads with P6 barcode (% recognized barcodes)', fractions)

        print(scheme.scheme_df)

        li = "Mapping statistics for raw fastq file, R1:\n"
        li += "Total reads: {0}\n".format(total)
        li += "Has expected P6 barcode {0} {1:.5}%, Unexpected {2} {3:.5}%\n".format(
            total_in_scheme, 100 * total_in_scheme/total, total_not_in_scheme, 100 * total_not_in_scheme/total)
        li += "=" * 50 + '\n'
        li += "P6 barcodes in scheme:\n"
        li += "P6_BC\tReads (R1)\t% all reads\t% reads with recognized barcodes\tProtein\n"

        sorted_bc = sorted(in_scheme.keys(), key=lambda x: in_scheme[x])[::-1]
        for bc in sorted_bc:
            val = in_scheme[bc]
            li += "{0}\t{1}\t{2:.6}%\t{3:.6}%\t{4}\n".format(bc, val, 100* val/total, 100*val/total_in_scheme,
                                                             ', '.join(scheme.p6_to_gene[bc]))
        li += "=" * 50 + '\n'
        li += "Top 10 P6 barcodes not in scheme:\n"
        li += "P6_BC\tReads (R1)\t% all reads\n"

        sorted_bc = sorted(not_in_scheme.keys(), key=lambda x: not_in_scheme[x])[::-1]
        for bc in sorted_bc[:10]:
            val = not_in_scheme[bc]
            li += "{0}\t{1}\t{2:.6}%\n".format(bc, val, 100* val/total)
        print(li)
    #print(not_in_scheme)


    def R1_fastqs(cls, dirname):
        for fname in glob.glob(dirname + '/*.fastq'):
            if re.search('_R2', os.path.basename(fname)):
                continue
            yield fname

    def fastq_numbers(
        cls, _scheme, initial_fastq, star_text_file=None, col_name='Reads', **kwargs):

        print(f"\n*\nDetermining fastq numbers ({initial_fastq})...")
        if 'use_top_n_reads' in kwargs:
            use_top_n_reads = kwargs['use_top_n_reads']
        else:
            use_top_n_reads = None
            
        bname = os.path.basename(initial_fastq)

        if not len(bname):
            print("fastq_numbers() passed a directory for fastq - assuming it's R1/R2 split reads.")
            p6p3_reads = {}

            for fname in cls.R1_fastqs(initial_fastq):
                p6p3 = os.path.splitext(os.path.basename(fname))[0]
                p6p3_reads[p6p3] = cls.count_reads_fastq(fname, use_top_n_reads=use_top_n_reads)

            if len(p6p3_reads) == 0:
                print("Looking in subfolders...")
                for fname in cls.R1_fastqs(initial_fastq + '/*/'):
                    p6p3 = os.path.splitext(os.path.basename(fname))[0]
                    p6p3_reads[p6p3] = cls.count_reads_fastq(fname, use_top_n_reads=use_top_n_reads)

            print(f"Read {len(p6p3_reads)} fastqs.")

            _scheme.fill_in_info_column_by_filename(col_name, p6p3_reads)

            total = sum(p6p3_reads.values())
            frac_counts = dict([(k, 100 * v/total) for k, v in p6p3_reads.items()])

            _scheme.fill_in_info_column_by_filename(col_name + " (% total)", frac_counts)

        if bname.split('.')[-1] == 'fastq':
            cls.count_r1_reads_with_scheme(_scheme, initial_fastq)#, use_top_n_reads=1E5)

        print('fastq_numers():\n', _scheme.scheme_df)

        return _scheme
