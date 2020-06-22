import collections, pandas, os, re, glob, sys, importlib, pickle, subprocess, time, json, pysam
from typing import List, Tuple, Union, Mapping
from pprint import pprint
from pathlib import Path
from argparse import Namespace

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# This RepEnrich2 has to be the modified python3.8+ one:
# https://github.com/dfporter/RepEnrich2/tree/py3
from RepEnrich2 import RepEnrich2
        
import sameRiver
import sameRiver.scheme
import sameRiver.makeBedDataFile
import sameRiver.scheme_signal_RNAs
import sameRiver.bedAndFastaStats
import sameRiver.bedgraphs
import sameRiver.clip_adapters
import sameRiver.combine_fastqs_for_mapping
import sameRiver.sam_to_bed_and_wig
import sameRiver.collapse_duplicates


importlib.reload(sameRiver.scheme_signal_RNAs)
importlib.reload(sameRiver.makeBedDataFile)
importlib.reload(sameRiver.bedAndFastaStats)
importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.bedgraphs)
importlib.reload(sameRiver.clip_adapters)
importlib.reload(sameRiver.combine_fastqs_for_mapping)
importlib.reload(sameRiver.sam_to_bed_and_wig)
importlib.reload(sameRiver.collapse_duplicates)

# Special repeats are those in the categories of ['rRNA', 'scRNA', 'snRNA', 'tRNA'].
special_repeats = [
    'tRNA-Arg-CGA', 'tRNA-Ser-TCG', 'HY4', 'tRNA-Met-i', 'tRNA-Gln-CAG', 'U1', 'tRNA-Val-GTA', 'tRNA-Ser-TCY', 'tRNA-Val-GTY', 'tRNA-Gln-CAA_', 'U13_', 'tRNA-Ser-TCA', 'tRNA-Pro-CCA', 'U5', 'tRNA-Arg-CGG',
    'tRNA-Arg-AGA', 'tRNA-Arg-CGA_', 
    'tRNA-Met_', 'tRNA-Leu-CTY', 'SSU-rRNA_Hsa', 'tRNA-Ile-ATT', 'tRNA-Ile-ATC', 'BC200', 'tRNA-Leu-TTG', 'LSU-rRNA_Hsa', 'tRNA-Gly-GGY', 'tRNA-Leu-CTA', 
    'tRNA-Trp-TGG', 'tRNA-Thr-ACG_', 'tRNA-His-CAY', 'tRNA-Ala-GCY_', 'CRP1', 'tRNA-Arg-AGG', 'tRNA-Gln-CAA', '5S', 'tRNA-Glu-GAG', 'LFSINE_Vert', 'U4', 'tRNA-Asn-AAT', 'tRNA-Thr-ACG', 'tRNA-Lys-AAA', 'U6', 
    'tRNA-Pro-CCG', 'tRNA-Ser-AGY', 'U13', 'tRNA-Arg-CGY_', 'tRNA-Ala-GCG', 'tRNA-Ala-GCA', 'tRNA-Thr-ACY_', 'U7', 'tRNA-Phe-TTY', 'tRNA-Leu-CTA_', 'tRNA-Asn-AAC', 'tRNA-Thr-ACA', 'tRNA-Lys-AAG', 'U14', 
    'tRNA-Pro-CCY', 'tRNA-Gly-GGA', 'tRNA-Ser-TCA_', 'tRNA-Met', 'tRNA-His-CAY_', 'tRNA-Glu-GAA', 'tRNA-Leu-CTG', 'HY1', 'tRNA-Glu-GAG_', 'U3', 'tRNA-Tyr-TAC', 'HY3', 'tRNA-Thr-ACY', 'tRNA-Ala-GCY', 'tRNA-Gly-GGG', 'U2', 'tRNA-Cys-TGY',
    'tRNA-Asp-GAY', 'U17', 'tRNA-Ile-ATA', 'tRNA-Val-GTG', 'tRNA-Tyr-TAT', 'tRNA-Leu-TTA', 'HY5', 'U8'
    #'tRNA-SeC(e)-TGA', 'tRNA-Leu-TTA(m)', 'tRNA-Ser-TCA(m)',
]
snRNA = ["U6", "U2", "U13_", "U7", "U3", "U1", "U5", "U13", "U4", "U8", "U17", "U14",]

class starCaller():

    def __init__(self):
        pass

    def star_cmd_to_genome(
        self, paths=None, read1_fname: str='', read2_fname: str=''):

        if paths is None:
            paths = self.file_paths

        cmd = paths['STAR']
        cmd += f" --genomeDir {paths['STAR index']}"
        cmd += ' --runThreadN 10'
        cmd += f' --readFilesIn {read1_fname} {read2_fname}'
        cmd += ' --outReadsUnmapped Fastx'
        # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
        # of a paired end read) reads in separate file(s).
        # output in separate fasta/fastq files, Unmapped.out.mate1/2
        cmd += f" --outFileNamePrefix {paths['sams']}/genome."

        return cmd

def fastaReader(fname: str):
    """Read a file with ONLY ONE fasta sequence and return a dict of {name: sequence}."""

    seq = []
    with open(fname) as f:
        for li in f:
            if li[0] == '>':
                continue
            else:
                seq.append(li.rstrip('\n'))
    return {os.path.splitext(os.path.basename(fname))[0]: ''.join(seq)}

class repeatsGenome():

    def __init__(self):
        pass

    def setup_genomes(
        self,
        repeats_fasta_directory='hg38re/',
        genomic_gtf='gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA',
        igv_output_directory='for_igv/',
        repeats_gtf='repeats_as_separate_chroms.gtf',
        combined_gtf='repeats_and_genome.gtf',
        ) -> None:

        # For the repeats, make the .fa and .fa.fai needed to build an IGV genome.
        self.make_files_for_igv_genome(
            repeats_fasta_directory=repeats_fasta_directory,
            output_directory=igv_output_directory)

        # Make a GTF for the repeats genome.
        self.make_gtf_of_repeats_genome(
            repeats_fasta_directory=repeats_fasta_directory,
            output_filename=repeats_gtf)

        # Simply concatenates the files.
        self.combine_genomic_and_repeats_gtf(
            repeats_gtf=repeats_gtf,
            genomic_gtf=genomic_gtf,
            combined_gtf=combined_gtf)

    def combine_genomic_and_repeats_gtf(
        self, repeats_gtf: str='repeats_as_separate_chroms.gtf',
        genomic_gtf: str='gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA',
        combined_gtf='repeats_and_genome.gtf',
        ) -> None:

        # Print help.
        if any(not os.path.exists(x) for x in [repeats_gtf, genomic_gtf]):
            print("""Input file not found. Arguments: repeats_fasta_directory: str='', output_filename: str=''""")

        # cat might not work everywhere, so we do this the slow way.
        with open(combined_gtf, 'w') as outf:

            # Write genomic gtf.
            with open(genomic_gtf) as f:
                for li in f:
                    outf.write(li)

            # Write repeats gtf.
            with open(repeats_gtf) as f:
                for li in f:
                    outf.write(li)
        
        print(f"Wrote the combined gtf file to {combined_gtf}.")

    def make_gtf_of_repeats_genome(
        self, repeats_fasta_directory: str='', output_filename: str='') -> None:

        # Print help.
        if any(x=='' for x in [repeats_fasta_directory, output_filename]):
            print("""Arguments: repeats_fasta_directory: str='', output_filename: str=''""")

        # Grab the fasta files from the input directory.
        fnames = [Path(x) for x in glob.glob(repeats_fasta_directory + '/*fa')]

        fastas = {}  # This could be more elegant.
        [fastas.update(fastaReader(fname)) for fname in fnames]

        # Make the output folder if needed.
        os.makedirs(os.path.dirname(os.path.realpath(output_filename)), exist_ok=True)

        # Write the GTF with gene/transcript/exon for each repeat.
        gtf_out = open(Path(output_filename), 'w')
        start = 1  # GTF is 1-based
        for name, seq in fastas.items():
            end = len(seq)
            li = f'{name}\tENSEMBL\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_support_level "NA"\n'
            li += f'{name}\tENSEMBL\ttranscript\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_type "repeat"; transcript_name "{name}; transcript_support_level "NA"; tag "basic";\n'
            li += f'{name}\tENSEMBL\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_type "repeat"; transcript_name "{name}; exon_number 1; exon_id "{name}.1"; transcript_support_level "NA"; tag "basic";\n'
            gtf_out.write(li)

        gtf_out.close()

    def make_files_for_igv_genome(
        self, repeats_fasta_directory: str='', output_directory: str='') -> None:
        """Write a .fa and .fa.fai file to load repeats into IGV.
        Builds the single .fa from a folder of fastas, namely the folder holding the
        Bowtie2 indexes."""

        # Print help.
        if any(x=='' for x in [repeats_fasta_directory, output_directory]):
            print("""Arguments: repeats_fasta_directory: str='', output_directory: str=''""")

        # Grab the fasta files from the input directory.
        fnames = [Path(x) for x in glob.glob(repeats_fasta_directory + '/*fa')]

        print(f"Loading {len(fnames)} fasta files from {repeats_fasta_directory}...")
        fastas = {}  # This could be more elegant.
        [fastas.update(fastaReader(fname)) for fname in fnames]

        # Make the output folder if needed.
        os.makedirs(output_directory, exist_ok=True)

        # Write the fastas to a single new output fasta.
        output_fasta = Path(output_directory, 'repeats_separate_chroms.fa')
        with open(output_fasta, 'w') as f:
            for name, seq in fastas.items():
                f.write(f'>{name}\n')
                [f.write(seq[i:i+100] + '\n') for i in range(0, len(seq), 100)]

        print(f"Wrote fasta to {output_fasta}.")

        # Make the fasta index.
        os.system(f'samtools faidx {output_fasta}')

    def igv_version_of_gtf(self, gtf_filename: str='', out_gtf_fname: str='') -> None:
        """Switches gene_id and gene_name tags in a gtf file so IGV will 
        label with the gene name.
        """

        # Print help.
        if any(x=='' for x in [gtf_filename, out_gtf_fname]):
            print("""Arguments: gtf_filename: str='', out_gtf_fname: str=''""")

        outfh = open(out_gtf_fname, 'w')

        with open(gtf_filename) as f:

            for li in f:
                gene_id = re.search('gene_id "([^"]+)"', li)#.group(1)
                transcript_id = re.search('transcript_id "([^"]+)"', li)#.group(1)
                gene_name = re.search('gene_name "([^"]+)"', li)#.group(1)
                transcript_name = re.search('transcript_name "([^"]+)"', li)

                if (gene_id is not None) and (gene_name is not None):
                    gene_id = gene_id.group(1)
                    gene_name = gene_name.group(1)
                    if transcript_name is not None:
                        transcript_name = transcript_name.group(1)
                        li = re.sub('transcript_id "[^"]+"', f'transcript_id "{transcript_name}"', li)
                    else:
                        gene_id_suff = gene_id.split('.')[-1]
                    #li = re.sub('gene_id "[^"]+"', f'gene_id "{gene_name}"', li)
                    #li = re.sub('gene_name "[^"]+"', f'gene_name "{gene_id}"', li)
                        li = re.sub('transcript_id "[^"]+"', f'transcript_id "{gene_name + gene_id_suff}"', li)

                outfh.write(li)

        outfh.close()

    def read_gtf(self, gtf_fname: str='') -> Mapping[str, List[int]]:
        """Return a dict of {gene_id: [start, end]} chromosome positions.
        Not used."""

        gene_positions = {}

        with open(gtf_fname) as f:
            for li in f:

                if li[0] == '#':
                    continue  # Skip comments.

                s = li.split('\t')  # Chr, ENSMBL/HAVANA, gene/exon/trxpt, start_pos, end_pos...

                if s[2] == 'gene':
                    gene_id = re.search('gene_id "([^"]+)";', s[-1]).group(1)
                    # GTF is 1-based.
                    gene_positions[gene_id] = [int(s[3]), int(s[4])]  # [start pos, end pos]

        return gene_positions

    def convert_sam_file_with_repeats_as_chroms_into_single_chrom(
        self, sam_fname:str='',
        gtf_fname:str='',
        out_sam_fname:str='',
        artificial_chom_fasta_index:str='',
        ) -> None:
        """Not used."""

        if any(x=='' for x in [sam_fname, gtf_fname, out_sam_fname, artificial_chom_fasta_index]):
            print("""Arguments: sam_fname:str='', gtf_fname:str='', out_sam_fname:str='', artificial_chom_fasta_index:str=''""")

        # Get the length of the artifical chromosome in bp.
        with open(artificial_chom_fasta_index) as f:
            _index = next(f).split('\t')  # Name, length in bp, offset...
            chrom_length = int(_index[1])

        gene_positions_in_artificial_chrom = self.read_gtf(gtf_fname)

        out_sam = open(out_sam_fname, 'w')
        out_sam.write(f'@SQ\tSN:repeats\tLN:{chrom_length}\n')

        with open(sam_fname) as f:
            for li in f:

                # Catch and write header.
                if li[0] == '@':

                    # Skip chromosome lines.
                    if li[:3] != '@SQ':
                        out_sam.write(li)
                    continue

                s = li.split('\t')  # Read name, flag, chrom, position...
                # SAM is 1-based, like GTF.
                gene_id = s[2]  # Chromosome is the gene id for repeats initially. 

                # If this is a repeat, edit and write. Don't write otherwise.
                if gene_id in gene_positions_in_artificial_chrom:
                    s[2] = 'repeats'
                    s[3] = str(gene_positions_in_artificial_chrom[gene_id][0] + int(3))

                    # We don't need to change anything for negative strand mapping.
                    out_sam.write('\t'.join(s))  # Write the new SAM line.

        out_sam.close()

    def make(
        self, repEnrichBowtieIndexFolder: str = 'AluOnly_RepEnrich_index/',
        repeats_artificial_chromosome_folder: str = None) -> None:
        """Synonym for make_artificial_chromosome."""

        self.make_artificial_chromosome(
            repEnrichBowtieIndexFolder=repEnrichBowtieIndexFolder,
            repeats_artificial_chromosome_folder=repeats_artificial_chromosome_folder)

    def make_artificial_chromosome(
        self, repEnrichBowtieIndexFolder: str = 'AluOnly_RepEnrich_index/',
        repeats_artificial_chromosome_folder: str = None) -> None:
        """
        Not using this! No good reason to use the artifical chromosome at present.

        Use the fasta files in the bowtie2 index for RepEnrich2 to build:
        1. A fasta with all repeats combined (~2.5GB), and fasta index.
        2. A gtf file with the locations annotated.

        Parameters:
            repEnrichBowtieIndexFolder: Folder containing *.fa files.
            repeats_artificial_chromosome_folder: The output folder.
        """

        if repeats_artificial_chromosome_folder is None:
            repeats_artificial_chromosome_folder = Path('repeats_artifical_chromosome/')

        output_fasta = Path(repeats_artificial_chromosome_folder, 'repeats_chrom.fa')
        output_gtf = Path(repeats_artificial_chromosome_folder, 'repeats.gtf')

        # Get filenames for all the fasta files used to build the index.
        fnames = [Path(x) for x in glob.glob(repEnrichBowtieIndexFolder + '/*fa')]
        # Get a dict of {seq_name: sequence}
        print(f"Loading {len(fnames)} fastas...")
        print(fnames)
        fastas = {}
        [fastas.update(fastaReader(fname)) for fname in fnames]

        print(fastas['AluSz6'][:500])
        print("Finished loading fastas.")
        repEnrichBowtieIndexFolder = Path(repEnrichBowtieIndexFolder)

        os.makedirs(repeats_artificial_chromosome_folder, exist_ok=True)

        seq_order = sorted(fastas.keys())

        print(f"Sorted seq_order: {seq_order}")

        Ns = 'N' * 500
        positions_in_chrom = []
        artificial_chrom = ''
        for n, name in enumerate(seq_order):
            positions_in_chrom.append(
                {'name': name, 'start': len(artificial_chrom), 
                'end': len(artificial_chrom) + len(fastas[name])})
            artificial_chrom += fastas[name] + Ns
            
        with open(Path(repeats_artificial_chromosome_folder, 'repeats_chrom.fa'), 'w') as f:
            f.write('>repeats\n')
            [f.write(artificial_chrom[i:i+250] + '\n') for i in range(0, len(artificial_chrom), 250)]

        with open(Path(repeats_artificial_chromosome_folder, 'repeats.gtf'), 'w') as f:
            for repeat in positions_in_chrom:
                chrom, name = 'repeats', repeat['name']

                # GTF is 1-based.
                start = repeat['start'] + 1
                end = repeat['end'] + 1

                li = f'{chrom}\tENSEMBL\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_support_level "NA"\n'
                li += f'{chrom}\tENSEMBL\ttranscript\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_type "repeat"; transcript_name "{name}; transcript_support_level "NA"; tag "basic";\n'
                li += f'{chrom}\tENSEMBL\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}"; gene_type "repeat"; gene_name "{name}"; transcript_type "repeat"; transcript_name "{name}; exon_number 1; exon_id "{name}.1"; transcript_support_level "NA"; tag "basic";\n'
                f.write(li)

        # Create the *.fa.fai fasta index.
        os.system(f'samtools faidx {output_fasta}')

class RepEnrichMethods():
    """Methods for mapping to repeats using RepEnrich2, and processing output.
    Typically, the order would be:
    map_to_repeats_with_bowtie2(fastq filenames)  # Writes data/pair1_/*sam
    combine_RepEnrich2_sams(data/pair1_/)  # Writes all_repeats.sam
    combine_genomic_and_repeats_bams(genomic_bam, repeats_sam)  # Writes a sams/all_reads.sam
    So this goes from fastq files to sams/all_reads.sam, and requires the genomic bam already
    exits.

    map_to_repeats_with_bowtie2(): fastqfile names -> multiple output sam files.
        Parameters:
            fastqfile1, fastqfile2
        Outputs:
            Path('data/pair1_/')

    combine_RepEnrich2_sams(): Combines the sams of mapping all reads to each repeat individually.
        Parameters:
            folder_of_sams: str = 'data/pair1_/',
            output_filename: str = 'sam_subset_to_assigned_reads/all_repeats.sam'
        Outputs:
            The output_filename sam file is created.

    combine_genomic_and_repeats_bams(): Combines the filtered genomic sam and single repeats sam.
        Parameters:
            genomic_bam: str = 'sams/genome.Aligned.out.bam',
            repeats_sam: str = 'sam_subset_to_assigned_reads/all_repeats.sam',
            merged_filename: str = 'sams/all_reads.bam'
        Outputs:
            repeats_after_removing_genomic_sam = Path(os.path.dirname(repeats_sam), 'reads_not_in_genomic_mapping.sam')
            repeats_after_removing_genomic_bam = os.path.splitext(repeats_after_removing_genomic_sam)[0] + ".bam"
            merged_filename ('sams/all_reds.bam') is created.

    count_genomic_instances_of_repeats(): Create list of repeat elements sorted by # of genomic instances.
        This is used to assign reads to repeats by priority when combining sam files.
        Input is a file of genomic instances for each repeat type. This file is generally ~500 MB.
        To save time, a file of a previous calculation result can be loaded if passed as a
        filename to 'saved_result'.
        Parameters:
            annotation_file: str = 'hg38_repeatmasker_clean.txt',
            saved_result: str=None
        Outputs:
            Returns [[repeat1_name, ...], [repeat_family_1_name, ...]]
    """

    def map_to_repeats_with_bowtie2(
        self, fastqfile1: str, fastqfile2: str, cpus: int=10):

        annotation_file = self.file_paths['repeats annotation']
        index_folder = self.file_paths['bowtie2 repeats index']
        outputfolder = Path(self.file_paths['sams'], 'RepEnrich_output/')
        genomic_mapping_bam = Path(self.file_paths['sams'], 'genome.filtered.bam')

        # Just FYI: an example command to call Bowtie to map to repeats.
        #cmd = "bowtie2 -q -p 1"  # -q = fastq input. -p is CPU number.
        #cmd += f" -x {index_folder}"
        #cmd += f" -1 {fastqfile1} -2 {fastqfile2}"
        #cmd += f" -S {self.file_paths['sams']}/repeats_bowtie.sam"

        os.makedirs(outputfolder, exist_ok=True)

        # Write an empty bam file with just the header and no reads.
        # RepEnrich2 expects to extract reads from the genomic mapping,
        # so it requires this file.
        header_from_genomic_bam = Path(self.file_paths['sams'], 'genomic_header.bam')
        cmd = f'samtools view -o {header_from_genomic_bam} -H {genomic_mapping_bam}'
        p = subprocess.check_output(cmd.split(' '))#, shell=True, stdout=stdout)

        # Use the argparse.Namespace class to call RepEnrich2.
        args = Namespace(
            annotation_file = str(annotation_file),
            outputfolder = str(outputfolder),
            outputprefix = 'sample_name_prefix',
            setup_folder = str(index_folder),
#            repeat_bed = setup_folder + os.path.sep + 'repnames.bed' 
            alignment_bam = str(header_from_genomic_bam),  # Empty.
            fastqfile = str(fastqfile1),
            fastqfile2 = str(fastqfile2),
            cpus = cpus,
            b_opt = "-k 1 -p 1 --quiet --fast",
            collapserepeat = 'Simple_repeat',
            pairedend = 'TRUE',
            allcountmethod = 'FALSE',
            is_bed = 'FALSE'
        )
        print("Mapping to repeats with bowtie:")
        pprint(args)

        RepEnrich2.run(args)

        """
        RepEnrich2 writes to:
        outputfolder + os.path.sep + outputfile_prefix 
        where we are using ./outputfolder/sample_name_prefix currently (defined in the Namespace args).
        Output files include ./outputfolder/{prefix}_fraction_counts.txt, which contains the estimated
        counts for each repeat. RepEnrich2 analysis uses those files as inputs to R analysis
        for differential enrichment.

        Also output in {outputfolder}/ are pair1_/, pair2_/ and sorted_/ folders of sam files.
        For our purposes, we want to separate reads that mapped to the genome as well,
        and create a single sam file of reads mapped to each element. Reads are of course
        counted multiple times if they map to multiple elements.

        Potential solutions include:
        1. For each read, find the set of 'best' mappings, take one of those and discard the rest.
        2. Establish a priority for each repeat class, and take the mapping of top priority.

        The second method is likely to be vastly faster than the first. A priority could be established
        based on instances in the genome.

        hg38_repeatmasker_clean.txt gives genomic instances.
        
        An example of a commmand line call to RepEnrich2.py:
        cmd = "python RepEnrich2.py "
        cmd += " hg38_repeatmasker_clean.txt /data/Sample_Output_Folder "
        cmd += " Sample_Name /data/RepEnrich2_setup_mm9 "
        cmd += " data/sample_name_multimap_R1.fastq --fastqfile2 /data/sample_name_multimap_R2.fastq "
        cmd += " /data/sample_name_unique.bam  --cpus 16 --pairedend TRUE"
        """

        # Sets self.sam_repeats_filename (output_filename) and creates that file.
        self.combine_RepEnrich2_sams(
            folder_of_sams=Path(args.outputfolder, 'pair1_'),
            output_filename=Path(self.file_paths['sams'], 'all_repeats.sam'))
        
    # Print and run a system command.
    def proclaim(self, cmd):
        print(cmd)
        os.system(cmd)

    def combine_genomic_and_repeats_bams(
        self,
        genomic_bam: Union[str, Path] = 'sams/genome.filtered.bam',
        repeats_sam: Union[str, Path] = 'sams/all_repeats.sam',
        merged_filename: Union[str, Path] = 'sams/all_reads.bam') -> None:
        
        # For a specific set of repeat elements, we use their mapping to repeats, not the genome.
        # Specifically, 'rRNA', 'scRNA', 'snRNA', 'tRNA'. special_repeats is a global variable
        # for this module.

        sf = self.file_paths['sams']

        reads_to_keep_in_repeats = set()
        for repeat in special_repeats:
            fname = Path(sf, 'assigned_read_names', repeat)
            if not os.path.exists(fname):
                continue
            reads_to_keep_in_repeats |= set([x.rstrip('\n') for x in open(fname).readlines()])

        
        bam = pysam.Samfile(genomic_bam, 'rb')
        
        bam_kept_genomic = pysam.AlignmentFile(f"{sf + '/kept_genomic.bam'}", "wb", template=bam)
        bam_removed_genomic = pysam.AlignmentFile(f"{sf + '/removed_as_repeats_genomic.bam'}", "wb", template=bam)
        
        for read in bam.fetch():
            if read.query_name in reads_to_keep_in_repeats:
                bam_removed_genomic.write(read)
            else:
                bam_kept_genomic.write(read)

        bam.close()
        bam_kept_genomic.close()
        bam_removed_genomic.close()

        # Some intermediate files.
        repeats_after_removing_genomic_sam = Path(os.path.dirname(repeats_sam), 'reads_not_in_genomic_mapping.sam')
        repeats_after_removing_genomic_bam = os.path.splitext(repeats_after_removing_genomic_sam)[0] + ".bam"

        # Compare which reads are mapped where.
        def cf_reads():

            # The -F 4 flag gets only mapped reads.
            # sort -u does the same thing as sort | uniq
            # cut -f1 gets the read name columns.
            cmd = f"samtools view -F 4 {sf + '/kept_genomic.bam'} | cut -f1 | sort -u > {sf}/genomic_reads.txt"
            self.proclaim(cmd)

            cmd = f"samtools view -F 4 {repeats_sam} | cut -f1 | sort -u > {sf}/repeats_reads.txt"
            self.proclaim(cmd)
   
            # These comm outputs are just a single column of read names.
            cmd = f"comm -12 {sf}/genomic_reads.txt {sf}/repeats_reads.txt > {sf}/AlignToBoth.txt"
            self.proclaim(cmd)
            cmd = f"comm -23 {sf}/genomic_reads.txt {sf}/repeats_reads.txt > {sf}/AlignToGenomicOnly.txt"
            self.proclaim(cmd)
            cmd = f"comm -13 {sf}/genomic_reads.txt {sf}/repeats_reads.txt > {sf}/AlignToRepeatsOnly.txt"
            self.proclaim(cmd)

        cf_reads()

        # All the reads mapped to repeats, regardless of whether they map to the genome as well.
        repeats = set([x.rstrip('\n') for x in open(f'{sf}/AlignToRepeatsOnly.txt').readlines()])

        # Open a sam file to write the filtered reads to (those not in the genomic bam).
        outf = open(repeats_after_removing_genomic_sam, 'w')
        
        # For each read in the repeats sam, write that read iff it's not in the genomic bam.
        with open(repeats_sam) as f:

            # Grab the header.
            header = []
            for li in f:
                if li[0] != '@':
                    break
                else:
                    header.append(li)

            # Write the header.
            outf.write(''.join(header))

            print(f"Wrote {len(header):,} header lines to {repeats_after_removing_genomic_sam}.")

            # Write out the first read to a sam.
            if li.split('\t')[0] in repeats:
                outf.write(li)

            # Write the rest of the reads.
            for li in f:
                if li.split('\t')[0] in repeats:
                    outf.write(li)

        outf.close()  # Finished writing the filtered reads.

        # Now we convert to having all repeats mapped to an artificial chromosome.
        #g = sameRiver.mapping.repeatsGenome()  # Init does nothing.
        #g.convert_sam_file_with_repeats_as_chroms_into_single_chrom(
        #    sam_fname='sams/all_reads.sam', gtf_fname='./repeats_artifical_chromosome/repeats.gtf',
        #    out_sam_fname='test.sam', artificial_chom_fasta_index='./repeats_artifical_chromosome/repeats_chrom.fa.fai')

        # Sort and covert to bam the filtered reads mapping to repeats.
        cmd = f'samtools sort {repeats_after_removing_genomic_sam} > {repeats_after_removing_genomic_bam}'
        self.proclaim(cmd)

        cmd = f"samtools index {repeats_after_removing_genomic_bam}"
        self.proclaim(cmd)

        # Combine the genomic reads and the reads mapping to repeats.
        # -f forces the merge even if merged_filename exists.
        cmd = f'samtools merge -f {merged_filename} {genomic_bam} {repeats_after_removing_genomic_bam}'
        self.proclaim(cmd)

    def combine_RepEnrich2_sams(
        self, folder_of_sams: str = 'data/pair1_/',
        output_filename: str = 'sam_subset_to_assigned_reads/all_repeats.sam'):
        """For a sam/*sam folder, combine into a single sam.
        """

        sf = self.file_paths['sams']  # For legibility.

        self.sam_repeats_filename = Path(output_filename)

        if not hasattr(self, 'repeats_listed_by_priority'):
            self.count_genomic_instances_of_repeats(
                annotation_file=self.file_paths['repeats annotation'])

        print(f"Combining sam files in {folder_of_sams} and outputting to {output_filename}")
        # Get a list of sam files output by Bowtie2 calls in RepEnrich2.
        sams = [Path(x) for x in glob.glob(str(folder_of_sams) + '/*.txt')]
        
        # Make output folders if needed.
        for path in [
            f'{sf}/just_read_names', f'{sf}/temp', f'{sf}/assigned_read_names', 
            os.path.dirname(output_filename)]:
            os.makedirs(path, exist_ok=True)
        
        # Make a file to output all of the mapped reads to.
        single_sam_fname = Path(os.path.dirname(output_filename), 'no_header_all_repeats.sam')
        single_sam_output = open(single_sam_fname, 'w')
        headers = open(f"{sf}/headers_temp.sam", 'w')
        
        paths_to_just_read_names = []
        all_previously_assigned_reads_fname = Path(sf, 'temp')
        already_assigned_reads = set()
        
        sum_of_new_reads, sum_reads_output_to_single_sam = (0, 0)

        # For each genome of a repetitive element, in order of decreasing priority:
        for n, repeat in enumerate(self.repeats_listed_by_priority):

            fname = Path(folder_of_sams, repeat + '.txt')
            if fname not in sams:  # Only use repeats with some matching sam file.
                continue
            
            just_names_filename = Path(sf, 'just_read_names', repeat)
            paths_to_just_read_names.append(just_names_filename)
            
            # Make a file with only the read names.
            os.system(f'cut -f1 {fname} | sort | uniq > {just_names_filename}')

            # From the file with only read names, make a set of all read names.
            with open(just_names_filename) as f:
                reads = set([x.rstrip('\n') for x in f])
            
            # Discard reads that were treated as mapped to a higher priority genome.
            new_reads = reads - already_assigned_reads
            
            # Make a file of the reads that are now assigned to this genome.
            with open(Path(sf, 'assigned_read_names', repeat), 'w') as f:
                f.write('\n'.join([x for x in new_reads]))
            
            # Add the newly assigned reads to the set of all assigned reads.
            already_assigned_reads |= new_reads

            sum_of_new_reads += len(new_reads)

            # Write the header for this same file to a file of headers, and the reads to
            # a sam file of all reads.
            with open(fname) as f:
                header = next(f)
                header = re.sub('SN:repname', f'SN:{repeat}', header)
                headers.write(header)
                for li in f:
                    s = li.split('\t')
                    if s[0] in new_reads:
                        s[2] = repeat
                        single_sam_output.write('\t'.join(s))
                        sum_reads_output_to_single_sam += 1
        
        print(f"Collected {sum_of_new_reads:,} reads and wrote {sum_reads_output_to_single_sam:,} to {output_filename}.")
        
        # Close the header file and sam reads file, and then add the header to the sam file to create
        # the final sam output.
        headers.close()
        single_sam_output.close()
        os.system(f"cat {sf}/headers_temp.sam {single_sam_fname} > {output_filename}")
    
    def count_genomic_instances_of_repeats(
        self, annotation_file: str = 'hg38_repeatmasker_clean.txt',
        saved_result: Union[str, None] = None) -> List[str]:
        """Create a list of repeat elements sorted by their number of genomic instances.
        This is used to assign reads to repeats by priority when combining sam files.
        Input is a file of genomic instances for each repeat type. This file is generally ~500 MB.
        To save time, a file of a previous calculation result can be loaded if passed as a
        filename to 'saved_result'.
        Returns [[repeat1_name, ...], [repeat_family_1_name, ...]]
        """

        # If no existing datafile, make one and then re-call this function to load it.
        if saved_result is None:

            # Set these as {repeat -> number of instances}
            repeat_instances = collections.defaultdict(int)  # By repeat, ~species.
            family_instances = collections.defaultdict(int)  # By repeat family, ~genus.

            with open(annotation_file) as fh:

                # The header is split into two lines for visualization, hence the weird column names.
                header = next(fh).split()
                header2 = next(fh)
                blank = next(fh)
                col_number_for_repeat = header.index('matching') + 1  # It's a weird format.
                col_number_for_family = header.index('repeat') + 1

                # Count instances.
                for li in fh:
                    s = li.split()
                    repeat_instances[s[col_number_for_repeat]] += 1
                    family_instances[s[col_number_for_family]] += 1

            # Export the data to a JSON file.
            os.makedirs(Path(self.file_paths['sams'], "data/"), exist_ok=True)

            datafilename = Path(
                self.file_paths['sams'], "data", "number_of_genomic_instances_for_repeats.txt")

            with open(datafilename, 'w') as f:
                json.dump([repeat_instances, family_instances], f, indent=4)

            # Also write all the repeat names to a text file. This is used to categorize
            # mapped reads in combined bed/bedgraphs as having mapped to repeats.
            with open(Path(self.file_paths['sams'], "data", "names_of_repeats.txt"), 'w') as fh:
                fh.write('\n'.join(list(repeat_instances.keys())))

            print(f"Wrote the number of genomic instances of each repeat to {datafilename}")

            return self.count_genomic_instances_of_repeats(
                annotation_file=annotation_file, saved_result=datafilename)

        # Loading an existing data file:
        print(f"Getting the number of genomic instances of each repeat from {saved_result}...")

        with open(saved_result, 'r') as f:
            [repeat_instances, family_instances] = json.load(f)

        # Make sorted lists to use to determine priority when combining sams.
        self.repeats_sorted_by_instance = sorted(
            repeat_instances, key=lambda x: repeat_instances[x])[::-1]
        self.repeat_family_sorted_by_instance = sorted(
            family_instances, key=lambda x: family_instances[x])[::-1]

        # Some repeats, (rRNA, snRNA, scRNA and tRNA) are not prioritized by the number 
        # of genomic instances. rRNA goes first, then snRNA, then the other repeats we are 
        # treating as genes (which are subsequently listed in order of genomic instances). 
        # special_repeats and snRNA are defined as global variables in this module.
        # Overall priority ends up being (bt2=bowtie2, s=star):
        # rRNA(bt2) > snRNA(bt2) > scRNA/tRNA(bt2) > unique genome mappings(s) > other repeats(bt2)

        special_repeats_sorted = [
            x for x in self.repeats_sorted_by_instance if x in special_repeats]

        special_repeats_order = ['LSU-rRNA_Hsa', 'SSU-rRNA_Hsa', '5S'] + [
            x for x in special_repeats_sorted if x in snRNA]

        special_repeats_order.extend([
            x for x in special_repeats_sorted if x not in special_repeats_order])

        # self.repeats_listed_by_priority ([LSU, SSU, ...]) is the essential product
        # of this method.
        self.repeats_listed_by_priority = special_repeats_order + [
            x for x in self.repeats_sorted_by_instance if x not in special_repeats_order]

        print(f"Loaded {len(self.repeats_listed_by_priority)} repeats from {saved_result}.")
        return self.repeats_listed_by_priority 


class mappingMethods(RepEnrichMethods):
    """General mapping methods for repeats and the genome.
    
    mapping(): Entry point method.
        Parameters:
            which_first: str='separate',
            fastqfile1: str = '', fastqfile2: str = ''
        Outputs:
            Calls the selected mapping method.

    map_all_reads_separately_to_genome_and_repeats():
        Parameters:
            fastqfile1: str, fastqfile2: str
        Outputs:
            Creates the merged bam/sam file self.file_paths['sams'] + /all_reads.sam

    map_to_genome_with_star(): fastqs to self.file_paths['sams'] + '/genome.Aligned.out.sam
        Parameters:
            read1_fname: str, read2_fname: str
        Outputs: The filtered alignment.
            self.file_paths['sams'] + '/genome.Aligned.out.sam
            Returns a dict of some filepaths.

    """

    def mapping(
        self, #fastq_directory: str = '',
        which_first: str='separate',
        fastqfile1: str = '', fastqfile2: str = '',
        cpus: int=10, clobber: Union[str, bool] = True):
        """Map the given fastq files in the manner indicated by which_first.
        """

        # Use default file paths if not given.
        if fastqfile1 == '':
            fastqfile1 = self.file_paths['cutadapt'] + '/concatenated_R1_reads_for_mapping.fastq'
        if fastqfile2 == '':
            fastqfile2 = self.file_paths['cutadapt'] + '/concatenated_R2_reads_for_mapping.fastq'

        # Initialize a log dict if needed.
        if not(hasattr(self, 'log')):
            self.log = {}

        # Define a sam folder if not already set.
        if 'sams' not in self.file_paths:
            self.file_paths['sams'] = self.file_paths['fastq'] + '/sam/'

        # Make any sam/bed directories that don't exist.
        for _dir in [self.file_paths['beds'], self.file_paths['sams'] + '/split/']:
            os.makedirs(_dir, exist_ok=True)

        # Mapping call.

        if which_first == 'repeats':
            self.map_to_repeats_then_genome(fastqfile1=fastqfile1, fastqfile2=fastqfile2)
        else:
            self.map_all_reads_separately_to_genome_and_repeats(
                fastqfile1=fastqfile1, fastqfile2=fastqfile2, cpus=cpus, clobber=clobber)


    def map_all_reads_separately_to_genome_and_repeats(
        self, fastqfile1: str, fastqfile2: str, cpus: int=10,
        clobber: Union[str, bool] = True) -> None:

        sf = self.file_paths['sams']  # To make this more concise.

        # If either clobber is True for the genome, or the output file (filtered)
        # from STAR mapping to the genome doesn't exist, map.
        if (clobber is True or clobber=='genome') or (
            not os.path.exists(Path(sf, 'genome.filtered.bam'))):

            output_fnames = self.map_to_genome_with_star(
                read1_fname=fastqfile1, read2_fname=fastqfile2)

        # If either clobber is True for bowtie2 mapping to repeats, or the output file 
        # (combined) doesn't exist, map.
        if (clobber is True or clobber=='repeats') or (
            not os.path.exists(Path(sf, 'all_repeats.sam'))):
    
            # Outputs to Path(self.file_paths['sams'], 'all_repeats.sam'):
            self.map_to_repeats_with_bowtie2(fastqfile1, fastqfile2, cpus=cpus)

        # If clobber is not boolean False, or the final output bam
        # doesn't exist, combine the outputs into all_reads.bam.
        # clobber=any_string will evaluate in the conditional below to True.
        if (clobber is not False) or (
            not os.path.exists(Path(sf, 'all_reads.bam'))):

            # Compare the outputs of the two independent mappings.
            # self.sam_repeats_filename: all repeats.
            # mapped_filtered_sam: genomic mappings (genome.filtered.bam).
            self.combine_genomic_and_repeats_bams(
                genomic_bam=Path(sf, 'genome.filtered.bam'),
                repeats_sam=Path(sf, 'all_repeats.sam'),
                merged_filename=Path(sf, 'all_reads.bam'))

            self.split_collapse_and_make_beds_and_bedgraphs_from_sam(   
                input_sam_file=sf + '/all_reads.bam')

    def map_to_genome_with_star(
        self, read1_fname: str, read2_fname: str) -> None:
        """Outputs to {self.file_paths['sams']}/genome/ and {self.file_paths['sams']}/.
        Makes a STAR call, then saves the raw output with suffix .not_mapq_filtered.
        Then uses samtools view to MAPQ>10 filter and output the filtered reads to 
        genome.Aligned.out.sam. Further filtering finally outputs to genome.filtered.sam.
        Returns a dict of file paths.
        """

        sf = self.file_paths['sams']

        if 'STAR' not in self.file_paths:
            self.file_paths['STAR'] = 'STAR'
            print("Assuming STAR is on PATH.")

        _starCaller = starCaller()

        cmd = _starCaller.star_cmd_to_genome(
            paths=self.file_paths,
            read1_fname=read1_fname, read2_fname=read2_fname)

        print(cmd, '\n')
        subprocess.check_output(cmd.split(' '))
        print("Finished mapping.")

        cmd = 'mv {gen} {gen}.not_mapq_filtered'.format(gen=sf + '/genome.Aligned.out.sam')
        subprocess.check_output(cmd.split(' '))

        # The cmd below excludes secondaries (-F 256) and unmapped reads (-F 4). 
        # It requires first mates (-f 64) and excludes low MAPQ (-q 10).
        # It prints the header (-h) and outputs (-o) to a file rather than STDOUT.
        cmd = 'samtools view -h -F 256 -f 64 -F 4 -q 10'
        cmd += f" -o {sf + '/genome.filtered.bam'}" # Output filename
        cmd += f" {sf + '/genome.Aligned.out.sam.not_mapq_filtered'}"  # Input.

        print(cmd)
        subprocess.check_output(cmd.split(' '))

        unmapped_fastq_fname1 = sf + '/genome.Unmapped.out.mate1'
        unmapped_fastq_fname2 = sf + '/genome.Unmapped.out.mate2'

        # Sort and index the bam file.
        cmd = f"samtools sort {sf + '/genome.filtered.bam'} > {sf + '/genome.filtered.bam.sorted'}"
        self.proclaim(cmd)
        cmd = f"mv {sf + '/genome.filtered.bam.sorted'} {sf + '/genome.filtered.bam'}"
        self.proclaim(cmd)
        cmd = f"samtools index {sf + '/genome.filtered.bam'}"
        self.proclaim(cmd)

        return {
            'unmapped_fastq1': unmapped_fastq_fname1,
            'unmapped_fastq2': unmapped_fastq_fname2,
            'mapped_filtered_bam': self.file_paths['sams'] + '/genome.filtered.bam'
            }

    def map_to_repeats_then_genome(self, fastqfile1: str, fastqfile2: str) -> None:
        """Maps the given filenames and outputs to {self.file_paths['sams']}/repeats.
        """
        # Requires bowtie2 to be on path and 'STAR'/'STAR index' to be in self.file_paths.
        # Requires samtools to be on path.

        #not_mapping = """
        sf = self.file_paths['sams']
        srm = f'{sf}/single_repeat_chromosome_mapping/'
        os.makedirs(srm, exist_ok=True)

        cmd = self.file_paths['STAR']
        cmd += f" --genomeDir {self.file_paths['STAR repeats index']}"
        cmd += ' --runThreadN 10'
        # This doesn't do anything here, but it has to be used when the genome is built:
        cmd += ' --genomeSAindexNbases 5'
        # Max intron size = 1. Setting to 0 causes the default behavior:
        cmd += ' --alignIntronMax 1'
        cmd += f' --readFilesIn {fastqfile1} {fastqfile2}'
        cmd += ' --outFilterMultimapNmax 50'
        cmd += ' --alignEndsType EndToEnd --outReadsUnmapped Fastx'
        # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
        # of a paired end read) reads in separate file(s).
        # output in separate fasta/fastq files, Unmapped.out.mate1/2
        cmd += f" --outFileNamePrefix {srm}/repeats."

        # Output sam file is always '{PREFIX}Aligned.out.sam'
        # Unmapped r1: {PREFIX}Unmapped.out.mate1
        # Unmapped r2: {PREFIX}Unmapped.out.mate2

        print(cmd)
        
        try:
            result = subprocess.check_output(cmd.split(' '))
        except subprocess.CalledProcessError as exc:
            result = exc.output
            print("\n\nError calling STAR:\n", result)
            sys.exit()

        self.map_to_genome_with_star(
            read1_fname=f'{srm}/repeats.Unmapped.out.mate1',
            read2_fname=f'{srm}/repeats.Unmapped.out.mate2',
        )
        #"""
        # Get the SQ and PG lines from the repeats header.
        seqlines = ''
        pg_line = ''
        with open(f'{srm}/repeats.Aligned.out.sam') as fh:
            for li in fh:
                if li[:3] == '@SQ':
                    seqlines += li
                if li[:3] == '@PG':
                    pg_line += li
                if li[0] != '@':
                    break

        # Make sam copy of genomic map (replace this with just using pysam to read the bam.)
        cmd = f'samtools view -h -o {sf}/genome.filtered.sam {sf}/genome.filtered.bam'
        print(cmd)
        os.system(cmd)

        # Get all the header lines from the genomic sam.
        hd_line = ''
        post_seqlines = ''
        with open(f'{sf}/genome.filtered.sam') as fh:
            for li in fh:
                if li[:3] == '@HD':
                    hd_line = li
                elif li[:3] == '@SQ':
                    seqlines += li
                elif li[0] == '@':
                    post_seqlines += li
                else:
                    break

        # Write a combined sam including the combined header and combined sam.
        with open(f'{sf}/all_reads.sam', 'w') as fh:
            fh.write(hd_line)
            fh.write(seqlines)
            fh.write(post_seqlines)
            fh.write(pg_line)

            with open(f'{sf}/genome.filtered.sam') as genomefh:
                for li in genomefh:
                    if li[0] == '@':
                        continue
                    fh.write(li)

            with open(f'{srm}/repeats.Aligned.out.sam') as repeatsfh:
                for li in repeatsfh:
                    if li[0] == '@':
                        continue
                    fh.write(li)

        self.split_collapse_and_make_beds_and_bedgraphs_from_sam(
            input_sam_file=f'{sf}/all_reads.sam')

    def split_collapse_and_make_beds_and_bedgraphs_from_sam(
        self, input_sam_file: str = None) -> None:

        if input_sam_file is None:
            input_sam_file = self.file_paths['sams'] + '/all.filtered.sam'

        if os.path.splitext(input_sam_file)[1] == '.bam':
            as_sam = os.path.splitext(input_sam_file)[0] + '.sam'
            cmd = f'samtools view -h -o {as_sam} {input_sam_file}'
            print(cmd)
            os.system(cmd)
            input_sam_file = as_sam

        # Split the giant sam file and make bed and bedgraph files.
        split_sam_filenames = sameRiver.sam_to_bed_and_wig.split_sam(
            input_sam_file,
            self.file_paths['sams'] + '/split/')

        # Collapse duplicates.
        sam_dups_dir = self.file_paths['sams'] + '/before_duplicates_collapsed/'
        os.makedirs(sam_dups_dir, exist_ok=True)

        for sam_fname in split_sam_filenames:
            dups_file = '{}/{}'.format(
                sam_dups_dir, re.sub('.sam$', '.with_dups.sam', os.path.basename(sam_fname)))
            os.system('mv {} {}'.format(sam_fname, dups_file))
            sameRiver.collapse_duplicates.collapse_sam(dups_file, sam_fname)

        # Make bed and bedgraph files.
        #for sam_fname in split_sam_filenames:
        for sam_fname in glob.glob(self.file_paths['sams'] + '/split/*sam'):
            sameRiver.sam_to_bed_and_wig.sam_to_bed_and_wig(
                sam_fname,
                '{}/{}.bed'.format(
                    self.file_paths['beds'],os.path.basename(sam_fname).split('.sam')[0]),
                self.file_paths['beds']  # Wigs directory.
                )

    #def filter_sam(self, sam_fname, out_fname):
    #    cmd = 'samtools view -h -F 256 -f 64 -F 4'
    #    cmd += ' -o {} {}'.format(
    #        out_fname, sam_fname)
    #    print(cmd)
    #    subprocess.check_output(cmd.split(' '))
        # The above excludes secondaries (-F 256) and unmapped reads (-F 4), 
        # and requires first mates (-f 64).
        # It prints the header (-h) and outputs (-o) to a file rather than STDOUT.

    def split_sam(self, sam_fname, out_dir):
        sameRiver.sam_to_bed_and_wig.split_sam(sam_fname, out_dir=out_dir)

    def cutadapt(self, fastq_input_dir=None, out_dir=None):
        
        if fastq_input_dir is None:
            fastq_input_dir = self.file_paths['fastq'] + '/ready_to_map/'
        
        if out_dir is None:
            out_dir = self.file_paths['fastq'] + '/ready_to_map/cutadapt/'
        
        sameRiver.clip_adapters.clip_using_cutadapt(fastq_input_dir, out_dir)

        self.file_paths['cutadapt'] = out_dir