# Readme

Code for the development of the easyCLIP method and the analysis of easyCLIP datasets.

A reproducible workflow intended to be easily usable is being constructed in snakemake in a separate repository:
https://github.com/dfporter/easyclip-snakemake

## Organization

[Genome build](doc/README_genome.md)

[Raw reads to bed/bedgraph/ann_counts files](doc/README_raw_fastq_to_beds_and_counts.md)

[ann_counts files to pvalues](doc/README_ann_counts_to_stats.md)

[Motif identification](doc/README_motif_identification.md)

[Biotypes](doc/README_biotypes.md)

[Adapter fluorescence calculations](doc/README_fluorescence.md)

[Datasets and GEO upload](doc/README_datasets.md)

Code is divided into jupyter notebooks stored in dataAndScripts/,
 scripts for passing commands on a server,
 and everything else, stored in sameRiver/.

Currently processing before analysis is built around submitting jobs on a server because of the
 high RAM required for STAR.
Analysis is built around Jupyter notebooks expected to be run locally on 
 (1) annotated counts files, (2) bedgraph files, and
 (3) a scheme file -
 all but the last of these are built by server commands. 

The current pipeline is:
(1) Raw fastq files, not split by barcode, and (2) a scheme file made by hand
 -> a main_x.py script that makes an exp object and processes 
 -> output the counts/bed/bedgraph file for analysis.

A scheme.xlsx file is included.

```bash
# Dependencies.
# Python3.6+
pip install numpy scipy cutadapt pandas statsmodels xlrd dill
pip install HTSeq openpyxl

git clone https://github.com/dfporter/easyCLIP
git clone https://github.com/dfporter/RepEnrich2.git --branch py3
```

Here's an example of a main_x.py file:

```python
import os, sys, re, glob, pandas, importlib
import sameRiver
import sameRiver.exp

#################################################################
# Read processing and mapping. 
#################################################################

importlib.reload(sameRiver.exp)

top_dir = '/oak/stanford/groups/khavari/users/dfporter/seq/v2all/'

ex = sameRiver.exp.exp(
    name='v2', file_paths=sameRiver.exp.exp.autogenerate_paths(top_dir))

# Longer form, with explicit path names:
"""
ex = sameRiver.exp.exp(name='hp', file_paths={
    'scheme': top_dir + '/scheme.xlsx',
    'beds': top_dir + '/beds/',
    'fastq': top_dir + '/fastq/',
    'R1_fastq': top_dir + '/fastq/raw/R1.fastq',
    'R2_fastq': top_dir + '/fastq/raw/R2.fastq',
    'sams': top_dir + '/sams/',
    'counts': top_dir + '/counts.txt',
    'STAR index': '/oak/stanford/groups/khavari/users/dfporter/genome/GRCh38.gencode.29/star_index',
    'STAR': 'STAR',
    'STAR repeats index': '/oak/stanford/groups/khavari/users/dfporter/genome/repeats/star_repeats',
    })
"""

ex.read_scheme()
ex.split_by_barcode()
ex.convert_split_r1r2_folder_to_long_filenames()
ex.preprocess_split_barcodes_folder()
ex.mapping() #top_dir + '/fastq/ready_to_map/cutadapt/')
ex.line_numbers()  # Optional: write the size of bed files.

#################################################################
# Analysis 
#################################################################

# Get a combined gtf. This function also sets up a folder for IGV to make a .genome file.
# The important output is repeats_and_genome.gtf, which holds the RNA information
# for both the genome and the repeats-as-chromosomes RepEnrich reads.
g = sameRiver.mapping.repeatsGenome()
g.setup_genomes(
        repeats_fasta_directory='hg38re/',
        genomic_gtf='gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA',
        igv_output_directory='for_igv/',
        repeats_gtf='repeats_as_separate_chroms.gtf',
        combined_gtf='repeats_and_genome.gtf',
        )

# Make data/rna.data file to assign reads to genes.
import sameRiver.rnaDataFileMaker

maker = sameRiver.rnaDataFileMaker.rnaDataFileMaker()
RNAs = maker.make_from_gtf_file(
    gtf_filename='repeats_and_genome.gtf',
    #gtf_filename='/oak/stanford/groups/khavari/users/dfporter/genome/repeats_and_ensembl_release94_GRCh38/combined.gtf',
#    gtf_filename='/opt/genome/repeats_and_ensembl_release94_GRCh38/combined.gtf'
)

# -> outputs data/rna.data. 
# This holds gtf object information - just regions, really.

#####
# Make counts file.
# Outputs a data/bed_x.data file that holds signal, w/o RNA information:
ex.make_signal_data_file()
# Assigns signal to genes and writes counts.txt:
ex.make_scheme_signal_RNA_data_files(
    rna_data_object=RNAs)

# Make ann_counts.txt file. This has simplified
# column names and biotypes.
ex.annotate_counts_file()


```

Here's a typical workflow for analyzing a set of positives against
 some negatives, both having come from experiments that genereated
 annotated counts files with exp.annotate_counts_file().

```python
import importlib
import sameRiver
import sameRiver.metadata
import sameRiver.metadata.negative_metadata as negative_metadata
import sameRiver.metadata.positive_metadata_pcbp1_190416 as positive_metadata
import sameRiver.negativeCounts
import sameRiver.positiveCounts
import sameRiver.scheme
import sameRiver.statsForCountsNB
importlib.reload(sameRiver.statsForCountsNB)
importlib.reload(sameRiver.scheme)

# If never run before:
negatives = sameRiver.negativeCounts.negativeCounts(
    negative_metadata, xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx')

# Optional: write_txt=True to write some txt's of the data.
negatives.save(write_object=True, write_txt=True)

# If loading:
negatives = sameRiver.negativeCounts.negativeCounts.load()

# If never run before:
positives = sameRiver.positiveCounts.positiveCounts(
    positive_metadata, xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx')
positives.save(write_object=True, write_txt=True)

# If loading:
positives = sameRiver.positiveCounts.positiveCounts.load()

#nb = sameRiver.statsForCountsNB.statsForCountsNB.load()
nb = sameRiver.statsForCountsNB.statsForCountsNB(
    negatives=negatives, positives=positives, data_dir=positive_metadata.top_dir + '/data/')
nb.calculate_pvalues(which='per_read')
nb.write_targets(which='per_read', outfname='default')
nb.calculate_pvalues(which='per_protein')
nb.write_targets(which='per_protein', outfname='default')
nb.save()


```

Here's a positive_metadata.py file:
```python


top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/Runs/hiseq_pcbp1_190416/'
scheme_file = top_dir + '/scheme.xlsx'
ann_counts_file = top_dir + '/ann_counts.txt'
bed_file_dir = top_dir + '/beds/'

positive_proteins = [
    'PCBP1', 'PCBP1:100P', 'PCBP1:100Q',
    'PCBP1:dKH', 'CELF1',
]
```

Here's a negative_metadata.py file:
```python
top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/all/'

scheme_file_with_random_proteins = top_dir + '/scheme.xlsx'
ann_counts_file = top_dir + '/ann_counts.txt'
bed_file_dir = top_dir + '/beds/'

random_proteins = [
    'CAPNS2', 'CCIN', 'CDK4', 'CHMP3',
    'DCTN6', 'EPB41L5', 'ETS2', 'IDE',
    'ITPA', 'TPGS2', 'UBA2',
    ]
```

After main_x.py and analysis_x.py are run, jupyter notebooks are use for everything 
 else, including all of the figure and table creations.



