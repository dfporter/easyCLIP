# easyCLIP: raw fastq to mapped bed/bedgraph files and annotated counts per RNA

The current pipeline is:
(1) Raw fastq files, not split by barcode, and (2) a scheme file made by hand
 -> a main_x.py script that makes an exp object and processes 
 -> output the counts/bed/bedgraph file for analysis.

Here's an example of a main_x.py file:

```python
import os, sys, re, glob, pandas, importlib
import sameRiver
import sameRiver.exp

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
ex.mapping(top_dir + '/fastq/ready_to_map/cutadapt/')
ex.line_numbers() # Optional: write the size of bed files.

#####
# Make data/rna.data file to assign reads to genes.
import sameRiver.rnaDataFileMaker

maker = sameRiver.rnaDataFileMaker.rnaDataFileMaker()
RNAs = maker.make_from_gtf_file(
    gtf_filename='/oak/stanford/groups/khavari/users/dfporter/genome/repeats_and_ensembl_release94_GRCh38/combined.gtf',
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
