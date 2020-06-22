"""Go from bed+bedgraph+scheme files to counts.txt and ann_counts.txt.

Usage:
    make_counts_file.py [--genome_gtf=<gtf>] [--mapping=<txt>] [--clobber] [--input=<dir>] [--out=<dir>] [--name=<str>]
    make_counts_file.py rna [--genome_gtf=FILE] 
    make_counts_file.py meta [--genome_gtf=<gtf>] [--mapping=<txt>] [--clobber] [--input=<dir>] [--out=<dir>] [--name=<str>]


Options:
    -h --help            Show this message and exit.
    --genome_gtf=FILE   Genomic [default: combined.gtf]
    --mapping=FILE   Repeats [default:  enst_transcript_id_name_biotype_map.txt]
    --input=STR   Input folder [default: .]
    --clobber   Clobber. [default: False]
    --out=STR   Output folder. [default: ./meta/]
    --name=STR   Exp name.  [default: meta]
"""


import os, sys, re
import importlib
import sameRiver
import sameRiver.rnaDataFileMaker
import sameRiver.set_of_named_mRNAs
import sameRiver.exp
import sameRiver.metaExp
from docopt import docopt

importlib.reload(sameRiver.exp)
importlib.reload(sameRiver.rnaDataFileMaker)
importlib.reload(sameRiver.set_of_named_mRNAs)

def xls_or_xlsx_scheme(_dir):
    scheme_path = _dir + 'scheme.xls'
    if not os.path.exists(_dir + '/scheme.xls'):
        if os.path.exists(_dir + '/scheme.xlsx'):
            scheme_path = _dir + '/scheme.xlsx'
    return scheme_path

args = docopt(__doc__)
print(args)

if (not os.path.exists('./data/rna.data')) or (args['--clobber']):

    maker = sameRiver.rnaDataFileMaker.rnaDataFileMaker()
    RNAs = maker.make_from_gtf_file(
    #gtf_filename='/labs/khavari/dfporter/genome/repeats_and_ensembl_release94_GRCh38/headed.gtf'
        gtf_filename=args['--genome_gtf']
    )

if args['rna']:
    sys.exit()


top_dir = args['--input']

ex = sameRiver.exp.exp(name=args['--name'], file_paths={
    'scheme': xls_or_xlsx_scheme(args['--input']),
    #'scheme': '/Users/dfporter/pma/dataAndScripts/clip/miseq/meta/scheme.xls',
    'beds': top_dir + '/beds/',
    'wigs': top_dir + '/beds/read_start/',
    #'wigs': '/Users/dfporter/pma/dataAndScripts/clip/miseq/Runs/180105/beds/read_start_filtered/',
    'counts': top_dir + '/counts.txt',
    'xl_rate_fname': top_dir + '/percentCrosslinked.xlsx',
    #'r1r2_clipped_fastq': '/Volumes/Seagate/hiseq/180105/fastq/r1r2_clipped/'
    'fastq': top_dir + '/fastq/',

    'sams': top_dir + '/sams/',
    #'r1r2_split': '/Volumes/Seagate/hiseq/to_upload/heads/r1r2_split/',

    'STAR index': '/labs/khavari/dfporter/genome/GRCh38.gencode.29/star_index',
    #'STAR index': '/opt/genome/star_mito_index/',
    
    'STAR': 'STAR',
    
    'STAR repeats index': '/labs/khavari/dfporter/genome/repeats/star_repeats/',
    #'STAR repeats index': '/opt/genome/star_repeats/',
    })

_top = args['--out']

if not os.path.exists(_top):
    os.system('mkdir ' + _top)



meta = sameRiver.metaExp.metaExp(
    name='meta', file_paths={
        'scheme': xls_or_xlsx_scheme(_top),
        'beds': _top + '/beds/',
        'wigs': _top + '/beds/read_start/',
        'counts': _top + '/counts.txt'
    })

meta.add_exp(ex)
meta.process(clobber=args['--clobber'], verbose=True)
