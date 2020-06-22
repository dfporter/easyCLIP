import re, sys, HTSeq, collections, pandas, os, glob
from typing import List, Tuple, Union

rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def rc(seq):
    seq = ''.join([rc_dict[base] for base in seq])
    return seq[::-1]


def read_scheme(filename, rev_comp=False):
    """label\tbarcode\n"""
    df = pandas.read_excel(filename, index_column=False)
    if 'Name' in df.columns:
        return df
    tups = df.to_dict('records')
    names = []
    for t in tups:
        name = "{0}-{1}-{2}".format(t['Experiment'], t['Gene'], t['Replicate'])
        while name in names:
            name += '.1'
        names.append(name)
    df['Name'] = names
    for col in ['Gene', 'std?', 'P6_BC']:
        if col in df.columns:
            df[col] = [x.decode('utf-8') for x  in df[col]]
    return df

    
def read_a_read(fh):
    try:
        name = next(fh).rstrip('\n')
    except:
        return False, False, False
    seq = next(fh).rstrip('\n')
    next(fh)
    qual = next(fh).rstrip('\n')
    return name, seq, qual

    
def move_r2_bc_to_name_in_r2(fastq_file_r1, fastq_file_r2, out_fname_r2):
    outstr = ''
    basename = os.path.basename(fastq_file_r2).split('.fastq')[0]
    fastq_fh_r2 = open(fastq_file_r2, 'r')
    fastq_fh_r1 = open(fastq_file_r1, 'r')
    basename = re.sub('_R2', '', basename)
    out_fh = open(out_fname_r2, 'w')
    while True:
        name, seq, qual = read_a_read(fastq_fh_r2)
        name1, seq1, qual1 = read_a_read(fastq_fh_r1)
        if not name:
            break
        bc = re.search('@([\w-]+rand=[ATCGN]+)-', name1).groups()[0]
        name =  '@' + bc +'-' + name.lstrip('@')
        seq = seq[3:]
        qual = qual[3:]
        out_fh.write('{0}\n{1}\n+\n{2}\n'.format(name, seq, qual))
    out_fh.close()
    fastq_fh_r1.close()
    fastq_fh_r2.close()
    
    
def move_r1_bc_to_name_in_r1(fastq_file_r1, out_fname_r1):
    outstr = ''
    
    # When outputing the L5_L3 barcode to the read name, use the file's
    # basename...
    basename = os.path.basename(fastq_file_r1).split('.fastq')[0]
    # Except if there's no recognized barcode, use "NNN" instead of
    # explicitly saying 'no_recognized_barcode', as this matches regex patterns later.
    basename = re.sub('no_recognized_barcode', 'NNN', basename)
    
    fastq_fh = open(fastq_file_r1, 'r')
    out_fh = open(out_fname_r1, 'w')

    while True:
        name, seq, qual = read_a_read(fastq_fh)
        if not name:
            break
        name = '@' + basename + 'rand=' + seq[6:13] + '-' + name.lstrip('@')
        seq = seq[13:]
        qual = qual[13:]
        out_fh.write('{0}\n{1}\n+\n{2}\n'.format(name, seq, qual))

    out_fh.close()
    fastq_fh.close()


def guess_file_pair(fname):
    basename = os.path.basename(fname)
    
    if re.search('_R2', fname):  # This is R2, guess R1.
        return os.path.dirname(fname) + '/' + re.sub('_R2', '', os.path.basename(fname))
    
    else:
        s = basename.split('_')
        if len(s) < 2:
            return None
        
        guess = os.path.dirname(fname) + '/' + '_'.join(s[:-1] + ['R2'] + [s[-1]])
        
        if not os.path.exists(guess):
            print("Couldn't find {}, trying a second filename guess...".format(guess))
            s = basename.split('.')
            guess = os.path.dirname(fname) + '/' + '.'.join(s[:-1]) + '_R2.' + s[-1]
            return guess
        else:
            return guess


def move_barcodes_to_name_in_a_directory(input_dir, output_dir):
    output_filenames = {}

    # Process R1 files.
    for fname in glob.glob(input_dir + '*.fastq'):

        # R1 only.
        if re.search('_R2', os.path.basename(fname)):
            continue

        # Find the R1 mate.
        outname = output_dir + '/' + os.path.basename(fname)
        pair = guess_file_pair(fname)

        if (not pair) or (not os.path.exists(pair)):
            print(f"Could not find the expected pair file {pair} for file {fname}, skipping")
            continue

        print(f"Moving barcode to name in R1 file {fname}...")
        move_r1_bc_to_name_in_r1(fname, outname)

        # Set output files.
        print("Output R1 file with barcode in name to {}.".format(outname))
        true_basename = os.path.basename(fname).split('.fastq')[0]
        output_filenames.setdefault(true_basename, {'R1': outname})
        output_filenames[true_basename]['R1'] = outname
        output_filenames[true_basename]['R2'] = os.path.dirname(outname) + '/' + os.path.basename(pair)

    # Process R2 files.
    for fname in glob.glob(input_dir + '*.fastq'):

        # R2 only
        if not re.search('_R2', os.path.basename(fname)):
            continue

        # Find the R1 mate.
        outname = output_dir + '/' + os.path.basename(fname)
        pair = guess_file_pair(fname)
        
        if (not pair) or (not os.path.exists(pair)):
            print("Could not find the expected pair file {0} for file {1}, skipping".format(pair, fname))
            continue

        print("R2: {0}".format(fname))
        r1_fname = re.sub('_R2', '', outname)
        print("\t...Assuming R1 is {0}".format(r1_fname))
        move_r2_bc_to_name_in_r2(r1_fname, fname, outname)

        # Set output files.
        true_basename = re.sub('_R2', '', os.path.basename(fname))
        true_basename = true_basename.split('.fastq')[0]
        output_filenames.setdefault(true_basename, {'R2': outname})
        output_filenames[true_basename]['R2'] = outname
        output_filenames[true_basename]['R1'] = os.path.dirname(outname) + '/' + os.path.basename(pair)

    return output_filenames


def make_mappable(top_dir, top_out_dir, scheme_fname):
    """Move barcodes to read names."""

    scheme = pandas.read_excel(scheme_fname)

    # Separate standards.
    scheme['std?'] = [((type(s)==type('')) and (s[0:3]=='STD')) \
                        for s in scheme.Gene]
    stds = scheme[scheme['std?']].copy()
    std_barcodes = stds['P6_BC'].tolist()

    os.system('mkdir ' + top_out_dir)

    output_filenames = dict()

    # For each file in the top directory.
    _output_filenames = make_a_directory_mappable(top_dir, top_out_dir, std_barcodes)

    if len(_output_filenames):
        output_filenames.update(_output_filenames)
    
    # For each L5-separated directory that might be under the top directory.
    for subdir in glob.glob(top_dir + '/*/'):
        print(subdir)
        print('----')

        _output_filenames = make_a_directory_mappable(subdir, top_out_dir, std_barcodes)

        if len(_output_filenames):
            output_filenames.update(_output_filenames)

    r1_fnames = [x['R1'] for x in list(output_filenames.values())]
    r2_fnames = [x['R2'] for x in list(output_filenames.values())]

    return r1_fnames, r2_fnames

def concatenate_fastqs(
    r1_fnames: List, r2_fnames: List, concat_file_r1: str, concat_file_r2: str):

    # Remove concatenated files and then insert an empty file.
    for concat_file in [concat_file_r1, concat_file_r2]:
        os.system('rm ' + concat_file)
        os.system('touch  ' + concat_file)

    # Concatenate the output files.
    cat_r1 = 'cat ' + ' '.join(r1_fnames) + ' >> ' + concat_file_r1
    print(cat_r1)
    os.system(cat_r1)

    cat_r2 = 'cat ' + ' '.join(r2_fnames) + ' >> ' + concat_file_r2
    print(cat_r2)
    os.system(cat_r2)

"""
STAR --genomeDir /seq/STAR-2.5.2b/index/hg38 --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn $file --alignEndsType EndToEnd --outReadsUnmapped Fastx
"""

def make_a_directory_mappable(subdir, top_out_dir, std_barcodes):

        if len(glob.glob(subdir + '/*.fastq')) < 1:
            print("...skipping (no fastqs).")
            return {}

        if re.search('no_recognized', subdir):
            print("skipping (no recognized barcode)")
            return {}

        # Set output directory from barcode.
        p6_bc = subdir.split('/')[-2]
        output_dir = top_out_dir + p6_bc + '/'

        if p6_bc in std_barcodes:
            print('...standard (ignoring).')
            return {}
        
        os.system('mkdir ' + output_dir)
        
        output_filenames = move_barcodes_to_name_in_a_directory(subdir, output_dir)
        
        return output_filenames

if __name__ == '__main__':
    top_dir = sys.argv[1]
    top_out_dir = sys.argv[2]
    scheme_fname = sys.argv[3]
    make_mappable(top_dir, top_out_dir, scheme_fname)

    
