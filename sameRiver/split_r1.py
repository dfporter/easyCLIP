import re, sys, HTSeq, collections, pandas, os, glob, pprint

def split_a_dir_by_R1(read1, read2=None, output_dir='./', scheme_fname=''):
    for fname in glob.glob(input_dir + '/*_R2.fastq'):
        R1 = os.path.dirname(fname) + \
            re.sub('_R2', '', os.path.basename(fname).split('.fastq')[-1])
        split_a_fastq_by_R1(R1, fname, output_dir=output_dir, scheme_fname=scheme_fname)

def fjoin(_list):
    return '/'.join([str(x) for x in _list])

def split_a_fastq_by_R1(read1, read2=None, output_dir='./', scheme_fname=''):

    print(f"Spliting by R1:\nRead1,2:{read1}, {read2}\nOutput:{output_dir}\n")

    scheme = read_scheme(scheme_fname, rev_comp=False)

    if not os.path.exists(output_dir):
        os.system('mkdir ' + output_dir)

    missingf = open(fjoin([output_dir, 'no_recognized_barcode.fastq']), 'w')
    if read2 is not None:
        missingf_r2 = open(fjoin([output_dir, 'no_recognized_barcode_R2.fastq']), 'w')
    
    (outfiles, outfiles_r2) = {}, {}

    # Open output folders.
    for barcode, name in zip(scheme.P6_BC, scheme.Name):

        if (type(name) != type('')) or (type(barcode) != type('')):
            print('Skipping row in scheme with wrong data type for name or barcode:')
            print(f'barcode {barcode} name {name}')
            continue

        # We often add a two letter prefix to the P6 barcode, using non-DNA letters.
        # We chop that off here.
        barcode = barcode.upper()
        barcode = ''.join([letter for letter in barcode if (
            letter in ['A', 'T', 'C', 'G'])])

        outfiles[barcode] = open(fjoin([output_dir, str(barcode)]) + '.fastq', 'w')

        if read2 is not None:
            outfiles_r2[barcode] = open(fjoin([output_dir, str(barcode)]) + '_R2.fastq', 'w')

    # Make copies of the L5 barcode with all possible 1 nt SNP that all
    # output to the same file. If two barcodes differ by only 1 nt then this fails,
    # but such a design should never be used anyway.
    bases = ['A', 'T', 'C', 'G', 'N']
    barcodes = set(outfiles.keys())
    for barcode in barcodes:
        bc_list = list(barcode)

        for pos in range(len(barcode)):
            for base in bases:
                new_bc = bc_list[:]
                new_bc[pos] = base
                new_bc = ''.join(new_bc)
                if new_bc not in outfiles:
                    outfiles[new_bc] = outfiles[barcode]
                    outfiles_r2[new_bc] = outfiles_r2[barcode]

    barcodes = set(outfiles.keys())
    
    n_read = 0

    # This works on .fastq.gz
    read1_iter = HTSeq.FastqReader(read1).__iter__()
    read2_iter = HTSeq.FastqReader(read2).__iter__()

    done = False
    while not(done):

        try:
            read1 = next(read1_iter)
            read2 = next(read2_iter)

        except:
            done = True
            continue

        n_read += 1

        if(not (n_read % 100000)):
            print(("R1 split: Read: {0:,} million".format(float(n_read)/(10**6))))

        bc = read1.seq[:6].decode('utf-8')

        if bc in barcodes:
            read1.write_to_fastq_file(outfiles[bc])
            read2.write_to_fastq_file(outfiles_r2[bc])
        else:
            read1.write_to_fastq_file(missingf)
            read2.write_to_fastq_file(missingf_r2)

    for bc in outfiles:
        outfiles[bc].close()
        outfiles_r2[bc].close()


rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'N': 'N'}
    
    
def rc(seq):
    seq = ''.join([rc_dict[base] for base in seq])
    return seq[::-1]


def read_scheme(filename, rev_comp=False):
    """label\tbarcode\n"""
    df = pandas.read_excel(filename, index_column=False,
        encoding='sys.getfilesystemencoding()')
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
    return df

if __name__ == '__main__':
    read1 = sys.argv[1]
    read2 = sys.argv[2]
    split_dir = sys.argv[3]
    scheme_fname = sys.argv[4]
    split_a_dir_by_R1(
        read1, read2=read2, output_dir=split_dir, scheme_fname=scheme_fname)
