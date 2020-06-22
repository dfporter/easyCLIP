import re, sys, HTSeq, collections, pandas, os, glob

rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'N': 'N'}


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
    
    return df


def count_barcodes(read1, read2=None):
    
    fastq_file = HTSeq.FastqReader(read1)
    
    barcode_counts = collections.defaultdict(int) 
    barcode_reg_counts = collections.defaultdict(int) 
    
    for n, read in enumerate(fastq_file):
        
        if (n > 1000000):
            break 
            
        barcode_reg = read.seq[:10] 
        barcode = read.seq[:3]
        barcode_counts[barcode] += 1 
        barcode_reg_counts[barcode_reg] += 1 
        
    with open('barcode_counts.txt', 'w') as f: 
        for barcode in sorted(barcode_counts, key=lambda x: barcode_counts[x]):
        
            if barcode_counts[barcode] < 100:
                continue
                
            f.write("{b}\t{c}\n".format(b=barcode, c=str(barcode_counts[barcode]))) 
            
    with open('barcode_region_counts.txt', 'w') as f: 
        for barcode in barcode_reg_counts: 
            f.write("{b}\t{c}\n".format(b=barcode, c=str(barcode_reg_counts[barcode]))) 

def split_a_dir_by_R2(input_dir, output_dir, scheme_fname):
    
    file_org = {}
    
    print(f"Spliting by R2:\nInput:{input_dir}\nOutput:{output_dir}\n")
    
    if not os.path.exists(output_dir):
        os.system('mkdir ' + output_dir)
        
    #scheme = read_scheme(scheme_fname, rev_comp=True)
    
    # Organize filenames in the directory.
    for fname in glob.glob(input_dir + '/*.fastq'):
        print(f"L5 split file: {fname}.")
        
        basename = os.path.basename(fname).split('.fastq')[0]
        no_suffix = fname.split('.fastq')[0]
        
        if basename[-3:] == '_R2':
            true_base_name = basename[:-3]
            file_org.setdefault(true_base_name, {})
            file_org[true_base_name]['R2'] = fname
            
        else:
            file_org.setdefault(basename, {})
            file_org[basename]['R1'] = fname
            
    print(file_org)
    
    # Split files.
    for p6_bc in file_org:
        
        output_sub_dir = f'{output_dir}/{p6_bc}'

        if not os.path.exists(output_sub_dir):
            os.system('mkdir ' + output_sub_dir)
        
        read1, read2 = file_org[p6_bc]['R1'], file_org[p6_bc]['R2']
        split_a_pair_of_files_by_R2(
            read1, read2, output_dir=output_sub_dir,
            fname_pattern_for_no_recognized_barcode=f'{p6_bc}_')


def split_a_pair_of_files_by_R2(
    fname_R1, fname_R2, output_dir=None, allow_mismatch=False,
    fname_pattern_for_no_recognized_barcode=None):
    
    if output_dir is None:
        output_dir = os.path.dirname(fname_R1)

    to_num = {
        'TCA': 9,
        'AGT': 10,
        'AAC': 11, 
        'CAG': 12}

    outfiles_R2 = {}
    
    if fname_pattern_for_no_recognized_barcode is not None:
        int_to_fileh_R1 = {-1: open(
            output_dir + f'/{fname_pattern_for_no_recognized_barcode}no_recognized_barcode.fastq', 'w')}
        missingf_R2 = open(output_dir + f'/{fname_pattern_for_no_recognized_barcode}no_recognized_barcode_R2.fastq', 'w')
    
    else:
        int_to_fileh_R1 = {-1: open(
            output_dir + '/no_recognized_barcode.fastq', 'w')}
        missingf_R2 = open(output_dir + '/no_recognized_barcode_R2.fastq', 'w')
    
    # Used to increase efficiency when building the list of what read
    # goes where in the R1 file based on R2.
    bc_to_int = {'None': -1}

    # Open output files.
    for barcode in set(['CAG', 'AGT', 'AAC', 'TCA']):
        
        outfiles_R2[barcode] = open(
             output_dir + '/' + os.path.basename(fname_R2).split('.fastq')[0] + \
             '_' + str(barcode) + '.fastq', 'w')
        
        new_int = max(bc_to_int.values()) + 1
        
        bc_to_int[barcode] = new_int
        
        int_to_fileh_R1[new_int] = open(
             output_dir + '/' + os.path.basename(fname_R1).split('.fastq')[0] + \
             '_' + str(barcode) + '.fastq', 'w')

    # Include all single base variants of the 3 bp barcodes.
    if allow_mismatch:
        bases = ['A', 'T', 'C', 'G', 'N']
        barcodes = set(outfiles_R2.keys())

        for barcode in barcodes:
            bc_list = list(barcode)

            for pos in range(len(barcode)):
                for base in bases:

                    new_bc = bc_list[:]
                    new_bc[pos] = base
                    new_bc = ''.join(new_bc)

                    if new_bc not in outfiles_R2:
                        outfiles_R2[new_bc] = outfiles_R2[barcode]
                        bc_to_int[new_bc] = bc_to_int[barcode]


    fastq_file_R2 = HTSeq.FastqReader(fname_R2)
    sort_read1 = []
    
    # Go through the input R2 file.
    for n_read, read in enumerate(fastq_file_R2):
        
        if(not (n_read % 100000)):
            print(("R2: Read: {0} million".format(float(n_read)/(10**6))))
            
        bc = read.seq.decode('utf-8')[:3]
        
        if bc in outfiles_R2:
            read.write_to_fastq_file(outfiles_R2[bc])
            sort_read1.append(bc_to_int[bc])
        
        else:
            read.write_to_fastq_file(missingf_R2)
            sort_read1.append(-1)
    
    # Close output R2 files.
    missingf_R2.close()
    for bc in outfiles_R2:
        outfiles_R2[bc].close()
    
    # Go through input R1 file.
    if fname_R1 is not None:
        fastq_file_R1 = HTSeq.FastqReader(fname_R1)
        
        for n_read, read in enumerate(fastq_file_R1):
            if(not (n_read % 100000)):
                print(("R1: Read: {0} million".format(float(n_read)/(10**6))))
                
            read.write_to_fastq_file(int_to_fileh_R1[sort_read1[n_read]])
            
        for _int in int_to_fileh_R1:
            int_to_fileh_R1[_int].close()


def count_diversity(bcs): 
    pass 
 

if __name__ == '__main__': 
    read1 = sys.argv[1]
    read2 = sys.argv[2]
    results_fname = sys.argv[3]
    
    for f in [read1, read2]:
        if not os.path.exists(f):
            
            print(("{0} does not exist.".format(f)))
            sys.exit()
            
    split_a_dir_by_R2(read1, read2, results_fname)

