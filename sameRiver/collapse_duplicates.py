import sys, collections, os, glob, re


def collapse_bed(in_filename, out_filename):
    print('Collapsing duplicates in {0} and writing to {1}...'.format(
        in_filename, out_filename))

    outfile = open(out_filename, 'w')
    reads = collections.defaultdict(dict)
    
    with open(in_filename, 'r') as fh:
        for m, li in enumerate(fh):

            s = li.split('\t')    # Barcode is 11 nt
            #  0: chrom, 1: pos, 5: strand, 3: barcode
            reads[tuple([s[0], s[1], s[5]])][s[3]] = li
            
    for bc_dicts in list(reads.values()):
        for li in list(bc_dicts.values()):
            outfile.write(li) 
            
            
def collapse_bed_dir(top_level_dir, out_dir):
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)

    for in_filename in glob.glob(top_level_dir + '/*.bed'):
        out_filename = out_dir + '/' + os.path.basename(in_filename)
        collapse_bed(in_filename, out_filename)
        
    #for subdir in glob.glob(top_level_dir + '/*/'):
    #    out_subdir = out_dir + '/' + subdir.split('/')[-2] + '/'
    #    
    #    os.system('mkdir ' + out_subdir)
    #    for in_filename in glob.glob(subdir + '/*.bed'):
    #        out_filename = out_subdir + os.path.basename(in_filename)
    #        collapse_bed(in_filename, out_filename)
            
            
def collapse_sam(in_filename, out_filename):
    print('Collapsing duplicates in {0} and writing to {1}...'.format(
        in_filename, out_filename))

#    infile = open(in_filename, 'r')
    outfile = open(out_filename, 'w')
    headerli = ''
    reads = collections.defaultdict(dict)

    with open(in_filename, 'r') as fh:
        for m, li in enumerate(fh):

            if li[0] == '@':
                headerli += li
                continue

            bc = re.search('rand=([ATCGN]+)', li)
            s = li.split('\t')    # Barcode is 11 nt

            # _chr = s[2]; _pos = s[3]
            reads[tuple([s[2], s[3]])][bc.group(1)] = li

    outfile.write(headerli)
    for bc_dicts in list(reads.values()):
        for li in list(bc_dicts.values()):
            outfile.write(li)

            
def collapse_sam_dir(top_level_dir, out_dir):

    os.system('mkdir ' + out_dir)

    for subdir in glob.glob(top_level_dir + '/*/'):

        out_subdir = out_dir + '/' + subdir.split('/')[-2] + '/'
        print(out_subdir)
        os.system('mkdir ' + out_subdir)

        for in_filename in glob.glob(subdir + '/*.sam'):
            out_filename = out_subdir + os.path.basename(in_filename)
            collapse_sam(in_filename, out_filename)

