import collections, os, HTSeq, sys, re
            
def read_sam_generator(fname):
    """SAM is 1-based and the end is included.
    Bed is 0-based and end is excluded.
    BAM is 0-based and end is excluded."""
    with open(fname, 'r') as f:
        for li in f:
            if li[0] == '@':
                continue

            s = li.rstrip('\n').split('\t')

            strand = '+'
            if int(s[1]) & 16: # flag = s[1]
                strand = '-'            

            read_len = 0
            for match in re.finditer('(\d+)[MDN]', s[5]):
                read_len += int(match.groups()[0])

            if read_len == 0:
                yield None

            if int(s[3]) <= 0:
                print("SAM file {0} contains left index below zero: {1}".format(fname, li))
                yield None

            # num_reads = int(re.match('.*n=(\d+).*', li).groups()[0])
            # Switch to 0-based and end excluded here.
            yield (
                ''.join(re.split('(rand=[^-]+)-', s[0])[:2]),
#                s[0].split('-')[0],
                HTSeq.GenomicInterval(s[2].split('.')[0],  int(s[3]) - 1, int(s[3]) - 1 + read_len, strand))

def split_sam(fname, out_dir='sams/'):

    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)

    header = ''
    out_sam = {}
    out_sam_fnames = []

    with open(fname, 'r') as f:
        for li in f:

            if li[0] == '@':
                header += li
                continue

            s = li.split('\t')    
            barcode = s[0].split('rand')[0]

            if barcode not in out_sam:
                out_sam.setdefault(
                    barcode, open('{}/{}.sam'.format(out_dir, barcode), 'w'))
                out_sam[barcode].write(header + '@SQ     SN:repeats      LN:4301367\n')
                out_sam_fnames.append('{}/{}.sam'.format(out_dir, barcode))

            out_sam[barcode].write(li)

    for barcode in out_sam:
        out_sam[barcode].close()

    out_sam_fnames = list(set(out_sam_fnames))  # Shouldn't need this line.

    return out_sam_fnames

def sam_to_bed_and_wig(fname, out_bed_fname, wigs_dir):

    out_bed = open(out_bed_fname, 'w')

    for _dir in [wigs_dir, wigs_dir + '/read_start', wigs_dir + '/entire_read_len_counted']:
        if not os.path.exists(_dir):
            os.system('mkdir ' + _dir)

    ga_read_len = HTSeq.GenomicArray('auto', stranded=True)
    ga_read_start = HTSeq.GenomicArray('auto', stranded=True)

    for line_num, result in enumerate(read_sam_generator(fname)):

        if result is None:
            continue
            
        (name, read_iv) = result
        
        if read_iv.start == read_iv.end:
            continue

        if read_iv.start < 0:
            print('how did this happen?', read_iv)
            continue
        
        ga_read_len[read_iv] += 1
        
        out_bed.write(
                '\t'.join([read_iv.chrom, str(read_iv.start), str(read_iv.end), name.split('rand=')[1], '1',
                           read_iv.strand]) + '\n')

        #read_iv.start = read_iv.start + int(float(read_iv.end - read_iv.start)/2.)
        
        if read_iv.strand == '+':
            ga_read_start[HTSeq.GenomicPosition(read_iv.chrom, read_iv.start, read_iv.strand)] += 1
        else:  # Adjust for the fact that the .bed for a read is always [a, b).
            ga_read_start[HTSeq.GenomicPosition(read_iv.chrom, read_iv.end-1, read_iv.strand)] += 1
            
    out_bed.close()

    basename = os.path.basename(out_bed_fname).split('.bed')[0]
    ga_read_len.write_bedgraph_file(
        '{}/entire_read_len_counted/{}_+.wig'.format(wigs_dir, basename), strand='+')
    ga_read_len.write_bedgraph_file(
        '{}/entire_read_len_counted/{}_-.wig'.format(wigs_dir, basename), strand='-')

    ga_read_start.write_bedgraph_file(
        '{}/read_start/{}_+.wig'.format(wigs_dir, basename), strand='+')
    ga_read_start.write_bedgraph_file(
        '{}/read_start/{}_-.wig'.format(wigs_dir, basename), strand='-')


def read_sam(fname):
    counts = collections.defaultdict(int)
    ga_read_len = {}
    ga_read_start = {}
    
    out_bed = {}

    for _dir in ['beds/', 'beds/read_start', 'beds/entire_read_len_counted']:
        if not os.path.exists(_dir):
            os.system('mkdir ' + _dir)

    for line_num, result in enumerate(read_sam_generator(fname)):

        if result is None:
            continue
            
        (name, read_iv) = result
        
        if read_iv.start == read_iv.end:
            continue

        if read_iv.start < 0:
            print('how did this happen?', read_iv)
            continue

        barcode = name.split('rand')[0]
        ga_read_len.setdefault(barcode, HTSeq.GenomicArray('auto', stranded=True))
        ga_read_start.setdefault(barcode, HTSeq.GenomicArray('auto', stranded=True))
        
        ga_read_len[barcode][read_iv] += 1
        
        if barcode not in out_bed:
            out_bed[barcode] = open('beds/' + barcode + '.bed', 'w')
            out_bed[barcode].write(
                '\t'.join([read_iv.chrom, str(read_iv.start), str(read_iv.end), name.split('rand=')[1], '1',
                           read_iv.strand]) + '\n')
        else:
            out_bed[barcode].write('\t'.join([read_iv.chrom, str(read_iv.start), str(read_iv.end), name.split('rand=')[1], '1',
                           read_iv.strand]) + '\n')

        #read_iv.start = read_iv.start + int(float(read_iv.end - read_iv.start)/2.)
        read_iv.end = read_iv.start + 1
        if read_iv.strand == '+':
            ga_read_start[barcode][HTSeq.GenomicPosition(read_iv.chrom, read_iv.start, read_iv.strand)] += 1
        else:
            ga_read_start[barcode][HTSeq.GenomicPosition(read_iv.chrom, read_iv.end, read_iv.strand)] += 1
            
    for barcode, fh in out_bed.items():
        fh.close()

    for barcode, ga in ga_read_len.items():
        ga.write_bedgraph_file('beds/entire_read_len_counted/' + barcode + '_+.wig', strand='+')
        ga.write_bedgraph_file('beds/entire_read_len_counted/' + barcode + '_-.wig', strand='-')
        
    for barcode, ga in ga_read_start.items():
        ga.write_bedgraph_file('beds/read_start/' + barcode + '_+.wig', strand='+')
        ga.write_bedgraph_file('beds/read_start/' + barcode + '_-.wig', strand='-')

if __name__ == '__main__':
    read_sam(sys.argv[1])
