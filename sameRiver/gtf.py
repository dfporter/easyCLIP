import re, sys

def subset_to_only_tsl1_and_NA(fname: str):
    outf = open(fname + '.exons_only_tsl1andNA', 'w')

    with open(fname) as f:
        for li in f:
            
            if (('transcript_support_level "1"' in li) or ('transcript_support_level "NA"' in li)):
                
                if li.split('\t')[2] == 'exon':# or s[2] == 'intron'
                    outf.write(li)

    outf.close()

if __name__ == '__main__':
    subset_to_only_tsl1_and_NA(sys.argv[1])