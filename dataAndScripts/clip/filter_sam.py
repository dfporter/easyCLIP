import sys, subprocess

# Get only primary alignments:
# samtools view -F 256 > out.sam
# Get only the first mate in a pair:
# samtools view -f 64 > out.sam
# Get only the second mate in a pair:
# samtools view -F 64 > out.sam
# -F exludes the bit, -f requires the bit.
# First in pair: 64
# Second in pair: 128
# Can require the primary alignment and first in pair:
# samtools -F 256 -f 64 > out.sam
# The above excludes secondaries (-F 256) and requires first mates (-f 64).

def by_samtools(sam_fname, out_fname):
    cmd = 'samtools view -F 256 -f 64'
    cmd += ' -o {} {}'.format(
        out_fname, sam_fname)
    subprocess.check_output(cmd.split(' '))

def filter_out_secondary_alignments_and_mates(sam_fname):
    outf = open('filtered.sam', 'w')
    with open(sam_fname) as fh:
        for n, li in enumerate(fh):
            if not (n % 1E6):
                print(n)
            if li[0] == '@':
                outf.write(li)
                continue
            s = li.split('\t')
            # s[1] is FLAG
            try:
                if int(s[1]) & 256:  # Secondary alignment.
                    continue
            except:
                print(li)
            if int(s[1]) & 128:  # Second in pair.
                continue
            outf.write(li)
    outf.close()

if __name__ == '__main__':
    filter_out_secondary_alignments_and_mates(sys.argv[1])

