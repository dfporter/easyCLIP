import HTSeq, sys, os
from pathlib import Path, PurePath


def chromosome_names_and_lengths_from_fasta(fasta_fname):
    sequences = dict( (s[1], s[0]) for s in HTSeq.FastaReader(fasta_fname, raw_iterator=True) )
    sequences = {chrom_name:len(seq) for chrom_name, seq in sequences.items()}
    with open(
        PurePath(os.path.dirname(fasta_fname), os.path.basename(fasta_fname) + '.chrom_lengths'),
        'w') as f:
        f.write(
        	'\n'.join( [f'{k}\t{v}' for k,v in sequences.items()] )
        	)

def chomosomes(fname='/opt/genomes/gencode.v29/GRCh38.primary_assembly.genome.fa.gz.chrom_lengths'):
	chrom = {}
	with open(fname) as f:
		for li in f:
			s = li.rstrip('\n').split('\t')
			chrom[str(s[0])] = int(s[1])

	return chrom

if __name__ == '__main__':
	chromosome_names_and_lengths_from_fasta(sys.argv[1])