import re, sys, os, glob, subprocess, collections

def read_a_read(fh):

	try:
		name = next(fh)
	except:
	    return False, False, False

	seq = next(fh).rstrip('\n')
	next(fh)
	qual = next(fh).rstrip('\n')

	return name, seq, qual


def remove_short_reads_in_paired_fastq(
	fname_r1, fname_r2, out_fname_r1, out_fname_r2, min_len=16):

	log = collections.defaultdict(int)
	log.update(locals())

	print(f"Removing short reads from {fname_r1} and {fname_r2}...")

	fastq_file_r1 = open(fname_r1, 'r')
	fastq_file_r2 = open(fname_r2, 'r')
	fh1 = open(out_fname_r1, 'w')
	fh2 = open(out_fname_r2, 'w')

	while True:
		name1, seq1, qual1 = read_a_read(fastq_file_r1)

		name2, seq2, qual2 = read_a_read(fastq_file_r2)

		if (not name1) or (not name2):	
			break
		
		if (len(seq1)>min_len) and (len(seq2)>min_len):
			fh1.write('{0}{1}\n+\n{2}\n'.format(name1, seq1, qual1))
			fh2.write('{0}{1}\n+\n{2}\n'.format(name2, seq2, qual2))
			log['Reads > min_len'] += 1
		else:
			log['Reads <= min_len'] += 1

	log['Sum'] = log['Reads > min_len'] + log['Reads <= min_len']

	fh1.close()
	fh2.close()
	fastq_file_r1.close()
	fastq_file_r2.close()

	print(f"Finished removing short reads and output into {out_fname_r1} and {out_fname_r2}.")

	return log