import re, pandas, HTSeq, glob
import numpy as np

def for_split_lines(fname, skip_first_line=False):
    with open(fname) as f:
        if skip_first_line:
            next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            yield s
            
def genes_with_snoRNA_annotation(gtf_file='/opt/genome/refseq/GRCh38_latest_genomic.gff',
                                 out_file='/opt/genome/snoRNA_lines.gtf'):
    outli = ''
    with open(gtf_file) as f:
        for li in f:
            if re.search('snoRNA', li):
                outli += li
    with open(out_file, 'w') as f:
        f.write(outli)


def refseq_ids_of_snoRNAs(gtf_file='/opt/genome/snoRNA_lines.gtf'):
    refseqs = set()
    with open(gtf_file) as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            refseq = re.search('Name=([^;]+)', li)
            if refseq is not None:
                refseqs.add(refseq.group(1).split('.')[0])
    return refseqs


def genomic_coordinates_of_snoRNAs(gtf_file='/opt/genome/snoRNA_lines.gtf'):
    coords = {}
    with open(gtf_file) as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            refseq = re.search('Name=([^;]+)', li)
            if refseq is not None:
                name = refseq.group(1).split('.')[0]
                coords[name] = HTSeq.GenomicInterval(s[0], int(s[3]), int(s[4]), s[6])
    return coords


def seqs_with_snoRNA_annotation(
        fa_file='/opt/indexes/refseq_index/refseq_rna_reformatted.fasta',
        out_fa_file='/opt/genome/refseq/snoRNA.fasta'):
    ids = refseq_ids_of_snoRNAs()
    seqs = {}
    name = ''
    with open(fa_file) as f:
        for li in f:
            if li[0] == '>':
                name = li[1:].rstrip('\n')
                name = name.split('.')[0]
                if name in ids:
                    seqs[name] = ''
            else:
                if name in ids:
                    seqs[name] += li.rstrip('\n')
    with open(out_fa_file, 'w') as f:
        for name, val in seqs.items():
            f.write('>{0}\n{1}\n'.format(name, val))
            

def fasta(fa_file='/opt/genome/refseq/snoRNA.fasta'):
    seqs = {}
    name = ''
    with open(fa_file) as f:
        for li in f:
            if li[0] == '>':
                name = li[1:].rstrip('\n')
                name = name.split('.')[0]
                seqs[name] = ''
            else:
                seqs[name] += li.rstrip('\n')
    return seqs
