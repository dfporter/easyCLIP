#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 15:47:58 2017

@author: dfporter
"""
import re, sys, os, glob, subprocess, collections

rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'N': 'N'}
    
    
def rc(seq):
    seq = ''.join([rc_dict[base] for base in seq])
    return seq[::-1]


def create_architecture(top_level_dir):
    os.system('mkdir {0}'.format(top_level_dir))
    for subdir in ['raw', 'r2_clipped', 'r1r2_clipped']:
        os.system('mkdir {0}/{1}'.format(top_level_dir, subdir))


def read_a_read(fh):
    
    try:
        name = next(fh)
    except:
        return False, False, False
    
    seq = next(fh).rstrip('\n')
    next(fh)
    qual = next(fh).rstrip('\n')
    
    return name, seq, qual
    
    
def clip_by_r2(fname_r1, fname_r2, out_dir):
    """Clips R2 using the random adapter in R1, and filters out short inserts.
    """

    log = collections.defaultdict(int)
    log.update({x:y for x,y in locals().items() if x!='log'})

    fastq_file_r1 = open(fname_r1, 'r')
    fastq_file_r2 = open(fname_r2, 'r')
    out_fname_r1 = f"{out_dir}/{os.path.basename(fname_r1)}"
    out_fname_r2 = f"{out_dir}/{os.path.basename(fname_r2)}"



    n = -1
    keep = []


    # Write a new, clipped, R2 file, discarding reads with short inserts.
    fh = open(out_fname_r2, 'w')
    while True:
        name1, seq1, qual1 = read_a_read(fastq_file_r1)

        name2, seq2, qual2 = read_a_read(fastq_file_r2)
        if (not name1) or (not name2):
            break

        n += 1
        to_search = rc(seq1[6:13])
        seq2 = seq2.split(to_search)[0]
        
        if len(seq2) <= 6:
            log['R2 <=6 nt after removing R1 hexamer, discarded'] += 1
            continue
        log['R2 kept'] += 1

        qual2 = qual2[:len(seq2)]
        keep.append(n)
        fh.write(f'{name2}{seq2}\n+\n{qual2}\n')

    fh.close()
    fastq_file_r1.close()
    fastq_file_r2.close()

    log['Reads input to clip_by_r2'] = log['R2 <=6 nt after removing R1 hexamer, discarded'] + log['R2 kept']
    if log['Reads input to clip_by_r2']>0:
        log['% Kept after removing hexamer'] = 100 * log['R2 kept']/log['Reads input to clip_by_r2']
    else:
        log['% Kept after removing hexamer'] = 100

    # Write a new R1 file (not clipped), keeping only the inserts in the new R1 file.
    fastq_file_r1 = open(fname_r1, 'r')
    n = -1
    on_element = 0
    fh = open(out_fname_r1, 'w')
    while len(keep)>(on_element):
        n += 1
        if (n == keep[on_element]):
            name1, seq1, qual1 = read_a_read(fastq_file_r1)

            if not name1:
                break

            fh.write('{0}{1}\n+\n{2}\n'.format(name1, seq1, qual1))
            on_element += 1
            continue

        elif(n<keep[on_element]):
            name1, seq1, qual1 = read_a_read(fastq_file_r1)
            continue

        else:
            break

    fh.close()
    fastq_file_r1.close()
    fastq_file_r2.close()

    #check_file_lengths(out_fname_r1, out_fname_r2)
    return log

def check_file_lengths(out_fname_r1, out_fname_r2):
    n_lines_r2 = len(open(out_fname_r2, 'r').readlines())
    n_lines_r1 = len(open(out_fname_r1, 'r').readlines())
    if n_lines_r1 != n_lines_r2:
        print("Lines in R1 and R2 don't match:\n{0}:{1}\n{2}:{3}\n".format(
            os.path.basename(out_fname_r1), n_lines_r1, 
            os.path.basename(out_fname_r2), n_lines_r2))
        sys.exit()
    #else:
    #    print("Line numbers:\n{0}:{1}\n{2}:{3}\n".format(
    #        os.path.basename(out_fname_r1), n_lines_r1, 
    #        os.path.basename(out_fname_r2), n_lines_r2))


def clip_by_r1(fname_r1, fname_r2, out_dir, bc_l3=None):
    
    log = collections.defaultdict(int)
    log.update({x:y for x,y in locals().items() if x!='log'})


    fastq_file_r1 = open(fname_r1, 'r')
    fastq_file_r2 = open(fname_r2, 'r')
    
    out_fname_r1 = "{0}/{1}".format(out_dir, os.path.basename(fname_r1))
    out_fname_r2 = "{0}/{1}".format(out_dir, os.path.basename(fname_r2))
    
    if bc_l3 is None:
        bc_l3 = fname_r1.split('_')[-1].split('.fastq')[0]
    
    #print("clip_by_r1 running on L3 barcode {}...".format(bc_l3))
    # rc TCA = TGA. Adapter read from R1 as TGAAGATC
    try:
        adapter_end = rc(bc_l3) + 'AGATC'
    except:
        print("""{}:
Looked for a L3 barcode nt sequence in {} but couldn't understand.
Using AGATC instead.""".format(
        fname_r1, bc_l3))
        adapter_end = 'AGATC'
    n = -1

    # Write a new R1 file removing the L3 adapter, and filter by insert length.
    out_fh = open(out_fname_r1, 'w')
    keep = []
    while True:
        name1, seq1, qual1 = read_a_read(fastq_file_r1)
        if not name1:
            break
        n += 1
        seq1 = seq1.split(adapter_end)[0]

        if len(seq1) < (13+7):
            log['R1 < 20 nt, discarded'] += 1
            continue
        log['R1 kept'] += 1

        keep.append(n)
        qual1 = qual1[:len(seq1)]
        out_fh.write('{0}{1}\n+\n{2}\n'.format(name1, seq1, qual1))

    out_fh.close()

    fastq_file_r1.close()
    fastq_file_r2.close()

    log['Reads input to clip_by_r1'] = log['R1 < 20 nt, discarded'] + log['R1 kept']
    if log['Reads input to clip_by_r1']>0:
        log['% Kept'] = 100 * log['R1 kept']/log['Reads input to clip_by_r1']
    else:
        log['% Kept'] = 100

    # Write a new R2 file using only the reads kept by R1.
    out_fh = open(out_fname_r2, 'w')
    fastq_file_r2 = open(fname_r2)
    n = -1
    on_element = 0
    while len(keep) > on_element:
        n += 1
        if (n == keep[on_element]):
            name2, seq2, qual2 = read_a_read(fastq_file_r2)
            out_fh.write('{0}{1}\n+\n{2}\n'.format(name2, seq2, qual2))
            on_element += 1
            continue
        elif(n<keep[on_element]):
            name2, seq2, qual2 = read_a_read(fastq_file_r2)
            continue
        else:
            break
    out_fh.close()
    fastq_file_r2.close()

    #check_file_lengths(out_fname_r1, out_fname_r2)
    
    return log
def file_list_paired(input_dir):
    file_org = {}
    for fname in glob.glob(input_dir + '/*.fastq'):
        basename = os.path.basename(fname).split('.fastq')[0]
        if re.search('_R2', basename):
            true_base_name = re.sub('_R2', '', basename)
            print(true_base_name)
            file_org.setdefault(true_base_name, {})
            file_org[true_base_name]['R2'] = fname
        else:
            file_org.setdefault(basename, {})
            file_org[basename]['R1'] = fname
    return file_org

def clip_using_cutadapt(input_dir, out_dir, min_length=17):
    """Cutadapt discards reads<min_length."""
    
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    if not os.path.exists(out_dir + '/too_short_cutadapt_output/'):
        os.system('mkdir ' + out_dir + '/too_short_cutadapt_output/')

    if input_dir.rstrip('/') == out_dir.rstrip('/'):
        print("Input and output directories cannot be the same.")
        return

    for fastq_r2 in glob.glob(input_dir + '/*_R2_*.fastq'):
        fastq_r1 = re.sub('_R2', '', fastq_r2)
        basename1 = os.path.splitext(os.path.basename(fastq_r1))[0]
        basename2 = os.path.splitext(os.path.basename(fastq_r2))[0]

        cmd = 'cutadapt -A CTATACTACCCTTCGCTTCACACACACAAG -a CTTGTGTGTGTGAAGCGAAGGGTAGTATAG'
        cmd += ' --pair-filter=any'
        cmd += ' --too-short-output {} --too-short-paired-output {}'.format(
            out_dir + f'/too_short_cutadapt_output/{basename1}_too_short.fastq',
            out_dir + f'/too_short_cutadapt_output/{basename2}_too_short_R2.fastq')

        cmd += ' --minimum-length {} -o {} -p {}'.format(
            min_length,            
            out_dir + '/' + os.path.basename(fastq_r1),
            out_dir + '/' + os.path.basename(fastq_r2))

        cmd += f' {fastq_r1} {fastq_r2}'
        
        res = subprocess.check_output(cmd.split(' '))
        res = res.decode()

def clip_a_file_pair_using_cutadapt(fastq_r1, fastq_r2, out_dir, min_length=17):
    """Cutadapt discards reads<min_length."""

    basename1 = os.path.splitext(os.path.basename(fastq_r1))[0]
    basename2 = os.path.splitext(os.path.basename(fastq_r2))[0]
    out_fname1 = out_dir + '/' + os.path.basename(fastq_r1)
    out_fname2 = out_dir + '/' + os.path.basename(fastq_r2)
    too_short_fname1 = out_dir + f'/{basename1}_too_short_cutadapt_output.fastq'
    too_short_fname2 = out_dir + f'/{basename1}_too_short_cutadapt_output_R2.fastq'

    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)

    cmd = 'cutadapt -A CTATACTACCCTTCGCTTCACACACACAAG -a CTTGTGTGTGTGAAGCGAAGGGTAGTATAG'
    cmd += ' --pair-filter=any'
    cmd += f' --too-short-output {too_short_fname1}'
    cmd += f' --too-short-paired-output {too_short_fname2}'
    cmd += f' --minimum-length {min_length} -o {out_fname1} -p {out_fname2}'
    cmd += f' {fastq_r1} {fastq_r2}'
    
    res = subprocess.check_output(cmd.split(' '))
    res = res.decode()

    # Not used at the moment:
    try:
        total_pairs = re.search('Total read pairs processed:\s+([\d,]+)\n', b)
        total_pairs = re.sub(',', '', total_pairs.group(1))
        total_pairs = int(total_pairs)
    except:
        total_pairs = 0.

    return out_fname1, out_fname2

def clip_a_dir_by_r2(top_level_dir, input_r1r2_split_dir):
    
    for p6_bc_dir in glob.glob(input_r1r2_split_dir + '/*'):  
        print("Clipping adapters by R2 in {0}...".format(p6_bc_dir))
        
        r1_bc = p6_bc_dir.split('/')[-1]
        r2_clipped_dir = top_level_dir + '/r2_clipped/' + r1_bc
        os.system('mkdir ' + r2_clipped_dir)
        
        file_org = file_list_paired(p6_bc_dir)
        for p6p3_split_dir, _fastqs in list(file_org.items()):
            if ('R1' not in _fastqs) or ('R2' not in _fastqs):
                print('Skipping {0} as files were missing.'.format(
                    p6p3_split_dir))
                continue
            clip_by_r2(_fastqs['R1'], _fastqs['R2'], r2_clipped_dir)


def clip_a_dir_by_r1(top_level_dir, r1r2_clipped_output_dir):
    
    for p6_bc_dir in glob.glob(top_level_dir + '/r2_clipped/*'):    
        print("Clipping adapters by R1 in {0}...".format(p6_bc_dir))
        
        r1_bc = p6_bc_dir.split('/')[-1] 
        r1r2_clipped_dir = r1r2_clipped_output_dir + '/' + r1_bc
        
        os.system('mkdir ' + r1r2_clipped_dir)
        
        file_org = file_list_paired(p6_bc_dir)
        
        for p6p3_split_dir, _fastqs in list(file_org.items()):
        
            if ('R1' not in _fastqs) or ('R2' not in _fastqs):
                print('Skipping {0} as files were missing.'.format(
                    p6p3_split_dir))
                continue
            
            clip_by_r1(_fastqs['R1'], _fastqs['R2'], r1r2_clipped_dir)


def clip_r1r2(top_level_dir, r1r2_clipped_output_dir='./', input_r1r2_split_dir=None):
    if input_r1r2_split_dir is None:
        input_r1r2_split_dir = top_level_dir + '/r1r2_split/'
    create_architecture(top_level_dir)
    clip_a_dir_by_r2(top_level_dir, input_r1r2_split_dir)
    clip_a_dir_by_r1(top_level_dir, r1r2_clipped_output_dir)


if __name__ == '__main__': 
    top_level_dir = sys.argv[1]
    input_r1r2_split_dir = sys.argv[2]
    clip_r1r2(top_level_dir, input_r1r2_split_dir)
