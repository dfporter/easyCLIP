used_sequencing_runs = {

    'l1': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R2.fastq.gz',
    },

    'l2': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R2.fastq.gz',
    },

    'm0': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R2.fastq.gz',
    },

    'm1': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R2.fastq.gz',
    },

    'm2': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R2.fastq.gz',
    },
}

import os, sys, re

import exp
import metaExp
from pprint import pprint as pp

def read_a_read(fh):
    
    try:
        name = next(fh).rstrip('\n')
    except:
        return False, False, False
    
    seq = next(fh).rstrip('\n')
    next(fh)  # The ever-useful '+' line in fastq files.
    qual = next(fh).rstrip('\n')

    return name, seq, qual

def concat(fastq_file_r1, fastq_file_r2, out_read1_fname, out_read2_fname):
    
    fastq_fh_r2 = open(fastq_file_r2, 'r')
    fastq_fh_r1 = open(fastq_file_r1, 'r')
        
    out_read1 = open(out_read1_fname, 'a')
    out_read2 = open(out_read2_fname, 'a')
    
    while True:
        name1, seq1, qual1 = read_a_read(fastq_fh_r1)
        name2, seq2, qual2 = read_a_read(fastq_fh_r2)
        
        if not name1:  # End of file.
            break
        
        #prefix = re.search('(@[\w-]+rand=[ATCGN]+)-', name1.split(' ')[0]).groups()[0]
        #name1 = prefix + ' ' + name1.split(' ')[1]
        #name2 = prefix + ' ' + name2.split(' ')[1]
        
        out_read1.write(f'{name1}\n{seq1}\n+\n{qual1}\n')
        out_read2.write(f'{name2}\n{seq2}\n+\n{qual2}\n')

def combine_medgenome_sequencing_and_miseq(meta, out_top_dir='combined/'):
    """lane1, lane2, miseq0, miseq1, miseq2"""
    
    to_combine = {}  # {protein: {basename -> exp_name}}
    filenames_to_concat = {}  # {basename -> (R1 file path, R2 file path)}
    
    for exp_name in ['l1', 'l2', 'm0', 'm1', 'm2']:
        df = meta.exps[exp_name].scheme.scheme_df
        
        for protein, long_fname_R1, long_fname_R2 in zip(
            df.Gene, df['long_fname_R1'], df['long_fname_R2']):
            
            to_combine.setdefault(protein, {})
            to_combine[protein].setdefault(long_fname_R1, [])
            filenames_to_concat.setdefault(long_fname_R1, [])
            
            to_combine[protein][long_fname_R1].append(exp_name)
            filenames_to_concat[long_fname_R1].append(
                (meta.exps[exp_name].file_paths['cutadapt'] + '/split/' + long_fname_R1,
                meta.exps[exp_name].file_paths['cutadapt'] + '/split/' + long_fname_R2)
                )
            
    pp(to_combine)

    os.makedirs(f"{out_top_dir}/fastq/raw", exist_ok=True)
    for fname in ['R1.fastq', 'R2.fastq']:
        if not os.path.exists(f"{out_top_dir}/fastq/raw/{fname}"):
            os.system(f"touch {out_top_dir}/fastq/raw/{fname}")
            
    out_dir = f'{out_top_dir}/fastq/ready_to_map/cutadapt/'
    os.makedirs(out_dir, exist_ok=True)
    
    out_read1_fname = f"{out_dir}/concatenated_R1_reads_for_mapping.fastq"
    out_read2_fname = f"{out_dir}/concatenated_R2_reads_for_mapping.fastq"

    os.system(f'touch {out_read1_fname}')
    os.system(f'touch {out_read2_fname}')
    
    print(filenames_to_concat)

    for list_of_fname_tuples in filenames_to_concat.values():
        for (r1_fname, r2_fname) in list_of_fname_tuples:
            print(r1_fname, r2_fname)
            concat(r1_fname, r2_fname, out_read1_fname, out_read2_fname)

def run():
    seq_dir = "/oak/stanford/groups/khavari/users/dfporter/seq/"
    if not os.path.exists(seq_dir):
        raise IOError(f"input seq directory does not exist: {seq_dir}")
        
    meta = metaExp.metaExp(
        file_paths={'top_dir': f'{seq_dir}/hiseq/200522_medgenome/combined_with_miseq/'})

    for name, paths in used_sequencing_runs.items():
        meta.add_exp(
            exp.exp(
                name=name,
                file_paths=exp.exp.autogenerate_paths(paths['location'])
            )
       )
        meta.exps[name].read_scheme()

    combine_medgenome_sequencing_and_miseq(
        meta, out_top_dir=f'{seq_dir}/hiseq/200522_medgenome/combined_with_miseq')

    meta.combine_scheme_files()

if __name__ == '__main__':
    run()
