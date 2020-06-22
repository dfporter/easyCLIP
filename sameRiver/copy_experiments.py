import os, sys, re, glob, collections
from pathlib import Path
from typing import List, Mapping, Union
from pprint import pprint

used_sequencing_runs = {

    'mg': {'Instrument': 'Illumina Hiseq 2500/Miseq', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/fastq/raw/R2.fastq.gz',
    },

    # l1, l2 incorporated into mg.
    #'l1': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R2.fastq.gz',
    #},

    #'l2': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R2.fastq.gz',
    #},

    'rb': {'Instrument': 'Illumina Hiseq 4000', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_rbfox_190418/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_rbfox_190418/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_rbfox_190418/fastq/raw/R2.fastq.gz',
    },

    'hp': {'Instrument': 'Illumina Hiseq 4000', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_pcbp1_190416/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_pcbp1_190416/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/hiseq_pcbp1_190416/fastq/raw/R2.fastq.gz'
    },

    '05': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180105_hiseq/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180105_hiseq/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180105_hiseq/fastq/raw/R2.fastq.gz'
    },

    '17': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180117_hiseq/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180117_hiseq/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/180117_hiseq/fastq/raw/R2.fastq.gz'
    },

    '24': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/170924_hiseq/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/170924_hiseq/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/170924_hiseq/fastq/raw/R2.fastq.gz'
    },

    'pc': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R2.fastq.gz',
    },

    'tc': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/181115/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/181115/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/181115/fastq/raw/R2.fastq.gz',
    },

    # m0, m1, m2 incorporated into mg.
    #'m0': {'Instrument': 'Illumina Miseq', 'read length': 76,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R2.fastq.gz',
    #},

    #'m1': {'Instrument': 'Illumina Miseq', 'read length': 76,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R2.fastq.gz',
    #},

    #'m2': {'Instrument': 'Illumina Miseq', 'read length': 76,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R2.fastq.gz',
    #},
}

def proclaim(cmd):
    print(cmd)
    os.system(cmd)

class copier():
    """
    paths: a dict formatted like:
        'pc': {'Instrument': 'Illumina Miseq', 'read length': 76,
        'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113',
        'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R1.fastq.gz',
        'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R2.fastq.gz',
        },
    """

    def __init__(self, paths: dict, new_top_dir: str):
        self.old_paths = paths
        new_paths = {}

        for run_name, run_dict in paths.items():
            basename = run_dict['location'].rstrip('/').split('/')[-1]

            new_paths[run_name] = {}
            new_paths[run_name].update(run_dict)
            new_paths[run_name]['location'] = str(Path(new_top_dir, basename + '/'))
            new_paths[run_name]['R1_fastq'] = str(Path(new_top_dir, basename, 'fastq', 'raw', 'R1.fastq.gz'))
            new_paths[run_name]['R2_fastq'] = str(Path(new_top_dir, basename, 'fastq', 'raw', 'R2.fastq.gz'))

            exp_top = new_paths[run_name]['location']
            os.makedirs(exp_top, exist_ok=True)
            os.makedirs(os.path.dirname(new_paths[run_name]['R1_fastq']), exist_ok=True)
            os.makedirs(f"{exp_top}/sams/", exist_ok=True)
            os.makedirs(f"{exp_top}/beds/", exist_ok=True)
            os.makedirs(f"{exp_top}/data/", exist_ok=True)
            os.makedirs(f"{exp_top}/tables/", exist_ok=True)
            os.makedirs(f"{exp_top}/fastq/ready_to_map/cutadapt/", exist_ok=True)

            os.system(f"touch {new_paths[run_name]['R1_fastq']}")
            os.system(f"touch {new_paths[run_name]['R2_fastq']}")

            r1 = run_dict['location'] + '/fastq/ready_to_map/cutadapt/concatenated_R1_reads_for_mapping.fastq'
            r2 = run_dict['location'] + '/fastq/ready_to_map/cutadapt/concatenated_R2_reads_for_mapping.fastq'

            #proclaim(f"rsync {r1} {exp_top}/fastq/ready_to_map/cutadapt/concatenated_R1_reads_for_mapping.fastq")
            #proclaim(f"rsync {r2} {exp_top}/fastq/ready_to_map/cutadapt/concatenated_R2_reads_for_mapping.fastq")
            proclaim(f"rsync {run_dict['location']}/scheme.xlsx {exp_top}/scheme.xlsx")

        pprint(new_paths)

if __name__ == '__main__':
    if input(f"Going to copy sequencing runs to a new top directory of {sys.argv[1]}. OK [Y/N]?:") == 'Y':
        c = copier(used_sequencing_runs, sys.argv[1])

#
