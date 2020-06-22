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

    # These were sequenced on a hiseq in 'hp'. Should combine these with hp, but for now
    # we are simply excluding them.
    #'pc': {'Instrument': 'Illumina Miseq', 'read length': 76,
    #'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113',
    #'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R1.fastq.gz',
    #'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/190113/fastq/raw/R2.fastq.gz',
    #},

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
