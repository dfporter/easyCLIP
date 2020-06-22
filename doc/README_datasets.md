# easyCLIP: datasets

## Assembly for GEO upload

The metadata sheet for GEO upload is doc/seq_template_v2.1.xls.

To assemble the files, we first output the metadata spreadsheet,
 which is done by parsing v2all/scheme.xlsx and three files (for the
 main table, raw files table, and insert size table). These files 
 are placed by hand into the excel template provided by GEO (in doc/).

The output of this parsing is then used to select the fastq files for
 ftp upload to GEO.

The GEO upload folder was oak/seq/geo_upload/.

## Used sequencing runs

```python
used_sequencing_runs = {

    'mg': {'Instrument': 'Illumina Hiseq 2500/Miseq', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/combined_with_miseq/fastq/raw/R2.fastq.gz',
    },

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

```
## Combining 2020 sequencing

There were two Hiseq 2500 Rapid Run flow cell lanes and three miseq runs that contained some
 samples that were present in all. Each of these runs had their usual read processing done
 separately (l1, l2, m0, m1, m2). These were combined by concatenating the cutadapt output fastq
 files together (from their respective ready_to_map/cutadapt/split/ files) into concatenated
 R1 and R2 fastq files in the "mg" run (placed in its ready_to_map/cutadapt/ folder).
 These combined fastqs were 70,987,971 reads each. The script that does the combining is run
 as:

```bash
# Script in sameRiver/.
python combine_2020_sequencing.py
```

Output from mapping the combined fastqs to the genome with STAR:
```txt
                          Number of input reads |   70987971
                      Average input read length |   90

==> mg/combined_with_miseq/sams/genome.Log.progress.out <==
           Time    Speed        Read     Read   Mapped   Mapped   Mapped   Mapped Unmapped Unmapped Unmapped Unmapped
                    M/hr      number   length   unique   length   MMrate    multi   multi+       MM    short    other
May 24 14:40:27    461.7    70531728       90    66.7%     82.5     0.7%    11.4%     2.0%     0.0%    19.0%     0.9%
ALL DONE!
```

## Insert sizes

Insert sizes and std deviations were determined by running insert_size.py:
```python
import glob, os, re, sys, collections

def insert_sizes_in_file(f):
    sizes = collections.defaultdict(int)
    not_done = True
    while not_done:
         try:
             next(f)
             li = next(f).rstrip('/n')
             sizes[len(li)] += 1
             next(f)
             next(f)
         except:
             not_done = False
    return sizes

def insert_sizes_in_dir(dirname):
    flist = glob.glob(dirname + '/*fastq')
    flist = [x for x in flist if '_R2_' not in x]

    print(flist)

    sizes = collections.defaultdict(dict)
    for fname in flist:
        print(fname)

        name = os.path.basename(fname)

        if 'rbfox' in fname:
            s = name.split('_')
            s[2] = re.sub('73', 'rb', s[2])
            name = os.path.dirname(fname) + '/' + '_'.join(s)

        with open(fname) as f:
            sizes[name] = insert_sizes_in_file(f)

    return sizes

# From /labs/khavari/dfporter/hits/seqs:
all_sizes = {}
for dirname in [
    './all/fastq/r1r2_clipped/', './hiseq_rbfox_190418/fastq/r1r2_clipped/',
    './hiseq_pcbp1_190416/fastq/r1r2_clipped']:
    all_sizes.update(insert_sizes_in_dir(dirname))

outli = ""
for k, v in all_sizes.items():
    outli += "{}\t{}\n".format(k, v)
with open('insert_sizes.txt', 'w') as f:
    f.write(outli)
```

The output insert_sizes.txt file was download from the server and parsed like so (this code is from prepare_fastq_for_upload.ipynb):

```python
sizes = {}
with open('tables/insert_sizes.txt') as f:
    for li in f:
        #print(li)
        s = li.rstrip('\n').split('\t')
        s[1] = '{' + s[1].split(' {')[1].rstrip(')')#re.sub(r"defaultdict(<class 'int'>, ", "", s[1])
        s[1] = dict(eval(s[1]))
        sizes[s[0]] = s[1]
#

rows = []

for k, _dict in sizes.items():
    # _dict = {length: number_of_instances}
    lengths = np.fromiter(_dict.keys(), dtype=float)
    num_obs = np.fromiter(_dict.values(), dtype=float)
    #plt.plot(lengths, num_obs, 'k.')
    #plt.title(k)
    #plt.show()
    #plt.clf()
    #print(lengths)
    total_obs = np.sum(num_obs)
    mu = np.sum([n_obs*size for size, n_obs in zip(lengths, num_obs)])/total_obs
    var = np.average((lengths - mu)**2, weights=num_obs)
    #print(k, mu, np.sqrt(var))
    k = os.path.basename(k)
    spk = k.split('_')
    rows.append({
        'file name 1': k + '.gz',
        'file name 2': '_'.join([spk[0], spk[1], spk[2], 'R2', spk[3]]) + '.gz',
        'average insert size': np.around(mu, decimals=3),
        'standard deviation': np.around(np.sqrt(var), decimals=3),
    })
df = pandas.DataFrame(rows)

df = df.loc[:,['file name 1', 'file name 2', 'average insert size', 'standard deviation']]

df.to_excel('tables/insert_sizes_for_geo_upload.xlsx')

if os.path.exists('tables/insert_sizes_for_geo_upload.xlsx'):
    inserts = pandas.read_excel('tables/insert_sizes_for_geo_upload.xlsx')
    inserts = inserts.loc[[x in geo['raw file'].tolist() for x in inserts['file name 1']], :]
    inserts.to_excel('tables/insert_sizes_for_geo_upload_subset_to_those_included_in_upload.xlsx')
```

## Used datasets

Used datasets get written to logs/. 
 From those files we get:

```python
datasets_used_for_heatmap = ['Exp73_Rbfox1_rbCTGATC_TCA', 'Exp73_hnRNPD_rbCGATTA_AAC', 'Exp73_Rbfox1_rbCTGATC_AAC', 'Exp28_SF3B1_17AGCTAG_AGT', 'Exp73_Rbfox2_rbGTCGTC_TCA', 'Exp15_FBL_24AGCTAG_TCA', 'Exp73_Rbfox1_rbCTGATC_AGT', 'Exp71_CELF1_hpGCGGAC_TCA', 'Exp73_hnRNPD_rbCGATTA_AGT', 'Exp61_PCBP1_hpCGATTA_TCA', 'Exp73_hnRNPD_rbCGATTA_TCA', 'Exp73_Rbfox2_rbGTCGTC_AGT', 'Exp16_hnRNPC_24TGAGTG_AAC', 'Exp15_hnRNPC_24TGAGTG_TCA', 'Exp28_SF3B1_17AGCTAG_TCA', 'Exp61_PCBP1_hpCGATTA_AGT', 'Exp71_CELF1_hpGCGGAC_AGT', 'Exp15_FBL_24AGCTAG_AGT', 'Exp73_Rbfox2_rbAGCTAG_TCA']

# These FBL replicates were discarded:
# Exp16_FBL_24AGCTAG_AAC
# Exp16_FBL_24AGCTAG_CAG

# For biotypes:
datasets_for_rbps = ['Exp73_Rbfox1_rbCTGATC_TCA', 'Exp73_hnRNPD_rbCGATTA_AAC', 'Exp73_Rbfox1_rbCTGATC_AAC', 'Exp28_SF3B1_17AGCTAG_AGT', 'Exp73_Rbfox2_rbGTCGTC_TCA', 'Exp15_FBL_24AGCTAG_TCA', 'Exp73_Rbfox1_rbCTGATC_AGT', 'Exp71_CELF1_hpGCGGAC_TCA', 'Exp73_hnRNPD_rbCGATTA_AGT', 'Exp61_PCBP1_hpCGATTA_TCA', 'Exp73_hnRNPD_rbCGATTA_TCA', 'Exp73_Rbfox2_rbGTCGTC_AGT', 'Exp16_hnRNPC_24TGAGTG_AAC', 'Exp15_hnRNPC_24TGAGTG_TCA', 'Exp28_SF3B1_17AGCTAG_TCA', 'Exp61_PCBP1_hpCGATTA_AGT', 'Exp71_CELF1_hpGCGGAC_AGT', 'Exp15_FBL_24AGCTAG_AGT', 'Exp73_Rbfox2_rbAGCTAG_TCA']
datasets_for_randos = ['Exp56_EPB41L5_tcAGCTAG_TCA', 'Exp32_ITPA_05GCCTAA_TCA', 'Exp33_CHMP3_17TGTTGG_TCA', 'Exp56_CCIN_tcCTGATC_AGT', 'Exp33_ETS2_17TGCCGA_TCA', 'Exp33_IDE_17GCTCAT_TCA', 'Exp32_ITPA_05GCCTAA_AGT', 'Exp56_CCIN_tcCTGATC_TCA', 'Exp31_TPGS2_05CCACTC_CAG', 'Exp33_UBA2_17TGAGTG_AGT', 'Exp31_UBA2_05TGAGTG_AAC', 'Exp33_CHMP3_17TGTTGG_AGT', 'Exp32_CDK4_05GCCATG_TCA', 'Exp33_ETS2_17TGCCGA_AGT', 'Exp56_EPB41L5_tcAGCTAG_AGT', 'Exp31_TPGS2_05CCACTC_AAC', 'Exp32_UBA2_05TGAGTG_TCA', 'Exp33_DCTN6_17GCTGTA_TCA', 'Exp31_CDK4_05GCCATG_AAC']
datasets_for_pcbp1 = [
    'Exp71_PCBP1-dKH_hpCTGATC_AGT', 'Exp71_PCBP1-100Q_hpGTCGTC_TCA', 
    'Exp71_PCBP1-dKH_hpCTGATC_TCA', 'Exp71_PCBP1-100P_hpCCACTC_TCA', 
    'Exp71_PCBP1-100Q_hpGTCGTC_AGT','Exp71_PCBP1-100P_hpCCACTC_AGT', 
    'Exp61_PCBP1_hpCGATTA_TCA',  'Exp61_PCBP1-100P_hpGCCATG_AGT', 
    'Exp61_PCBP1_hpCGATTA_AGT', 'Exp61_PCBP1-100Q_hpAGCTAG_AGT', 
    'Exp61_PCBP1-100Q_hpAGCTAG_TCA', 'Exp61_PCBP1-100P_hpGCCATG_TCA', 
    'Exp61_PCBP1-dKH_hpATCGTG_TCA', 'Exp61_PCBP1-dKH_hpATCGTG_AGT',
    ]

# For stats:
datasets_in_negativeCounts = ['gene_name', 'Gene type', 'Exp32_ITPA_05GCCTAA_TCA', 'Exp33_CHMP3_17TGTTGG_TCA', 'Exp33_ETS2_17TGCCGA_TCA', 'Exp33_IDE_17GCTCAT_TCA', 'Exp32_ITPA_05GCCTAA_AGT', 'Exp31_TPGS2_05CCACTC_CAG', 'Exp33_UBA2_17TGAGTG_AGT', 'Exp31_UBA2_05TGAGTG_AAC', 'Exp33_CHMP3_17TGTTGG_AGT', 'Exp32_CDK4_05GCCATG_TCA', 'Exp33_ETS2_17TGCCGA_AGT', 'Exp31_TPGS2_05CCACTC_AAC', 'Exp32_UBA2_05TGAGTG_TCA', 'Exp33_DCTN6_17GCTGTA_TCA', 'Exp31_CDK4_05GCCATG_AAC']
datasets_in_positiveCounts = ['gene_name', 'Gene type', 'Exp56_EPB41L5_tcAGCTAG_TCA', 'Exp32_ITPA_05GCCTAA_TCA', 'Exp61_HCT116_hpGTAGCC_TCA', 'Exp73_Rbfox1_rbCTGATC_TCA', 'Exp33_CHMP3_17TGTTGG_TCA', 'Exp56_CCIN_tcCTGATC_AGT', 'Exp33_ETS2_17TGCCGA_TCA', 'Exp71_PCBP1-dKH_hpCTGATC_AGT', 'Exp33_IDE_17GCTCAT_TCA', 'Exp73_hnRNPD_rbCGATTA_AAC', 'Exp32_ITPA_05GCCTAA_AGT', 'Exp71_PCBP1-100Q_hpGTCGTC_TCA', 'Exp61_HCT116_hpGTAGCC_AGT', 'Exp73_Rbfox1_rbCTGATC_AAC', 'Exp56_CCIN_tcCTGATC_TCA', 'Exp71_PCBP1-dKH_hpCTGATC_TCA', 'Exp28_SF3B1_17AGCTAG_AGT', 'Exp71_PCBP1-100P_hpCCACTC_TCA', 'Exp73_Rbfox2_rbGTCGTC_TCA', 'Exp31_TPGS2_05CCACTC_CAG', 'Exp61_PCBP1-100Q_hpAGCTAG_TCA', 'Exp33_UBA2_17TGAGTG_AGT', 'Exp71_PCBP1-100Q_hpGTCGTC_AGT', 'Exp31_UBA2_05TGAGTG_AAC', 'Exp15_FBL_24AGCTAG_TCA', 'Exp73_Rbfox1_rbCTGATC_AGT', 'Exp71_CELF1_hpGCGGAC_TCA', 'Exp73_hnRNPD_rbCGATTA_AGT', 'Exp61_PCBP1_hpCGATTA_TCA', 'Exp33_CHMP3_17TGTTGG_AGT', 'Exp73_hnRNPD_rbCGATTA_TCA', 'Exp61_PCBP1-dKH_hpATCGTG_AGT', 'Exp32_CDK4_05GCCATG_TCA', 'Exp73_Rbfox2_rbGTCGTC_AGT', 'Exp16_hnRNPC_24TGAGTG_AAC', 'Exp61_PCBP1-100P_hpGCCATG_AGT', 'Exp33_ETS2_17TGCCGA_AGT', 'Exp15_hnRNPC_24TGAGTG_TCA', 'Exp28_SF3B1_17AGCTAG_TCA', 'Exp56_EPB41L5_tcAGCTAG_AGT', 'Exp61_PCBP1_hpCGATTA_AGT', 'Exp61_PCBP1-100Q_hpAGCTAG_AGT', 'Exp71_CELF1_hpGCGGAC_AGT', 'Exp71_PCBP1-100P_hpCCACTC_AGT', 'Exp61_PCBP1-100P_hpGCCATG_TCA', 'Exp31_TPGS2_05CCACTC_AAC', 'Exp32_UBA2_05TGAGTG_TCA', 'Exp33_DCTN6_17GCTGTA_TCA', 'Exp61_PCBP1-dKH_hpATCGTG_TCA', 'Exp31_CDK4_05GCCATG_AAC', 'Exp15_FBL_24AGCTAG_AGT', 'Exp73_Rbfox2_rbAGCTAG_TCA']
```
## seq_template_generated_sample_list.xlsx

```python
import collections

pcbp1s = [
    'PCBP1',
    'PCBP1 100P', 'PCBP1 100Q', 'PCBP1 dKH',
    'PCBP1:100P', 'PCBP1:100Q','PCBP1:dKH',
    'PCBP1-100P', 'PCBP1-100Q','PCBP1-dKH',
]
randos = [
    'UBA2', 'ETS2', 'CAPNS6', 'TPGS2', 
    'EPB41L5', 'CHMP3', 'CDK4', 'IDE', 
    'DCTN6', 'CCIN', 'ITPA', 
]
rbps_endog = ['FBL', 'hnRNPC', 'Rbfox2']
rbps_uORF = ['CELF1', 'Rbfox1',  'hnRNPD',]
#rbps_pLEX = ['SF3B1']

def to_expression_type(_type):
    if _type in pcbp1s:
        return 'FLAG-HA-HIS tagged ORF CRISPR targetted to AAVS1 safe harbor locus with Puro selection'
    if _type in randos:
        return 'FLAG-HA-HIS tagged ORF transiently expressed from pLEX vector 2 days without selection'
    if _type in rbps_endog:
        return 'Immunopurified endogenous protein (unmodified cells)'
    if _type in rbps_uORF:
        return 'FLAH-HA-HIS tagged ORF expressed from pLEX with a uORF in front of CDS to lower expression levels'
    if _type == 'HCT116':
        return 'Unmodified cells, immunopurified with anti-HA but without epitope'
    if re.sub('[:-]', ' ', _type) in ['SF3B1', 'SF3B1 K700E']:
        return '3X FLAG tagged SF3B1 expressed from pLEX vector for 2 days without selection'


os.chdir('/Users/dfporter/pma/dataAndScripts/clip/miseq/v2all/')

geo_upload_columns = [
    'Sample name', 'title', 'source name', 'organism',
    'characteristics: Expression',
    'characteristics: Exp',
    'characteristics: L3-BC', 
    'characteristics: L5-BC',
    'molecule', 'description', 'processed data file', 'raw file']


scheme_fname = './scheme.xlsx'
geo_fname = '/Users/dfporter/pma/doc/seq_template_v2.1.xls'

scheme = pandas.read_excel(scheme_fname)
scheme['Gene'] = [re.sub(':', '-', x) for x in scheme.Gene]

geo = scheme.copy()
geo['Sample name'] = range(1, len(geo.index)+1)

geo['source name'] = [
    {True: 'HCT116 cells', False: 'HEK293T cells'}[
        bool('PCBP1' in protein or 'HCT116' in protein)] for protein in geo['Gene']]
geo['organism'] = 'H. sapiens'
geo['characteristics: Expression'] = [to_expression_type(x) for x in geo['Gene']]
geo['characteristics: Exp'] = geo['Experiment']
geo['characteristics: Protein'] = geo['Gene']
geo['characteristics: L3-BC'] = geo['P3_BC']
geo['characteristics: L5-BC'] = [x[2:] for x in geo['P6_BC']]
geo['molecule'] = 'Library from RNA UV cross-linked to the indicated protein'
geo['processed data file'] = 'ann_counts.txt'
geo['basename'] = ['_'.join([str(_) for _ in x]) for x in zip(
    geo['Experiment'], geo['Gene'], geo['P6_BC'], geo['P3_BC'])]
geo['basename_r2'] = ['{}_{}_{}_R2_{}'.format(*x) for x in zip(
    geo['Experiment'], geo['Gene'], geo['P6_BC'], geo['P3_BC'])]

def purged(name):
    black_list = [
        'Exp16_FBL_24AGCTAG_CAG',  # Very small dataset.
        'Exp16_FBL_24AGCTAG_AAC',  # I believe this was also discarded because it was too small.
        'Exp16_hnRNPC_24TGAGTG_AGT', # Unknown why this was blacklisted.
        'Exp16_hnRNPC_24TGAGTG_CAG', # Unknown why this was blacklisted.
        'Exp28_hnRNPC_17CGATTA_AAC', # Unknown why this was blacklisted.
        'Exp31_UBA2_15TGAGTG_AAC',# Too correlated with CDK4 AAC Exp31.
        'Exp31_UBA2_15TGAGTG_CAG', # Empty.
        'Exp31_UBA2_05TGAGTG_CAG',
        'Exp31_CDK4_15GCCATG_AAC', # Too correlated with UBA2 AAC Exp31.
        'Exp31_CDK4_15GCCATG_CAG', # Empty.
        'Exp33_CDK4_17GCCATG_AGT',
        'Exp31_CDK4_05GCCATG_CAG',
        'Exp33_CAPNS6_17CACTGT_TCA', # Empty.
        'Exp61_PCBP1-100P_hpGCCATG_TCA',  # Too small.
        'Exp61_PCBP1-100P_hpGCCATG_AGT',  # Too small.
        ]
    unused_datasets = used_datasets.loc[[x<1 for x in used_datasets['places_used']], :]['long_basename'].tolist()
    #unused_datasets = [x for x in unused_datasets if not re.search('No[_ ]vector', x, re.IGNORECASE)]

    black_list.extend(unused_datasets)

    if name in black_list:
        return False
    if 'STD' in name:
        return False
    return True

geo = geo[[purged(name) for name in geo['basename']]].copy()

rep_n = {}
def add_replicate(name):
    if name in rep_n:
        rep_n[name] += 1
    else:
        rep_n[name] = 1
    return 'easyCLIP {} Replicate {}'.format(name, rep_n[name])

geo['title'] = [add_replicate(protein) for protein in geo['Gene']]
geo['raw file'] = ['{}.fastq.gz'.format(x) for x in geo.basename]
geo['raw file read mate'] = ['{}.fastq.gz'.format(x) for x in geo.basename_r2]

geo_upload_columns += ['raw file read mate']
geo = geo.loc[:, [(x in geo_upload_columns) for x in geo.columns]]
geo = geo.loc[:, geo_upload_columns]

geo.index = geo['Sample name']
geo.to_excel('tables/seq_template_generated_sample_list.xlsx')
    
df = pandas.read_excel('tables/seq_template_generated_sample_list.xlsx')

```


## "Raw files" table, including md5 checksum

```bash
 md5sum geo_upload/*

# Generates this format, which is separated by a space not a tab:
#bb52b9439e0b8aaf32f7cd6edbded693  geo_upload/Exp15_FBL_24AGCTAG_AGT.fastq.gz
#368e7533aaecfa0c9d7939356197c4d6  geo_upload/Exp15_FBL_24AGCTAG_R2_AGT.fastq.gz
```
The output from md5sum was pasted by hand into a text file, given a header and 
 added to the raw files table:

```python
md5 = pandas.read_csv('tables/md5values.txt', sep=' ')
md5['file_checksum'] = md5.index
md5['raw_file'] = [os.path.basename(x) for x in md5.raw_file]
to_md5 = dict(zip(md5['raw_file'], md5['file_checksum']))
```

The rest of the code to make the raw files table, which is
 in prepare_fastq_for_upload.ipynb:
```python
run_code_to_instrument_and_read_len = {
    'rb': {'Instrument': 'Illumina Hiseq 4000', 'read length': 76},
    'hp': {'Instrument': 'Illumina Hiseq 4000', 'read length': 76},
    '05': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125},
    '17': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125},
    '24': {'Instrument': 'Illumina Hiseq 2500', 'read length': 125},
    'pc': {'Instrument': 'Illumina Miseq', 'read length': 76},
    'tc': {'Instrument': 'Illumina Miseq', 'read length': 76},
}

def fname_to_sequencing_run(fname):
    s = os.path.basename(fname).split('_')
    run_code = s[2][:2]
    return run_code_to_instrument_and_read_len[run_code]['Instrument']

def fname_to_read_length(fname):
    s = os.path.basename(fname).split('_')
    run_code = s[2][:2]
    return run_code_to_instrument_and_read_len[run_code]['read length']

df = pandas.read_excel('tables/seq_template_generated_sample_list.xlsx', index_col=0)

rows = []
for row in df.to_dict('records'):
    rows.append({
        'file name': row['raw file'],
        'file type': 'fastq',
        'instrument model': fname_to_sequencing_run(row['raw file']),
        'read length': fname_to_read_length(row['raw file']),
        'single or paired-end': 'paired-end',
        'file checksum': to_md5.get(row['raw file'], '')
    })
    rows.append({
        'file name': row['raw file'],
        'file type': 'fastq',
        'instrument model': fname_to_sequencing_run(row['raw file read mate']),
        'read length': fname_to_read_length(row['raw file']),
        'single or paired-end': 'paired-end',
        'file checksum': to_md5.get(row['raw file read mate'], '')
    })
df = pandas.DataFrame(rows)
df = df.loc[:, ["file name", "file type", "file checksum", "instrument model", "read length", "single or paired-end"]]

df.to_excel('tables/raw_files_table_for_geo_upload.xlsx')
```

