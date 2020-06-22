import pandas, os, sys, re, time, collections, Bio, pprint, random, glob

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from typing import List, Mapping, Union

class currated_set_of_nonredundant_studies():

    def __init__(
        self, fname: str='./all_TCGA_data/currated_set_of_nonredundant_studies_list.txt',
        study_desc_fname: str='cancerLists/tcga_study_ids_and_descriptions.do',
        verbose=True,
        ):
        """
        Top of the file tcga_study_ids_and_descriptions.do:
        cancer_study_id name    description
        paac_jhu_2014   Acinar Cell Carcinoma of the Pancreas (JHU, J Pathol 2014)  Whole exome sequencing of 23 surgically resected pancreatic carcinomas with acinar differentiation and their matched normals.
        
        study_ids are e.g. cancer_study_id
        study_names are e.g. Acinar Cell Carcinoma of the Pancreas (JHU, J Pathol 2014)
        """

        print(f"Getting study categories from {fname}.")

        lines = open(fname).readlines()

        self.studies = set()
        (cat1, cat2, cat3) = ('', '', '')
        self.study_categories = {}

        def clean_name(li):
            cat = re.sub('Deselect All', '', li).rstrip('\n')
            cat = re.sub('CAT\d', '', cat)
            cat = re.sub('[\t ]', '', cat)
            return cat

        for n, li in enumerate(lines):
            if 'CAT1' in li:
                cat1 = clean_name(li)
                cat2, cat3 = '', ''
            if 'CAT2' in li:
                cat2 = clean_name(li)
                cat3 = ''
            if 'CAT3' in li:
                cat3 = clean_name(li)
                
            if ' samples' in li:
                study_name = lines[n-1].rstrip('\n')
                study_name = re.sub('\t', '', study_name)
                study_name = re.sub('  ', '', study_name)
                self.studies.add(study_name)

                # The skin samples we keep separate.
                if cat1 == 'Skin':
                    self.study_categories[study_name] = [cat2, cat2, cat3]
                else:
                    self.study_categories[study_name] = [cat1, cat2, cat3]
                

        self.nonredundant_studies = self.studies
        #print("__", self.study_categories)
        study_desc = pandas.read_csv(study_desc_fname, sep='\t')
        #scc = [x for x, desc in zip(study_desc.cancer_study_id, study_desc.name) \
        #      if True]#re.search('Squamous cell carcinoma', desc, re.IGNORECASE)]# or \
        self.name_to_id = dict(zip(study_desc.name, study_desc.cancer_study_id))
        self.study_id_to_study_name = dict(zip(study_desc.cancer_study_id, study_desc.name))

        with_id = dict([(name, self.name_to_id[name]) for name in self.nonredundant_studies \
                  if name in self.name_to_id])

        without_id = [name for name in self.nonredundant_studies if name not in self.name_to_id]
        
        if verbose:
            print(f"Nonredundant studies with study id found: {len(with_id)}")
            print(f"Nonredundant studies without study id: {len(without_id)}")

        self.nonredundant_study_ids = list(with_id.values())
        self.nonredundant_study_names = list(with_id.keys())

    def get_study_ids(self):
        return self.nonredundant_study_ids

    def get_study_names(self):
        return self.nonredundant_study_names

    def get_name_to_id(self, name: str=None):
        if name is None:
            return self.name_to_id
        return self.name_to_id.get(name, f'No study ID found for study name {name}')

    def get_study_id_to_from_study_name(self, name: str=None):
        return self.name_to_id(self, name=name)

    def get_study_name_from_study_id(self, study_id: str=None) -> str:
        if study_id is None:
            return self.study_id_to_study_name
        if study_id in self.study_id_to_study_name:
            return self.study_id_to_study_name[study_id]

        # Maybe it's a filename:
        _id = study_id
        if '/' in _id:
            _id = os.path.basename(_id)
        _id = _id.split('.')[0]
        _id = re.sub('our_rbps_', '', _id)
        _id = re.sub('all_rbps_', '', _id)
        _id = re.sub('census_', '', _id)

        return self.study_id_to_study_name.get(
            _id, f'No study name found for study ID {_id}')

    def get_study_cat_from_study_name(self, name: str=None):
        if name is None:
            return self.study_categories
        return self.study_categories.get(name, f'No study category found for study name {name}.')

    def get_study_cat_from_study_id(self, study_id):
        study_name = self.get_study_name_from_study_id(study_id=study_id)
        return self.get_study_cat_from_study_name(study_name)


class tcgaConvolve():

    @staticmethod
    def smooth(y: np.array, box_pts: int) -> np.array:
        #https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    def for_each_mut(
        self,
        by_mutation: dict):

        gene_array = {}
        for (gene, mut), num_patients in by_mutation.items():
            pos = re.search('(\d+)', mut)
            if pos is None:
                continue
            if gene not in gene_array:
                gene_array[gene] = np.zeros(pos+1)
            gene_array[gene][pos] += num_patients


class dataLoader():
    """Load TCGA mutation data from folders of .txt files downloaded via their API.
    Put that data into some dicts of (gene, mutation) -> set of patient IDs.
    (And other dicts of similar information).
    """

    def __init__(
        self, 
        study_info: currated_set_of_nonredundant_studies,
        missense_only=False):

        self.study_info = study_info
        self.missense_only = missense_only

        self.studies = []
        self.by_mutation = collections.defaultdict(set)
        self.mutation_to_studies = collections.defaultdict(dict)

        self.by_gene = collections.defaultdict(set)
        self.n_by_gene = collections.defaultdict(int)
        self.n_by_mutation = collections.defaultdict(int)

    def get_total_patient_num(self, study_fnames: list) -> int:

        all_patients = set()
        print("Counting total patients", end='')
        for fname in study_fnames:
            print('.', end='')
            if not os.path.exists(fname):
                continue
            fh = open(fname, 'r')
            got_header = False

            n_lines = 0
            while not got_header:
                li = fh.readline()
                n_lines += 1
                
                if n_lines > 200:
                    print("Giving up on finding patients in {}...".format(fname))
                    got_header = True

                if 'GENE_ID' in li:
                    got_header = True
                    header = li.rstrip('\n').split('\t')
                    header = set(header) - set(['GENE_ID', 'COMMON'])
                    patients = [self.sample_id_to_patient_id(x) for x in header]
                    all_patients |= set(patients)
            fh.close()

        return len(all_patients)

    def wipe_all_studies(self):
        self.studies = []
        self.by_mutation = collections.defaultdict(set)
        self.mutation_to_studies = collections.defaultdict(dict)  

    def do_a_set_of_studies(self, study_fnames: List[str]):
        """
        self.by_mutation[(gene, mut)] = set of patient IDs
        self.mutation_to_studies[(gene, mut)][study_desc] = patient_ids
        """

        self.n_total_patients = self.get_total_patient_num(study_fnames)
        self.add_studies(study_fnames=study_fnames)

    def add_studies(self, study_fnames: List[str]):
        # This sets self.by_mutation, self.studies, and self.mutation_to_studies.
        [self.add_study(fname) for fname in study_fnames if os.path.exists(fname)]
        self.refresh_counters()  

    def refresh_counters(self):
        self.by_gene = collections.defaultdict(set)
        self.n_by_gene = collections.defaultdict(int)
        self.n_by_mutation = collections.defaultdict(int)

        for (gene, mut) in self.by_mutation.keys():
            self.n_by_mutation[(gene, mut)] = len(self.by_mutation[(gene, mut)])
            self.by_gene[gene] |= self.by_mutation[(gene, mut)]

        for gene, set_of_patient_ids in self.by_gene.items():
            self.n_by_gene[gene] = len(set_of_patient_ids)

    def add_study(self, study_fname: str) -> None:
        """
        study_ids are eg cancer_study_id
        study_names are eg Acinar Cell Carcinoma of the Pancreas (JHU, J Pathol 2014)
        study filenames are {somepath}/{prefix}{study_id}.txt

        self.studies: List of pandas.DataFrame objects.
        self.by_mutation: Dict of (gene, mut) -> set of patient IDs
        self.study_info: currated_set_of_nonredundant_studies object
        """

        studies = self.studies
        by_mutation = self.by_mutation

        try:
            studies.append(
                pandas.read_csv(study_fname, sep='\t', comment='#', dtype='str'))#, index_col='GENE_ID'))
        except:
            print(f"Failed with file {study_fname}",)

        if 'GENE_ID' not in studies[-1].columns:
            print(f"No GENE_ID columns for {study_fname}. Skipped.")
            return studies, by_mutation
        
        studies[-1] = studies[-1].loc[[bool(x!='GENE_ID') for x in studies[-1].GENE_ID], :]
        df = studies[-1]

        # Format: GENE_ID  COMMON Patient1
        patient_cols = [x for x in df.columns if (x not in ['GENE_ID', 'COMMON'])]

        by_mutation_this_study = collections.defaultdict(set)

        for n in range(0, len(df.index)):  # For each protein.
            gene = df.iloc[n]['COMMON']

            # Are there any mutations?
            mutations = [x for x in df.iloc[n][patient_cols] if type(x) == type('') and x != 'NaN']
            
            if self.missense_only:
                mutations = [x for x in mutations if self.is_missense(x)]
            
            # If there are mutations:
            if len(mutations):
                if random.randint(0, 2000) == 1:
                    print(f"study {len(studies)}, {gene}->{mutations})")
                    
                for patient, muts in zip(patient_cols, df.iloc[n][patient_cols]):

                    if type(muts) != type(''):
                        continue

                    for mut in muts.split(','):
                        # For each mutation in this gene in this patient:

                        if self.missense_only and not self.is_missense(mut):
                            continue

                        # Record the patient ID for this (gene, mutation) pair.
                        self.by_mutation[(gene, mut)].add(
                            self.sample_id_to_patient_id(patient))

                        by_mutation_this_study[(gene, mut)].add(
                            self.sample_id_to_patient_id(patient))

        for (gene, mut), patient_ids in by_mutation_this_study.items():
            if self.study_info is None:  # Use the study ID.
                self.mutation_to_studies[(gene, mut)][study_fname] = patient_ids

            else:  # Use the study name.
                study_name = self.study_info.get_study_name_from_study_id(study_fname)
                if study_name in self.mutation_to_studies[(gene, mut)]:
                    self.mutation_to_studies[(gene, mut)][study_name] |= patient_ids
                else:
                    self.mutation_to_studies[(gene, mut)][study_name] = patient_ids

    @staticmethod
    def sample_id_to_patient_id(sample_id):
        if sample_id.count('-') <= 1:
            return sample_id
        if 'TCGA' in sample_id:
            return '-'.join(sample_id.split('-')[:3])
        if sample_id.count('-') > 1:
            return '-'.join(sample_id.split('-')[:2])

    @staticmethod
    def get_nonrbps():
        nonrbps = [] #['CDKN2A', 'HIST1H1C', 'HLA-A', 'HLA-B']
        return nonrbps

    @classmethod
    def remove_non_rbps(cls, by_mutation):
        nonrbps = cls.get_nonrbps()
        to_del = set()
        for k in by_mutation.keys():
            if k[0] in nonrbps:
                to_del.add(k)
        for k in to_del:
            del by_mutation[k]
        return by_mutation

    @staticmethod
    def is_missense(mut):
        if ('*' in mut) or ('?' in mut) or ('splice' in mut):
            return False
        if not re.match('\w+(\d+)\w+', mut):
            return False        
        return True

    @staticmethod
    def get_study_info(study_fname, study_info):
        if '/' in study_fname:
            study_fname = os.path.basename(study_fname)
        study_fname = study_fname.split('.')[0]
        study_name = re.sub('our_rbps_', '', study_fname)

        if study_name in study_info:
            #print(study_name, study_info.loc[study_name]['name'])
            return study_info[study_name]
        else:
            print("Couldn't find study id for", study_fname)
            return 'Unknown'
            
    def map_mutations_to_the_number_of_patients_by_study(self) -> None:
        
        self.mutation_to_studies_n = collections.defaultdict(dict)
        self.by_study_cat1_n = collections.defaultdict(dict)

        # ONLY FOR MSK-IMPACT Clinical Sequencing Cohort (MSKCC, Nat Med 2017):
        # This is mixed, so patient IDs have to be looked up.
        patient_data_fname = "all_TCGA_data/msk_impact_2017_clinical_data.tsv"
        df = pandas.read_csv(patient_data_fname, sep='\t')
        df['Cancer Type'] = [x.replace(' ', '') for x in df['Cancer Type']]

        match_types = {
            'PancreaticCancer': 'Pancreas',
            'ColorectalCancer': 'Bowel',
            'Non-SmallCellLungCancer': 'Lung',
            'BladderCancer': 'Bladder/Urinary',
            'EndometrialCancer': 'Uterus',
            'ProstateCancer': 'Prostate',
            'GastrointestinalNeuroendocrineTumor': 'Other',
            'SoftTissueSarcoma': 'Soft Tissue',
            'HepatobiliaryCancer': 'Liver',
            'BreastCancer': 'Breast',
            'SalivaryGlandCancer': 'Other',
            'MatureB-CellNeoplasms': 'Lymphoid',
            'CancerOfUnknownPrimary': 'Other', 
            'Melanoma': 'Melanoma',
            'Non-HodgkinLymphoma': 'Lymphoid',
            'RenalCellCarcinoma': 'Kidney',
        }

        df['Cancer Type'] = [match_types.get(x, x) for x in df['Cancer Type']]
        patient_to_cancer_type = dict(zip(df['Patient ID'], df['Cancer Type']))

        # 1. Make mutation_to_studies_n dict of {(protein, mutation): {study_name: number of patients with mutation}}
        # 2. Make mutation_to_studies_n dict of {(protein, mutation): {study category: number of patients with mutation}}
        for mut, studies in self.mutation_to_studies.items():

            for study_name in studies:
                #print(all_dl.mutation_to_studies[mut])
                self.mutation_to_studies_n[mut][study_name] = len(
                    self.mutation_to_studies[mut][study_name])

                #study_id = study_fname.split('/')[-1].split('.txt')[0]

                study_cat = self.study_info.get_study_cat_from_study_name(study_name)

                # This study is of mixed cancers, so it is assigned by patient ID.
                if study_name == "MSK-IMPACT Clinical Sequencing Cohort (MSKCC, Nat Med 2017)":

                    for patient_id in self.mutation_to_studies[mut][study_name]:
                        cancer_type = patient_to_cancer_type[patient_id]
                        if cancer_type not in self.by_study_cat1_n:
                            self.by_study_cat1_n[mut][cancer_type] = 1
                        else:
                            self.by_study_cat1_n[mut][cancer_type] += 1

                    continue

                if type(study_cat) == type([]):
                    study_cat = study_cat[0]

                if study_cat not in self.by_study_cat1_n[mut]:
                    self.by_study_cat1_n[mut][study_cat]= len(
                        self.mutation_to_studies[mut][study_name])
                else:
                    self.by_study_cat1_n[mut][study_cat] += len(
                        self.mutation_to_studies[mut][study_name])


