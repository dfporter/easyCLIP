import pandas, re, os, collections, random

class scheme():
    
    def __init__(self, scheme_fname=''):
        if scheme_fname != '':
            self.read_scheme(scheme_fname=scheme_fname)

    def __str__(self):
        li = "sameRiver.scheme.scheme object with the following data frame (self.scheme_df):"
        li += str(self.scheme_df)
        return li
    
    def fill_in_info_column_by_filename(self, col_name, dict_by_fname):
        basenames = [os.path.splitext(os.path.basename(x))[0] for x in self.scheme_df['long_fname_R1']]
        self.scheme_df[col_name] = [
            dict_by_fname.get(x, '') for x in basenames]
        self.generate_mappings()

    def fill_in_info_column_by_p6(self, col_name, p6_info):
        for p6, v in p6_info.items():
            print(p6, 'aa')
            for p3 in ['TCA', 'AGT', 'AAC', 'CAG']:
                p6p3 = '_'.join([str(p6), str(p3)])
                if p6p3 in self.p6p3_string_to_scheme_dict:
                    self.p6p3_string_to_scheme_dict[p6p3].update({col_name: v})
        
        
        self.scheme_df = pandas.DataFrame.from_dict(self.p6p3_string_to_scheme_dict, orient='index')
        self.scheme_df.index = range(len(self.scheme_df.index))
        self.generate_mappings()
        
    def fill_in_info_column_by_p6p3(self, col_name, p6p3_info):
        for p6p3, v in p6p3_info.items():
            if type(p6p3) != type(''):
                p6p3 = '_'.join([str(p6p3[0]), str(p6p3[1])])
                
            if p6p3 in self.p6p3_string_to_scheme_dict:
                self.p6p3_string_to_scheme_dict[p6p3].update({col_name: v})
        
        
        self.scheme_df = pandas.DataFrame.from_dict(self.p6p3_string_to_scheme_dict, orient='index')
        self.scheme_df.index = range(len(self.scheme_df.index))
        self.generate_mappings()
        
    def generate_mappings(self):
        df = self.scheme_df
        
        self.p6p3_to_gene_replicate = {}
        self.p6p3_string_to_gene_replicate = {}
        
        for exp, p6, p3, gene, replicate in zip(
            df.Experiment, df.P6_BC, df.P3_BC, df.Gene, df.Replicate):
            self.p6p3_to_gene_replicate[(exp, p6, p3)] = (gene, replicate)
            self.p6p3_string_to_gene_replicate['_'.join([str(p6), str(p3)])] = (gene, replicate)

        self.p6p3_to_scheme_dict = {}
        self.p6p3_string_to_scheme_dict = {}
        
        for row in df.to_dict('records'):
            (exp, p6, p3) = (row['Experiment'], row['P6_BC'], row['P3_BC'])
            self.p6p3_to_scheme_dict[(exp, p6, p3)] = row
            self.p6p3_string_to_scheme_dict['_'.join([str(p6), str(p3)])] = row
            
        self.p6_to_gene = collections.defaultdict(set)
        self.p6_to_list_of_dicts = collections.defaultdict(list)
        
        for row in df.to_dict('records'):
            self.p6_to_gene[row['P6_BC']] |= set([row['Gene']])
            self.p6_to_list_of_dicts[row['P6_BC']].append(row)
            
    def read_scheme(self, scheme_fname=''):
        self.scheme_df = pandas.read_excel(scheme_fname, index=False)

        # We often add a two letter prefix to the P6 barcode, using non-DNA letters.
        # We chop that off here.
        self.scheme_df['P6_DNA'] = [x.upper() for x in self.scheme_df['P6_BC']]
        self.scheme_df['P6_DNA'] = [''.join([letter for letter in barcode if (
            letter in ['A', 'T', 'C', 'G'])]) for barcode in self.scheme_df['P6_DNA']]
        
        if 'Gene' in self.scheme_df.columns:
            self.scheme_df['Gene'] = [x.replace(' ', ':') for x in self.scheme_df['Gene']]
        if 'Protein' in self.scheme_df.columns:
            self.scheme_df['Protein'] = [x.replace(' ', ':') for x in self.scheme_df['Protein']]

        if 'Percent XL' not in self.scheme_df.columns:
            self.scheme_df['Percent XL'] = 0
        
        self.generate_mappings()
        self.set_long_form_filenames()

        if 'Gene' in self.scheme_df.columns:
            self.proteins = set(self.scheme_df['Gene'].tolist())
        if 'Protein' in self.scheme_df.columns:
            self.proteins = set(self.scheme_df['Proetin'].tolist())

    def p6p3s(self, as_type='str'):
        if as_type == 'str':
            return list(self.p6p3_string_to_gene_replicate.keys())
        
    def gene_from_fname(self, fname):

        # First just look for something we recognize.
        longest_match = ''
        fname = re.sub(':', '-', fname)
        for gene in self.proteins:
            if re.sub(':', '-', gene) in fname:
                if len(gene) > len(longest_match):
                    longest_match = gene

        if len(longest_match) > 1:
            return longest_match

        # Failing that, guess by position.
        s = os.path.basename(fname).split('.')[0].split('_')
        if len(s) < 2:
            return ''

        _str = '_'.join(s[:2])
        if _str in self.p6p3_string_to_gene_replicate:
            return self.p6p3_string_to_gene_replicate[_str][0]

    def fname_from_gene(self, gene):
        df = self.scheme_df

        if 'Gene' in self.scheme_df.columns:
            s = df[df['Gene']==gene].copy()
        if 'Protein' in self.scheme_df.columns:
            s = df[df['Protein']==gene].copy()

        if len(s.index) == 0:
            print("Looked for a filename for {} but could not find it in the scheme file.".format(gene))
            return ['']

        print(s)

    def p6p3_from_filename(self, bed_fname):
        basename = os.path.basename(bed_fname).split('.')[0]
        sp = basename.split('_')
        return '_'.join(sp[:2])
    
    def percent_xl_of_fname(self, bed_fname):
        p6p3 = self.p6p3_from_filename(bed_fname)
        if p6p3 not in self.p6p3_string_to_scheme_dict:
            return ''
        return self.p6p3_string_to_scheme_dict[p6p3]['Percent XL']

    def percent_xl_of_p6p3(self, p6p3):
        if p6p3 not in self.p6p3_string_to_scheme_dict:
            return ''
        return self.p6p3_string_to_scheme_dict[p6p3]['Percent XL']
    
    def percent_xl_of_protein_from_file(
        self, protein,
        xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx',
        reload_file=False):

        if reload_file or not(hasattr(self, 'xl_by_prot_from_file')):
            xl = pandas.read_excel(xl_rate_fname)
            xl = xl[xl['Label']=='% XL (minimal region)'].copy()
            self.xl_by_prot_from_file = xl.groupby('Protein').mean().to_dict()['Value']

            for prot in list(self.xl_by_prot_from_file.keys()):
                self.xl_by_prot_from_file[re.sub(' ', '', prot)] = self.xl_by_prot_from_file[prot]
                self.xl_by_prot_from_file[re.sub(':', '-', prot)] = self.xl_by_prot_from_file[prot]

            print('loaded {}:'.format(xl_rate_fname))
            #print(m)
#            print(self.xl_by_prot_from_file)

        xl = self.xl_by_prot_from_file

        # Always use capitals.
        protein = protein.upper()

        for col in list(xl.keys()):
            xl[col.upper()] = xl[col]

        try:
            xl['PCBP1:100P'] = xl['PCBP1 L100P']
            xl['PCBP1:100Q'] = xl['PCBP1 L100Q']
            xl['PCBP1:dKH'] = xl['PCBP1 ∆KH2']
        except:
            pass

        if protein in self.xl_by_prot_from_file:
            return self.xl_by_prot_from_file[protein]
        
        else:
            _protein = self.gene_from_fname(protein)

            if _protein in self.xl_by_prot_from_file:
                return self.xl_by_prot_from_file[_protein]

            try:
                gene = protein.split('_')[1]
                if gene in self.xl_by_prot_from_file:
                    return self.xl_by_prot_from_file[gene]
                else:
                    return 0.
            except:
                return 0.
    
    def set_long_form_filenames(self):
        scheme_df = self.scheme_df
        
        scheme_df['long_fname_R1'] = [
        '{exp}_{gene}_{p6}_{p3}.fastq'.format(
            exp=exp, gene=re.sub(':', '-', gene), p6=p6, p3=p3) for (exp, gene, p6, p3) in zip(
            scheme_df['Experiment'], scheme_df['Gene'], scheme_df['P6_BC'], scheme_df['P3_BC'])
        ]
        scheme_df['long_fname_R2'] = [
            '{exp}_{gene}_{p6}_R2_{p3}.fastq'.format(
                exp=exp, gene=re.sub(':', '-', gene), p6=p6, p3=p3) for (exp, gene, p6, p3) in zip(
                    scheme_df['Experiment'], scheme_df['Gene'], scheme_df['P6_BC'], scheme_df['P3_BC'])
        ]
        self.long_fname_to_info = {}
        self.long_basename_to_info = {}
        for row in scheme_df.to_dict('records'):
            self.long_fname_to_info[row['long_fname_R1']] = row
            self.long_fname_to_info[row['long_fname_R2']] = row
            self.long_basename_to_info[os.path.splitext(row['long_fname_R1'])[-2]] = row

        self.p6p3_to_long_filename_r1 = {}
        self.p6p3_to_long_filename_r2 = {}

        for (p6, p6_dna, p3, long_fname_R1, long_fname_R2) in zip(
            scheme_df['P6_BC'], scheme_df['P6_DNA'], scheme_df['P3_BC'],
            scheme_df['long_fname_R1'], scheme_df['long_fname_R2']):

            self.p6p3_to_long_filename_r1[(p6, p3)] = long_fname_R1
            self.p6p3_to_long_filename_r1[(p6_dna, p3)] = long_fname_R1
            self.p6p3_to_long_filename_r2[(p6, p3)] = long_fname_R2
            self.p6p3_to_long_filename_r2[(p6_dna, p3)] = long_fname_R2

