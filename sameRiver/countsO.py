import HTSeq, collections, pandas, os, re, dill, importlib
import numpy as np

import sameRiver
import sameRiver.scheme
import sameRiver.biotypeUtils
importlib.reload(sameRiver.scheme)
importlib.reload(sameRiver.biotypeUtils)


class countsO():
    
    def __init__(
        self, filename=None, scheme=None, indexing_col_name='gene_name', index_col=None,
        exclude_unmapped=True, ignore_biotypes=[]):
        """The counts files passed to countsO (via filename) are expected to be raw read counts,
        without normalization. counts.txt and ann_counts.txt files are always raw read counts.

        load_excel() sets: 
            (1) the pandas.DataFrame self.raw_counts_df
            (2) and the dict self.raw_counts.

        self.counts_per_read holds reads per million:
            key=RNA, value: {column in input file->reads per million}
        self.raw_counts:
            key=RNA, value: {column in input file->raw reads}

        self.counts_per_million_df:
            dataframe version of counts_per_read and per_million_by_protein.

        self.per_million_by_protein:
            key=RNA, value: {protein->reads per million}
        self.raw_counts_by_protein:
            key=RNA, value: {protein->raw reads}
        self.counts_per_protein:
            key=RNA, value: {protein->crosslinks per protein molecules}

        """
        self.log = f"Instantiated with {locals()}\n"
        os.makedirs('./data/', exist_ok=True)

        self.file_dfs = {}
        self.dfs = {}
        self.schemes = {}
        self.indexing_col_name = indexing_col_name
        self.filename = filename

        if filename is not None:
            self.load_excel(filename, scheme=scheme, indexing_col_name=indexing_col_name,
                           index_col=index_col, exclude_unmapped=exclude_unmapped,
                           ignore_biotypes=ignore_biotypes)

                    
    def set_blacklist(self, black_list=None):
        if black_list is None:
            self.black_list = [
                'Exp16_FBL_24AGCTAG_CAG',  # Very small dataset.
                'Exp16_hnRNPC_24TGAGTG_AGT',
                'Exp16_hnRNPC_24TGAGTG_CAG',
                'Exp28_hnRNPC_17CGATTA_AAC',
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
                'YBX',
                'AURKA'
            ]
        else:
            self.black_list = black_list
            
    def blacklist(self):
        return self.black_list
    
    def in_blacklist(self, col_name):
        for col in self.black_list:
            if re.search(col, col_name):
                return True
        return False
    
    def apply_blacklist(self, black_list=None):
        """Removes columns.
        Resets self.raw_counts_df and self.raw_counts. These are initially set by
        self.load_excel.
        """

        print("Applying blacklist ({} columns)...".format(len(self.raw_counts_df.columns)))
        
        if black_list is None:
            
            if not hasattr(self, 'black_list'):
                self.set_blacklist()
            black_list = self.black_list
        
        if type(black_list) == type(''):
            black_list = [black_list]
            
        for blackl_pat in black_list:
            
            if blackl_pat in self.raw_counts_df.columns:
                del self.raw_counts_df[blackl_pat]
                
            matches = [x for x in self.raw_counts_df.columns if re.search(blackl_pat, x)]

            for match in matches:
                
                if match == 'gene_name':
                    print("Asked to blacklist the column gene_name. Refusing.")
                    continue
                    
                del self.raw_counts_df[match]
                
        self.raw_counts = {}
        for row in self.raw_counts_df.to_dict('records'):
            self.raw_counts[row[self.indexing_col_name]] = row
        
        print("...{} columns after blacklisting.".format(len(self.raw_counts_df.columns)))
    
    @staticmethod
    def vprint(verbose, *args):
        if verbose:
            if len(args) == 1:
                print(args[0])
            else:
                print(args)
        
    def edit(self, black_list=None, verbose=False, save=False,
        save_to=None, xl_rate_fname='percentCrosslinked.xlsx',
        mapping_fname='enst_transcript_id_name_biotype_map.txt',
        include_percent_xl=True):

        input_file_had_biotypes = False
        if 'Gene type' in self.raw_counts_df.columns:
            name_to_type = dict(zip(self.raw_counts_df.gene_name, self.raw_counts_df['Gene type']))
            input_file_had_biotypes = True

        self.vprint(verbose, 'counts.edit(), input slice:')
        self.vprint(verbose, self.raw_counts_df.iloc[:2, :2])
        self.simplify_column_names(verbose=verbose)#include_gene=True, include_exp=True, verbose=verbose)
        self.vprint(verbose, 'counts.edit(), simplified column names:')
        self.vprint(verbose, self.raw_counts_df.iloc[:2, :2])
        self.set_blacklist(black_list=black_list)
        self.apply_blacklist()  # Resets raw_counts, raw_counts_df.

        self.vprint(verbose, "counts.edit(), purging/filling in counts...")
        self.purge_counts()
        self.fill_in_counts()
        print('counts.edit(), purged, blacklisted, and filled in (still raw read counts):')
        print(self.raw_counts_df.iloc[:2, :2])

        self.set_counts_by_protein()

        if input_file_had_biotypes:
            self.raw_counts_df['Gene type'] = [name_to_type[x] for x in self.raw_counts_df.gene_name]
        else:
            self.raw_counts_df = self.add_biotypes_column(
                df=self.raw_counts_df, mapping_fname=mapping_fname, reload=True)

        # Update the dict.
        self.raw_counts = {}
        for row in self.raw_counts_df.to_dict('records'):
            self.raw_counts[row[self.indexing_col_name]] = row
        
        if save or (save_to is not None):
            if save_to is None:
                save_to = self.filename
            try:
                self.raw_counts_df.to_csv(save_to, sep='\t')
                print("Saved edited counts to {}.".format(save_to))
            except:
                print("Asked to save edited counts but failed. Tried to save to {}".format(save_to))
    
    def simplify_column_names(
        self, include_gene=False, include_exp=False, discard_unknown=False, verbose=False):

        print("Simplifying column names...")
        new_cols = []
        to_discard = []

        print("Using these values to look up column names: {}".format(
            self.scheme.long_basename_to_info.keys()))
        for col in self.raw_counts_df.columns:

            basename = os.path.splitext(os.path.basename(col))[0].rstrip('_+').rstrip('_-')

            if basename in self.scheme.long_basename_to_info:
                info = self.scheme.long_basename_to_info[basename]
                new_cols.append(basename)
                continue

            else:
                print("Failed to find column {} in scheme. Looked for somthing matching {}.".format(
                    col, basename))
                new_cols.append(col)
                to_discard.append(col)
                continue

            p6p3 = self.scheme.p6p3_from_filename(col)
            
            if p6p3:
                new_cols.append(p6p3)
                if include_gene and (p6p3 in self.scheme.p6p3_string_to_scheme_dict):
                    new_cols[-1] += '_' + self.scheme.p6p3_string_to_scheme_dict[p6p3]['Gene']
                if include_exp and (p6p3 in self.scheme.p6p3_string_to_scheme_dict):
                    new_cols[-1] += '_' + self.scheme.p6p3_string_to_scheme_dict[p6p3]['Experiment']
                if discard_unknown and include_gene and (p6p3 not in self.scheme.p6p3_string_to_scheme_dict):
                    to_discard.append(new_cols[-1])
            else:
                new_cols.append(col)
                to_discard.append(col)

        print('Old cols', self.raw_counts_df.columns)
        print('new cols', new_cols)
        self.raw_counts_df.columns = new_cols
        for col in set(to_discard) - set(['gene_name']):
            del self.raw_counts_df[col]

        self.raw_counts = {}
        for row in self.raw_counts_df.to_dict('records'):
            self.raw_counts[row[self.indexing_col_name]] = row

        self.vprint(verbose, "Simplified column names:", list(self.raw_counts_df.columns))
    
    def load_excel(self, filename, scheme=None, indexing_col_name=None, index_col=0,
                  exclude_unmapped=True, ignore_biotypes=[]):
        """Set self.raw_counts_df and self.raw_counts.
        self.raw_counts is key: gene name, value: dict of dataframe row (protein->raw counts)
        """

        print(f"Loading counts file {filename}")
        
        if filename.lower().endswith(('xls', 'xlsx')):
            self.raw_counts_df = pandas.read_excel(filename, index_col=index_col)
        else:
            self.raw_counts_df = pandas.read_csv(filename, sep='\t', index_col=index_col)
        
        if exclude_unmapped:
            if '_no_feature' in self.raw_counts_df.index:
                self.raw_counts_df.drop(['_no_feature'], axis=0, inplace=True)

            if '_ambiguous' in self.raw_counts_df.index:
                self.raw_counts_df.drop(['_ambiguous'], axis=0, inplace=True)
        
        if scheme is not None:
            self.scheme = scheme
        else:
            self.scheme = sameRiver.scheme.scheme() 
        
        #print(f"Indexing column={indexing_col_name}")
        
        if (indexing_col_name is not None) and (index_col is not None):
            self.raw_counts_df[indexing_col_name] = self.raw_counts_df.index.tolist()
        
        if (indexing_col_name is None):
            if 'gene_name' in self.raw_counts_df.columns:
                indexing_col_name = 'gene_name'
            elif 'transcript_id' in self.raw_counts_df.columns:
                indexing_col_name = 'transcript_id'
            self.indexing_col_name = indexing_col_name
        
        #print(f"Indexing column={indexing_col_name}")

        print(f'Indexing column ({indexing_col_name}): {self.raw_counts_df[indexing_col_name].head()}')
        
        if 'Gene type' in self.raw_counts_df.columns:
            print("Removing these biotypes before any processing: {}".format(ignore_biotypes))
            init_len = len(self.raw_counts_df.index)
            self.raw_counts_df = self.raw_counts_df[[
                (x not in ignore_biotypes) for x in self.raw_counts_df['Gene type'] ]].copy()
            print("After biotype removal, went from {} to {} RNAs ".format(
                init_len, len(self.raw_counts_df.index)) + \
                  "(introns and exons counted as separate RNAs.)")

        self.raw_counts = {}
        for row in self.raw_counts_df.to_dict('records'):
            self.raw_counts[row[indexing_col_name]] = row
       
    def fill_in_counts(self):
    
        self.raw_counts_df.replace(to_replace=np.nan, value=0., inplace=True)

        bed_names = set()
        for dict_of_beds in self.raw_counts.values():
            bed_names |= set(dict_of_beds.keys())

        for gene_type, dict_of_beds in self.raw_counts.items():
            for bed_fname in bed_names:
                if (bed_fname not in dict_of_beds) or np.isnan(self.raw_counts[gene_type][bed_fname]):
                    self.raw_counts[gene_type][bed_fname] = 0.
    
    
    def set_counts_by_protein(self):
        """Make dicts ordered by protein.
        """
    
        self.raw_counts_by_protein = collections.defaultdict(dict)
        
        for gene_type, dict_of_beds in self.raw_counts.items():
            for bed_fname, value in dict_of_beds.items():
                
                protein = self.scheme.gene_from_fname(bed_fname)
                
                self.raw_counts_by_protein[gene_type].setdefault(protein, [])
                self.raw_counts_by_protein[gene_type][protein].append(value)

    def reorder_by_protein(self, _counts=None):
        
        if _counts is None:
            _counts = self.raw_counts

        _counts_by_protein = collections.defaultdict(dict)

        for gene_type, dict_of_beds in _counts.items():
            for bed_fname, value in dict_of_beds.items():
                
                protein = self.scheme.gene_from_fname(bed_fname)
                
                _counts_by_protein[gene_type].setdefault(protein, [])
                _counts_by_protein[gene_type][protein].append(value)

        return _counts_by_protein
    
    def by_rank(self, _counts):
        _by_rank = collections.defaultdict(dict)

        bed_names = set()
        for dict_of_beds in _counts.values():
            bed_names |= set(dict_of_beds.keys())

        max_rank = len(bed_names)
        for gene_type, dict_of_beds in _counts.items():
            ranked = sorted(dict_of_beds.keys(), key=lambda x: dict_of_beds[x], reverse=True)

            for n, bed_fname in enumerate(ranked):
                _by_rank[gene_type][bed_fname] = n

            for bed_fname in bed_names:
                if bed_fname not in _counts[gene_type]:
                    _by_rank[gene_type][bed_fname] = max_rank

        return _by_rank

    def counts_as_rank(self):

        self.counts_per_protein_ranked = self.by_rank(self.counts_per_protein)
        self.counts_per_read_ranked = self.by_rank(self.counts_per_read)
        
    def purge_counts(self):
        """Remove columns we can't match to proteins in the scheme."""
        purged = collections.defaultdict(dict)

        for gene_type, dict_of_beds in self.raw_counts.items():
            for bed_fname, num_counts in dict_of_beds.items():
        
                gene_name = self.scheme.gene_from_fname(bed_fname)
                
                if (gene_name is None) or (not gene_name) or (str(gene_name) == 'None'):
                    continue
                
                purged[gene_type][bed_fname] = self.raw_counts[gene_type][bed_fname]
        
        self.raw_counts = purged
        
    def counts_cols_only(self, _dict=None, _df=None):
    
        if _dict is not None:
            dfy = pandas.DataFrame.from_dict(_dict, orient='index')
        elif _df is not None:
            dfy = _df
        else:
            print("Pass either dict of dicts or DataFrame.")
            return pandas.DataFrame()

        out_cols = [x for x in dfy if self.scheme.gene_from_fname(x)]
        out_df = dfy[out_cols].copy()

        try:
            out_df.index = dfy.transcript_id
        except:
            pass

        return out_df

    def add_biotypes_column(
        self, mapping_fname='enst_transcript_id_name_biotype_map.txt', #'/opt/genome/Homo_sapiens.GRCh38.83/enst_transcript_id_name_biotype_map.txt',
        df=None, lookup_column='gene_name', lookup_column_in_mapping_file=None, reload=False,
        gtf_filename='combined.gtf'):
        """
        Expect a mapping file with this format:
        transcript_id	gene_name	transcript_biotype
        ENST00000456328	DDX11L1	processed_transcript
        ENST00000450305	DDX11L1	transcribed_unprocessed_pseudogene
        ENST00000488147	WASH7P	unprocessed_pseudogene
        ENST00000619216	MIR6859-1	miRNA
        ENST00000473358	RP11-34P13.3	lincRNA
        ENST00000469289	RP11-34P13.3	lincRNA
        ENST00000607096	MIR1302-2	miRNA
        ENST00000417324	FAM138A	lincRNA
        ENST00000461467	FAM138A	lincRNA"""
        
        if not os.path.exists(mapping_fname):
            
            sameRiver.biotypeUtils.make_enst_gene_name_biotype_map_file(
                gtf_filename=gtf_filename,
                out_filename='enst_transcript_id_name_biotype_map.txt')
            mapping_fname = 'enst_transcript_id_name_biotype_map.txt'

        if (not hasattr(self, 'to_biotypes')) or reload:
            to_biotypes_df = pandas.read_csv(mapping_fname, sep='\t')

            #print("Adding biotypes column. Biotypes mapping file {0} describes {1} biotypes:\n{2}".format(
            #    mapping_fname, len(set(to_biotypes_df.transcript_biotype)),
            #    to_biotypes_df.transcript_biotype.value_counts()))
            
            if lookup_column_in_mapping_file is None:
                lookup_column_in_mapping_file = [
                    x for x in ['Gene name', 'gene_name', 'transcript_id'] if x in to_biotypes_df.columns]

            type_col = [x for x in ['Biotype', 'transcript_biotype'] if x in to_biotypes_df.columns][0]

            to_biotypes = {}
            for a_col in lookup_column_in_mapping_file:
                to_biotypes.update(
                    dict(zip(to_biotypes_df[a_col].tolist(), 
                        to_biotypes_df[type_col].tolist())))

        print("Adding biotypes.\n\n")
        #print(to_biotypes['PCBP1'])
        with open('temp.txt', 'w') as f:
            for k,v in to_biotypes.items():
                f.write('{}\t{}\n'.format(k, v))

        def fix_snhg(name, biotype):
            if re.search('\ASNHG\d*::intron', name):
                return 'snoRNA'
            return biotype
            
        if df is None:
            self.raw_counts_df['Gene type'] = [
                fix_snhg(x, to_biotypes.get(x.split('::')[0], 'Unknown')) for x in self.raw_counts_df[lookup_column].tolist()]
            
        else:
            df['Gene type'] = [
                fix_snhg(x, to_biotypes.get(x.split('::')[0], 'Unknown')) for x in df[lookup_column].tolist()]
            #print("Added biotypes:", df['Gene type'])

        return df
    
    def proteins(self):
        _proteins = set()

        for col in self.raw_counts_df.columns:
            _proteins.add(self.scheme.gene_from_fname(col))

        _proteins -= set(['', None])

        return _proteins
    
    def cols_of_protein(self, protein):
        return [col for col in self.raw_counts_df.columns if re.search(protein, col)]
    
    def only_columns_with_known_proteins(self, output_fname=None):
        
        def keep(col):
            if (len(col.split('_'))>2):
                return True
            elif re.search('Gene', col, re.IGNORECASE) or re.search('type', col, re.IGNORECASE) or re.search('[Nn]ame', col):
                return True
            elif len(col.split('_'))==2:
                return False
            return True
        
        self.raw_counts_df = self.raw_counts_df[[x for x in self.raw_counts_df.columns if keep(x)]].copy()
        
        if output_fname is not None:
            self.write_to_file(output_fname)
        
    def write_to_file(self, fname):
        print("Writing raw counts dataframe to {}.".format(fname))
        self.raw_counts_df.to_csv(fname, sep='\t')

    def save(self):
        print("Saving to data/counts.dill...")
        with open('data/counts.dill', 'wb') as f:
            dill.dump(self, f)
        print("...Saved.")

    @staticmethod
    def load():
        print("Loading {}...".format('data/counts.dill'))
        with open('data/counts.dill', 'rb') as f:
            counts = dill.load(f)
        print("...Loaded.")
        return counts

    @staticmethod
    def find_hits_columns(df, verbose=False):

        hits_columns = [col for col in df.columns if ((len(col.split('_'))>2) and (col != 'gene_name'))]
                
        if verbose:
            print("Found {0} HITS columns out of {1} total columns (had '_' pattern).".format(
                len(hits_columns), len(df.columns)))
            
        return hits_columns

    @staticmethod
    def drop_odds(df):
        if '_ambiguous' in df.index:
            df.drop('_ambiguous', inplace=True)
        if '_no_feature' in df.index:
            df.drop('_no_feature', inplace=True)    