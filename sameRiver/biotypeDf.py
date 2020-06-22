import pandas, re, collections, os
import sameRiver
import sameRiver.biotypeUtils
import sameRiver.biotyper

class biotypeDf():

    def __init__(
        self, df=None, gene_col='gene_name',
        value_col=None,
        type_col='Gene type'):

        self.df = df.copy()
        (self.gene_col, self.value_col, self.type_col) = (
            gene_col, value_col, type_col)

        if value_col is None:
            self.value_col = self.numeric_columns()

    def numeric_columns(self):
        return [x for x in self.df.columns if (
            self.df[x].dtype.kind in 'bifc')]

    #def add_biotype(self):

        # Sets df['Gene type']
#        sameRiver.countsO.countsO().add_biotypes_column(
#            df=self.df, lookup_column=self.gene_col,
#            mapping_fname='/opt/genome/Homo_sapiens.GRCh38.83/gene_name_to_biotype.txt')
        
#        return self.df

    def add_biotypes_column(
        self, mapping_fname='enst_transcript_id_name_biotype_map.txt', #'/opt/genome/Homo_sapiens.GRCh38.83/enst_transcript_id_name_biotype_map.txt',
        df=None, lookup_column='gene_name', lookup_column_in_mapping_file=None,
        reload=False,
        gtf_filename='combined.gtf'):
        """
        Expect a mapping file with this format:
        transcript_id   gene_name   transcript_biotype
        ENST00000456328 DDX11L1 processed_transcript
        ENST00000450305 DDX11L1 transcribed_unprocessed_pseudogene
        ENST00000488147 WASH7P  unprocessed_pseudogene
        ENST00000619216 MIR6859-1   miRNA
        ENST00000473358 RP11-34P13.3    lincRNA
        ENST00000469289 RP11-34P13.3    lincRNA
        ENST00000607096 MIR1302-2   miRNA
        ENST00000417324 FAM138A lincRNA
        ENST00000461467 FAM138A lincRNA"""
        
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

        self.to_biotypes = to_biotypes

        print("Adding biotypes.\n\n")
        with open('temp.txt', 'w') as f:
            for k,v in to_biotypes.items():
                f.write('{}\t{}\n'.format(k, v))

        def fix_snhg(name, biotype):
            if re.search('\ASNHG\d*::intron', name):
                return 'snoRNA'
            return biotype
        
        df['Gene type'] = [
            fix_snhg(x, to_biotypes.get(x.split('::')[0], 'Unknown')) for x in df[lookup_column].tolist()]

        return df

    def only_biotypes_with_enough_reads(self, minimum_reads=5000):
        ok_biotypes = sameRiver.biotyper.biotyper.biotypes_with_enough_reads(
        	self.df, cutoff=minimum_reads)
        self.df = self.df[[
            (x in ok_biotypes) for x in self.df[self.type_col].tolist()
        ]].copy()
        return ok_biotypes

    def only_biotypes_with_a_large_enough_fraction_of_some_dataset(self, minimum_percentage=5):
        res = self.fraction_of_values_by_biotype()
        v = res.max(axis=1)
        ok = v[[x>(minimum_percentage) for x in v]].copy()
        ok_biotypes = set(ok.index)

        print(f"Keeping biotypes {ok_biotypes} and removing {set(v.index) - ok_biotypes}")
        self.df = self.df[[
            (x in ok_biotypes) for x in self.df[self.type_col].tolist()
        ]].copy()

        return ok_biotypes

    def fraction_of_values_by_biotype(self):

        for_df = []

        _df = self.df
        
        if len(_df.index) == 0:
            if type(self.value_col) == type([]):
                return pandas.DataFrame([{'Gene type': 'protein_coding', self.value_col[0]: 0}])
            return pandas.DataFrame([{'Gene type': 'protein_coding', self.value_col: 0}])
        
        biotypes = set(_df[self.type_col].tolist())
        total = _df[self.value_col].sum()
                
        res = 100*_df.groupby(by='Gene type')[self.value_col].sum()/total

        return pandas.DataFrame(res)

    def fraction_of_value_by_biotype_flat(self):
        res = self.fraction_of_values_by_biotype()

        d = res.to_dict('dict')
        rows = []
        for prot_col, dict_by_type in d.items():
            for _type, percent in dict_by_type.items():
                rows.append({
                    'Protein': prot_col, 'Gene type': _type,
                    '%': percent
                    })
        df = pandas.DataFrame(rows)
        return df

    def number_of_rnas_by_biotype(self):
        _df = self.df

        if len(_df.index) == 0:
            return None #pandas.DataFrame([{'Gene type': 'protein_coding', '# RNAs': 0}])

        total = len(_df.index)

        res = _df[self.type_col].value_counts()
        res = pandas.DataFrame(res)
        res['# RNAs'] = res[self.type_col]
        res['Gene type'] = res.index

        return res

    def biotypes_with_enough_reads(self, cutoff=50000):

        biotypes_to_keep = set()
        for biotype in set(self.df['Gene type'].tolist()):
            
            all_for_type = 0
            sub = self.df[self.df['Gene type']==biotype]
            
            for col in self.numeric_columns():
                all_for_type += sub[col].sum()
            if all_for_type >= cutoff:
                biotypes_to_keep.add(biotype)
        
        print("Biotypes with at least {} reads across all datasets (Ignoring 'Unknown'): {}.".format(
            cutoff, len(biotypes_to_keep - set(['Unknown']))))
        
        return biotypes_to_keep - set(['Unknown'])
