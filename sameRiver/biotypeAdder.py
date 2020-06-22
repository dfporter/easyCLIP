import dill, pandas, os

class biotypeAdder():

    def save(self, fname='data/biotypeAdder.dill'):
        if not os.path.exists('./data/'):
            os.system('mkdir ./data/')
        
        print("Saving to {}...".format(fname))
        with open(fname, 'wb') as f:
            dill.dump(self, f)

    @staticmethod
    def load(fname='data/biotypeAdder.dill'):
        with open(fname, 'rb') as f:
            counts = dill.load(f)
        return counts

    def create(self, map_fname='enst_transcript_id_name_biotype_map.txt',):
        print("Reading {} to create biotypeAdder...".format(map_fname))
        to_biotypes = pandas.read_csv('enst_transcript_id_name_biotype_map.txt', sep='\t')
        self.to_biotypes = dict(zip(to_biotypes['gene_name'].tolist(),
                              to_biotypes['transcript_biotype'].tolist()))

    def add_biotypes_column_from_gene_name(
        self, df, map_fname='enst_transcript_id_name_biotype_map.txt',):

        if not hasattr(self, 'to_biotypes'):
            self.create(map_fname=map_fname)

        def repeat(_str):
            if _str[0] == '>':
                return _str.split('::')[0].split('#')[-1]
            else:
                return _str
        
        print("Adding biotypes...")
        for col in ['Gene name', 'gene_name']:
            if col in df.columns:
                df.loc[:, 'Gene type'] = [self.to_biotypes.get(x.split('::')[0], x) for x in df[col]]
                df.loc[:, 'Gene type'] = [repeat(x) for x in df['Gene type']]
                return df
        
        print("Did not find gene_name/Gene name column, so using index:", df.index)
        df.loc[:,'Gene type'] = [self.to_biotypes.get(x.split('::')[0], x) for x in df.index]
        df.loc[:,'Gene type'] = [repeat(x) for x in df['Gene type']]

        return df