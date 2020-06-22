import pandas, importlib, os, collections
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#importlib.reload(utils)

class heatmapMaker:
    
    def __init__(self):
        pass

    @staticmethod
    def drop_odds(df):
        if '_ambiguous' in df.index:
            df.drop('_ambiguous', inplace=True)
        if '_no_feature' in df.index:
            df.drop('_no_feature', inplace=True)    

    def heatmap(self, _df, table_filename='spearman_correlations_table.xlsx',
                         fig_filename='heatmap_spearman_correlations.pdf',
               edit_column_names=True, cols=None,
               cutoff=5, **kwargs):

        df = _df.copy()
        self.drop_odds(df)

        if cols is None:

            n_cols = len(df.columns)

            if edit_column_names:
                replicate = collections.defaultdict(int)
                new_names = []
                for col in df.columns:
                    simple_col = col.split('_')[2]
                    replicate[simple_col] += 1
                    new_names.append(simple_col + ' R' + str(replicate[simple_col]))
                df.columns = new_names
                
        else:
            cols = [x for x in cols if x in df.columns]
            df = df[cols]
            n_cols = len(df.columns)
        
        print(df.head())
        df['sum'] = df.sum(axis=1, numeric_only=True)
        df['average'] = [x/max(n_cols, 1) for x in df['sum'].tolist()]

        df = df[df['average']>cutoff].copy()
        print("Clustering using {0} RNAs.".format(len(df.index)))

        del df['average']
        del df['sum']

        plt.clf()

        def _log(_):
            try:
                if _ <= 0:
                    return 0
                else:
                    return np.log10(_)
            except:
                    return 0

        #fig, ax = plt.subplots()
        m = df.copy()#df[hits_columns].copy()
        m = m.dropna(axis=0, how='any')

        corr_m = m.corr(method='spearman')

        if not os.path.exists(os.path.dirname(table_filename)):
            os.system('mkdir ' + os.path.dirname(table_filename))
            
        corr_m.to_excel(table_filename)
        
        #print("Correlation matrix for heatmap: {}".format(corr_m))
        
        for col in corr_m.columns:
            corr_m[col] = [np.nan_to_num(x) for x in corr_m[col].tolist()]

        fig = plt.figure()

        sns.set(font_scale=0.5)

        hm = sns.clustermap(corr_m, square=True, yticklabels=1, xticklabels=1,
                         #cmap=sns.light_palette('black', as_cmap=True),
                        #cmap=sns.cubehelix_palette(100, start=.5, rot=-2, as_cmap=True)
                        **kwargs)
                        
        #fig.tick_params(labelsize=1)

        #utils.draw(fig, '../figs/heatmap_FBL_hnRNPC.pdf')
        #plt.tight_layout()
        plt.savefig(fig_filename)
        plt.show()
        plt.clf()
        plt.close()
