import pandas, re
import numpy as np
from typing import List
import scipy
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt

# Has functions to calculate P values:
import xlSignificance

_pma_dir = '/Users/dp/pma/'

class xlLoader(xlSignificance.xlSignificance):
    """Load % cross-linking data for missense mutations.
    """

    def __init__(self, fname: str=f'{_pma_dir}/percentCrosslinked.xlsx'):
        self.fname = fname

    def load(self, fname: str='') -> List[pandas.DataFrame]:
        """Returns the whole dataframe and the missense proteins.
        returns [whole xl dataframe, recurrently mutated proteins dataframe]
        """
        if fname == '':
            fname = self.fname

        df = pandas.read_excel(
            self.fname,
            sheet_name='XL'
            )

        df = df.loc[[not pandas.isna(x) for x in df.Exp], :]
        df = df.loc[[pandas.isna(x) for x in df.Discard], :]

        df['Rep'] = [s.split(':')[-1] for s in df['Sample']]
        # Group_ID denotes every sample in the Exp sharing the same WT protein and
        # replicate number.
        df['Group_ID'] = [(exp, p, r) for exp,p,r in zip(
            df.Exp, df.Protein, df.Rep)]

        recurrent = df.loc[[(type(x) == type('') and 'Recurrent' in x) \
            for x in df['Category']], :]
        recurrent.loc[:,'Mut/WT'] = [{True: 'WT', False: 'Mut'}[bool('WT' in x)] \
            for x in recurrent.Category]

        recurrent = self.add_order_column(recurrent)  # Sorts.

        self.df = df
        self.recurrent = recurrent
        
        return df, recurrent

    def proteins_ordered_by_xl_rate(self, fname: str='') -> List[str]:

        if fname != '':
            (df, recurrent) = self.load(fname)

        xl = recurrent[recurrent['Label']=='% XL (minimal region)'].copy()
        self.as_fraction_of_wt(xl)

        self.proteins_ordered_by_xl_rate = list(dict.fromkeys(xl['Protein'].tolist()))

        return self.proteins_ordered_by_xl_rate

    @staticmethod
    def add_order_column(
        recurrent: pandas.DataFrame,
        values_label: str = r'% XL (minimal region)') -> pandas.DataFrame:
        
        xl = recurrent.loc[[x==values_label for x in recurrent['Label']], :]

        protein_order = xl.loc[xl['Mut/WT']=='WT', ['Protein', 'Value']]
        protein_order.index = protein_order.Protein
        protein_order.index.name = 'prt'
        _order = protein_order.groupby(by=['Protein'])['Value'].mean()
        
        protein_order = sorted(protein_order.index, key=lambda x: protein_order.loc[x, 'Value'].mean(), reverse=True)

        def get_order(name: str) -> int:
            if name in protein_order:
                return protein_order.index(name)
            elif name.split(' ')[0] in protein_order:
                return protein_order.index(name.split(' ')[0]) + 0.5
            return -1

        recurrent.loc[:,'order'] = [get_order(x) for x in recurrent.loc[:,'Protein']]
        recurrent.sort_values(by=['order'], inplace=True)
        print(f"Number of RBPs {len(set(protein_order))} in the order {protein_order}")
        return recurrent

    @staticmethod
    def exp_cat(exp):
        _ = exp.split('Exp')[-1]
        if int(_) < 70:
            return 'R3'
        elif exp == 'Exp81' or exp == 'Exp84':
            return 'R1'
        elif exp == 'Exp77' or exp == 'Exp80':
            return 'R2'
        else:
            return exp

    def get_wt_rep_mean(self, protein, group_id, df, values_col):
        #sub = df.loc[[protein==x for x in df.Protein], :]+
        group_id_of_wt = (group_id[0], group_id[1].split(' ')[0], group_id[2])
        sub = df.loc[[x==group_id_of_wt for x in df.Group_ID], :]
        # Get WT:
        sub = sub.loc[[protein==x for x in sub.Protein], :]
        if len(sub.index):
            
            #subrep = sub.loc[[rep==x for x in sub.Exp], :]
            subrep = sub.loc[[x==values_col for x in sub.Label], :]
            return subrep['Value'].mean()
        else:
            _group_id = (
                group_id[0], group_id[1].split(' ')[0], 
                re.sub('\d+', '1', group_id[2]))

            print(f"Trying {_group_id}")
            if str(_group_id) != str(group_id):
                return self.get_wt_rep_mean(protein, _group_id, df, values_col)
            else:
                print(f"No WT for {protein}, {group_id}")
                return None

    @staticmethod
    def get_n_replicates(protein, df, values_col):
        sub = df.loc[[protein==x for x in df.Protein], :]

        if len(sub.index):
            
            subrep = sub.loc[[x==values_col for x in sub.Label], :]
            n_reps = len(subrep.index)
            return n_reps
        return 0

    def get_matching_samples(df, group_id):
        # Sample: Exp:Protein:Replicate
        sub = df.loc[[x==group_id for x in df.Group_ID], :]

    def get_total_fmol(self, df, group_id):

        sub = df.loc[[x==group_id for x in df.Group_ID], :]
        pmol_protein = sub.loc[[x=='pmol protein' for x in sub.Label], 'Value']
        xl_lane = sub.loc[[x=="% XL (whole lane)" for x in sub.Label], 'Value']
        if len(pmol_protein) == 0 or len(xl_lane) == 0:
        #    print(f"Problem: {group_id}, sub={sub}")

            return np.nan
        pmol_RNA_whole_lane = 0.01 * float(pmol_protein) * float(xl_lane)
        fmol_RNA_whole_lane = 1000 * pmol_RNA_whole_lane

        return fmol_RNA_whole_lane

    def add_total_fmol_column(self, df=None, inplace=True):
        if df is None:
            df = self.df

        lab = 'fmol RNA (whole lane)'

        d = df.to_dict('records')
        new_rows = {}
        for row in d:
            _r = dict(zip(row.keys(), row.values()))
            _r['Value'] = self.get_total_fmol(df, row['Group_ID'])
            _r['Label'] = lab
            new_rows[row['Group_ID']] = _r

        for row in new_rows.values():
            d.append(row)

        _df = pandas.DataFrame(d, columns=df.columns)
        _df.sort_index(by=['Exp', 'Protein', 'Label'], inplace=True)

        if inplace:
            self.df = _df
        return _df

    def as_fraction_of_wt(
        self,
        df: pandas.DataFrame,
        values_label: str = '% XL (minimal region)'):
    
        df.loc[:,'WT protein'] = [x.split(' ')[0] for x in df.loc[:, 'Protein']]

        df.loc[:,'WT rep mean'] = [
            self.get_wt_rep_mean(protein, group_id, df, values_label) \
            for (protein, group_id, sample) in zip(
                df['WT protein'], df['Group_ID'], df['Sample'])]
        print('---')
        df.loc[:,'Value vs WT'] = [x/max([wt, 10**-6]) \
            for x,wt in zip(df['Value'], df['WT rep mean'])]
        df.loc[:,'log2 Value vs WT'] = np.log2(df.loc[:, 'Value vs WT'].tolist())
        df.loc[:,'# XL replicates'] = [
            self.get_n_replicates(protein, df, values_label) \
            for protein in df.loc[:,'WT protein']]

    def effect_sizes(self, xl):

        def effect_size(arr):
            theta = (np.mean(arr) - 1) # arr = [MUT(batch 1)/WT(batch 1), MUT(batch 2)/WT(batch 2), ...]
            # If null hypothesis (MUT=WT) is true, np.mean(arr) ~ 1.
            theta = theta/np.std(arr)
            return theta

        # Get arrays of [MUT(batch 1)/WT(batch 1), MUT(batch 2)/WT(batch 2), ...]
        value_vs_wt_by_prot = xl.groupby('Protein')['Value vs WT'].apply(np.array)

        # Put effect sizes: (mean(MUT/WT) - 1)/std_deviation(MUT/WT),
        # where there are n MUT/WT values for n batches.
        mut_to_effect_size = {}
        for name, val in zip(value_vs_wt_by_prot.index, value_vs_wt_by_prot):
            mut_to_effect_size[name] = effect_size(val)
            print(f"{name}: {val} effect_size= {mut_to_effect_size[name]}")

        xl['Effect size (XL)'] = [mut_to_effect_size.get(name, 0) for name in xl.Protein]

        return mut_to_effect_size
