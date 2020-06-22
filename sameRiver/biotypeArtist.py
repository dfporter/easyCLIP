import pandas, importlib
import sameRiver
import numpy as np
import sameRiver.biotypeUtils
import sameRiver.intronExonSplit
import sameRiver.stacked_bargraph
import sameRiver.biotypeDf
importlib.reload(sameRiver.biotypeUtils)
importlib.reload(sameRiver.stacked_bargraph)


def combine_protein_replicates(df, rep_col, type_col, val_col, scheme):
    df['Protein'] = [
        scheme.gene_from_fname(x) for x in df[rep_col]]
    df.groupby(['Protein', type_col])[val_col].mean().reset_index()


def dataframe_of_target_biotypes(
    df, pvals, scheme, min_read_cutoff=None, which='per_read', p_val_to_define_targets=0.01):

    fraction_reads_in_targs = pandas.DataFrame(columns=['Protein', 'Gene type', '% reads (target RNAs)'])
    fraction_targets = pandas.DataFrame(columns=['Protein', 'Gene type', '# RNAs',])
    fraction_all_reads = pandas.DataFrame(columns=['Protein', 'Gene type', '% reads (all RNAs)',])
    
    df_raws = []
    df = df.copy()
    
    #df = nb.positives.counts_per_million_df.copy()
    #df = ba.add_biotypes_column_from_gene_name(df)

    get_numeric_columns = lambda _df: [x for x in _df.columns if (_df[x].dtype.kind in 'bifc')]

    numeric_columns = get_numeric_columns(df)
    from_fname = scheme.gene_from_fname

    bio = sameRiver.biotypeDf.biotypeDf(df=df[numeric_columns + ['Gene type']].copy())
    ok_biotypes = bio.only_biotypes_with_a_large_enough_fraction_of_some_dataset(
        minimum_percentage=5)
    df = df[[(x in ok_biotypes) for x in df['Gene type']]].copy()

    print(f'pizza: now at {df["Gene type"].value_counts()}')
    done = set()
    for protein in set([from_fname(col) for col in numeric_columns]):
        #print(f'pizza: {protein} in biotypeArtist:')
        if protein not in pvals.columns:
            print(f"No column for {protein} found in pvalues dataframe. Skipping.")
            continue

        cols = [col for col in numeric_columns if (from_fname(col)==protein)]

        bio = sameRiver.biotypeDf.biotypeDf(df=df[cols + ['Gene type']].copy())

        frac_values_all_rnas = bio.fraction_of_value_by_biotype_flat()
        combine_protein_replicates(frac_values_all_rnas, 'Protein', 'Gene type', '%', scheme)
        frac_values_all_rnas['% reads (all RNAs)'] = frac_values_all_rnas['%']
        #del frac_values_all_rnas['%']
        
        # Subset our dataframe to targets.
        #print(f'pizza: pvals[protein]={pvals[protein]}')
        targets = set(pvals.loc[[x<p_val_to_define_targets for x in pvals[protein]], :].index)
        #print(f'pizza: targets = {targets}')
        #print(f"pizza: df.head() = {df.head()}")
        target_df = df.loc[[(x in targets) for x in df.index],:].copy()
        #print(f"pizza: target_df.head() = {target_df.head()}")
        #print(f"pizza: target_df['Gene type'].value_counts()={target_df['Gene type'].value_counts()}")
        bio = sameRiver.biotypeDf.biotypeDf(df=target_df[cols + ['Gene type']].copy())
        frac_values_target_rnas = bio.fraction_of_value_by_biotype_flat()
        combine_protein_replicates(frac_values_all_rnas, 'Protein', 'Gene type', '%', scheme)
        frac_values_target_rnas['% reads (target RNAs)'] = frac_values_target_rnas['%']

        frac_targs = bio.number_of_rnas_by_biotype()
        if frac_targs is None:
            continue
        frac_targs['Protein'] = protein

        frac_targs = frac_targs[['Protein', 'Gene type', '# RNAs',]]
        #print(f"pizza: frac_targs={frac_targs}")
        del frac_values_target_rnas['%']
        del frac_values_all_rnas['%']
        
        if len(fraction_reads_in_targs.index) == 0:
            fraction_reads_in_targs = frac_values_target_rnas.copy()
        else:
            fraction_reads_in_targs = fraction_reads_in_targs.append(frac_values_target_rnas)
        if len(fraction_targets.index) == 0:
            fraction_targets = frac_targs.copy()
        else:
            fraction_targets = fraction_targets.append(frac_targs)#.extend(frac_targs.to_dict('records'))
        if len(fraction_all_reads.index) == 0:
            fraction_all_reads = frac_values_all_rnas.copy()
        else:
            fraction_all_reads = fraction_all_reads.append(frac_values_all_rnas)
    
    for df in [fraction_reads_in_targs, fraction_targets, fraction_all_reads]:
        df.index = df['Protein'].tolist()
    n = """
    #print(f'pizza: fraction_all_reads={fraction_all_reads}')
    fraction_reads_in_targs = pandas.DataFrame.from_dict(fraction_reads_in_targs)
    fraction_reads_in_targs.replace(np.nan, 0, inplace=True)
    fraction_targets = pandas.DataFrame.from_dict(fraction_targets)
    fraction_targets.replace(np.nan, 0, inplace=True)
    fraction_all_reads = pandas.DataFrame.from_dict(fraction_all_reads)
    fraction_all_reads.replace(np.nan, 0, inplace=True)
    fraction_reads_in_targs.columns = ['Protein', 'Gene type', '% reads (target RNAs)']#, 'Fraction']
    fraction_targets.columns = ['Protein', 'Gene type', '# RNAs',]
    fraction_all_reads.columns = ['Protein', 'Gene type', '% reads (all RNAs)',]#, 'Fraction']"""
    #print("pizza biotypeArtist(): ", "\n*\n".join(
    #    map(str, [fraction_reads_in_targs, fraction_targets, fraction_all_reads])))
    return [fraction_reads_in_targs, fraction_targets, fraction_all_reads]