
import os, re, matplotlib, pandas, functools, operator, collections
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def _plot(
        df, hue_order=[], hue_label='', sample_label='', value_label='',
        sample_order=[], fname='biotypes_bargraph.pdf',
        title='', plot_kwargs={}, group_to_color=None, fig_height=4):

    # Value initiations.
    plt.clf()
    fig = plt.figure()

    sns.set_style('ticks')

    if group_to_color is None:
        group_to_color = dict(zip(
            hue_order, sns.cubehelix_palette(len(hue_order), start=.5, rot=-2)))
    
    width = 0.75
    to_val = {}

    for sample in sample_order:  # Typically these are the proteins.
        sub = df[df[sample_label]==sample].copy()
        to_val[sample] = dict(zip(sub[hue_label].tolist(), sub[value_label].tolist()))

    ind = np.arange(0.2, len(sample_order), (len(sample_order)-0.2)/len(sample_order))
    
    bars = {}
    max_value_with_hue_across_all_samples = collections.defaultdict(float)
    print(f"pizza: hue_order before sorting={hue_order}")
    print(f"pizza: group_to_color={group_to_color}")
    for y_pos, sample in enumerate(sample_order):  # For each protein.
        hue_order = sorted(hue_order, key=lambda x: to_val[sample].get(x, 0))[::-1]

        for n, hue in enumerate(hue_order):
            in_group = df[df[hue_label]==hue].copy()
            #vals = [float(x) for x in in_group[hue_order].tolist()]
            if hue not in to_val[sample]:
                print(f"pizza: {hue} not in to_val[{sample}]. to_val[sample]={to_val[sample]}")
            vals = to_val[sample].get(hue, 0)
            bars[n] = plt.bar(
                ind[y_pos], 
                vals, 
                width, color=group_to_color[hue],# alpha=0.5,
                linewidth=0,
                #orient='h',
                label=hue)
            max_value_with_hue_across_all_samples[hue] = max([vals,
                max_value_with_hue_across_all_samples[hue]])
#            print('vals, ind: {0}, {1}'.format(vals, ind[y_pos]))
#            print('label {0}'.format(ylabel_order))
    patches = []
    plt.xticks(ind, sample_order)
    for name, color in group_to_color.items():
        if max_value_with_hue_across_all_samples[name] > 0:
            patches.append(matplotlib.patches.Patch(
                color=color, linewidth=0, label=name))
            plt.legend(handles=patches,
                bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                ncol=2, mode="expand", borderaxespad=0.
            )
    #plt.figure(figsize=(1,1))
    plt.ylabel(str(value_label))
    plt.xlabel(str(sample_label))
    plt.xticks(rotation=90)

#    topbar = plt.Rectangle((0,0),1,1,fc="red", edgecolor = 'none')
#    bottombar = plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor = 'none')
#    sns.factorplot(x='Reads (% total)', y='P6', hue='Gene type', kind='bar', data=df)

    if len(sample_order) > 30:
        fig.set_figwidth(6)
    else:
        fig.set_figwidth(0.3 * len(sample_order))  # 0.3 for publication
    fig.set_figheight(fig_height)
    
    if title:
        plt.title(title)
    
    plt.gcf().savefig(fname)
    plt.show()


def make_xvalues_cumulative(value_label, hue_label, hue_order, _df):
    as_cumulative = [0]
    cumulative_dict = {}
    group_val = dict(zip(_df[hue_label].tolist(), _df[value_label].tolist()))
    for group in hue_order:
        as_cumulative.append(as_cumulative[-1] + group_val[group])
        cumulative_dict[group] = as_cumulative[-1]

    _df[value_label] = [cumulative_dict[x] for x in _df[hue_label].tolist()]
    #print(f'pizza: make_xvalues_cumulative ({value_label}, {hue_label}, {hue_order}),')
    #print(f'\t_df=\n{_df}')
    return _df


def stacked_bargraph(
        df=pandas.DataFrame(), sample_label='', value_label='', hue_label='',
        sample_order=[], title='',
        plot_kwargs={}, fname='biotypes_bargraph.pdf', group_to_color=None, **kwargs):

    hues = list(set(df[hue_label].tolist()))  # Biotypes, for example.
    sample_set = set(df[sample_label].tolist())  # Proteins, for example.
    
    usable_proteins_ordered = [x for x in sample_order if x in sample_set]

    print('Groups {0}\n proteins (ordered) {1}'.format(hues, usable_proteins_ordered))
    
    hue_order = hues
    print(f"hue_label={hue_label}. hue_order={hues}.")

    df = df[[sample_label, hue_label, value_label]].copy()

    df_of_y = []
    
    hue_values_observed = set()
    for sample in usable_proteins_ordered:
        
        sub = df[df[sample_label]==sample].copy()
        to_val = dict(zip(sub[hue_label].tolist(), sub[value_label].tolist()))

        # Fill in missing values.
        to_add = []
        for hue in [x for x in hues if x not in to_val]:
            to_val[hue] = 0
            sub = sub.append({
                hue_label: hue,
                sample_label: sample,
                value_label: 0
            }, ignore_index=True)
        
        # Add missing values.
        #_df = pandas.DataFrame(to_add)
        #print(f'pizza: sub={sub}')#' to_add={to_add}')
        #for _dict in to_add:
        #    sub = sub.append(_dict, ignore_index=True)
        hue_order = sorted(hue_order, key=lambda x: to_val[x])
        
        # Using the group order, make the values for each group cumulative.
        df_of_y.append(make_xvalues_cumulative(
            value_label, hue_label, hue_order, sub.copy()))
        #print(f'df_of_y[-1]={df_of_y[-1]}')

    _df = pandas.concat(df_of_y)
    #print(f'concat: {_df}')
    """
    # Put the list of lists of dicts [[{}, {}..], [{},...]] together into one table.
    if len(df_of_y) > 1:
        _df = pandas.DataFrame(functools.reduce(operator.add, df_of_y, []))
    elif len(df_of_y) == 1:
        _df = pandas.DataFrame(df_of_y[0])
    else:
        print("---\nEmpty data frame to plot:", df_of_y)
        return
    """
    # Label columns correctly. Yes, this is terrible.
    new_names = []
    for col in _df.columns:
        x = _df[col].tolist()[0]
        if type(x) != type(''):
            new_names.append(value_label)
        elif x in sample_set:
            new_names.append(sample_label)
        else:
            new_names.append(hue_label)

    _df.columns = new_names#[sample_label, hue_label, value_label]
    #print(f"pizza: _df[sample_label, hue_label, value_label] = {_df.head()}")
    _plot(
        _df, hue_order=hue_order, hue_label=hue_label,
        sample_label=sample_label, value_label=value_label,
        sample_order=usable_proteins_ordered, fname=fname,
        title=title, group_to_color=group_to_color,
        plot_kwargs=plot_kwargs, **kwargs)
    
    
