import os, re, matplotlib, pandas, collections, importlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
import pandas, re
import seaborn as sns
import scipy as sp
from typing import List, Union, Mapping, Tuple


def subset_by_label_and_convert_to_matrix(
    input_df: pandas.DataFrame, label: str='% XL (minimal region)',
    xcol: str='Protein', ycol: str='Value', order_by: Union[List, None]=None,
    verbose=False,
    ) -> pandas.DataFrame:
    
    subset_by_label_df = input_df.loc[input_df['Label']==label, :].copy()
    subset_by_label_df[label] = subset_by_label_df[ycol]
    
    proteins_as_means_df = convert_to_matrix(subset_by_label_df, xcol=xcol, ycol=label)
    
    if verbose:
        print(f"label: {label}, xcol: {xcol}, ycol: {ycol}. proteins_as_means_df:")
        print(proteins_as_means_df.head())

    if order_by is not None:
        proteins_as_means_df = proteins_as_means_df[order_by]
        
    return proteins_as_means_df

def convert_to_matrix(
    input_df: pandas.DataFrame, xcol: str='Protein', ycol: str='Value') -> pandas.DataFrame:

    g = input_df.groupby(['Protein'])[ycol].mean()

    return pandas.DataFrame(g).T

def make_categorical_column_numeric(
    series: List[str], return_reverse_mapping=False) -> List[int]:

    text_to_num = {}
    num_to_text = {}
    numeric_col = []
    
    #if 'No' in list(series):
    #    text_to_num['No'] = 0
    series = list(series)
    for text in series:
        if text not in text_to_num:
            if 1 not in text_to_num.values():
                text_to_num[text] = 1
                num_to_text[1] = text
                lowest_int, text_with_lowest = (1, text)
            else:
                n = len(text_to_num) + 1
                text_to_num[text] = n
                num_to_text[n] = text
                highest_int, text_with_highest = (n, text)
    
    #if 'Yes' in list(series):
    #    former_yes_val = text_to_num[text]
    #    max_value = max(text_to_num.values())
    #    key_at_max_value = [x for x in text_to_num if text_to_num[x]==max_value][0]
    #    text_to_num['Yes'] = max_value
    #    text_to_num[key_at_max_value] = former_yes_val
    
    if 'No' in series and text_with_lowest != 'No':
        hold_num = text_to_num['No']
        text_to_num['No'] = lowest_int
        text_to_num[text_with_lowest] = hold_num
        num_to_text[lowest_int] = 'No'
        num_to_text[hold_num] = text_with_lowest

    if 'Yes' in series and text_with_highest != 'Yes':
        hold_num = text_to_num['Yes']
        text_to_num['Yes'] = highest_int
        text_to_num[text_with_highest] = hold_num
        num_to_text[highest_int] = 'Yes'
        num_to_text[hold_num] = text_with_highest

    #print(series, text_to_num)
    for text in series:
        numeric_col.append(text_to_num[text])
    
    #print(text_to_num)
    if return_reverse_mapping:
        return numeric_col, dict(zip(text_to_num.values(), text_to_num.keys()))
    else:
        return numeric_col

def nan_to_no(series: List) -> List:
    out_list = [x := (item if not pandas.isna(item) else 'No') for item in series]
    return out_list
    for item in series:
        if pandas.isna(item):
            out_list.append('No')
        else:
            out_list.append(item)
    return out_list