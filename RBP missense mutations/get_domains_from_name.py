import pandas, os, sys, re, time, collections, Bio, pprint, random, glob

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

id_mappings = pandas.read_csv('./old_pma_analysis_inputs/id_mappings.txt', sep='\t')
domains = pandas.read_csv('./old_pma_analysis_inputs/ensg_to_protein_domain.txt', sep='\t')

ensg_to_domains = collections.defaultdict(set)
for ensg, desc in zip(domains['Ensembl Gene ID'].tolist(), domains['Interpro Short Description'].tolist()):
    ensg_to_domains[ensg].add(desc)

symbol_to_domains = collections.defaultdict(set)
for symbol, ensg in zip(id_mappings['Gene Symbol'], id_mappings['Ensembl Gene ID']):
    symbol_to_domains[symbol] |= ensg_to_domains.get(ensg, set())

def get_domains_from_name(name: str) -> str:
	return symbol_to_domains.get(name, set())
	#if name in ensg_to_domains:
	#	return ensg_to_domains[name]