import pandas, re, collections, os
import sameRiver
import sameRiver.biotypeLookupFileMaker

def mask_low_raw_read_counts(nb, scheme, which, raw_read_cutoff=10):
    above_cutoff = raw_read_counts_above_cutoff(nb, scheme, cutoff=raw_read_cutoff)
    for prot in nb.positives.numeric_columns(nb.pvals_dfs[which]):
        print(prot)

        if prot not in nb.pvals_dfs[which]:
            continue
        nb.pvals_dfs[which][prot] = [{True: p, False: 1}[bool(rna in above_cutoff[prot])] for p, rna in
                             zip(nb.pvals_dfs[which][prot], nb.pvals_dfs[which].index)]
        print('---', prot, len(nb.targets(prot, which=which)))
    return above_cutoff


def raw_read_counts_above_cutoff(nb, scheme, cutoff=10, protein_list=None):
    
    df = nb.positives.raw_counts_df
    numeric_columns = nb.positives.numeric_columns(df)
    
    if protein_list is None:
        protein_list = list(set([scheme.gene_from_fname(x) for x in numeric_columns]))
        
    above_cutoff = {}
    for protein in protein_list:
        a = df[[x for x in numeric_columns if (scheme.gene_from_fname(x) == protein)]].copy()
        a['max'] = a.max(axis=1)
        a = a[a['max']>cutoff].copy()
        a.sort_values(by='max', inplace=True)
        print(protein, a.shape, end=', ')
        above_cutoff[protein] = list(a.index)
    return above_cutoff


def biotype_filter(df):

    _types = set(df['Gene type'])

    df.loc[:,'Gene type'] = [{True:'7SK', False:x}[bool(x=='RNA')] for x in df['Gene type']]
    _types = set(df['Gene type'])

    df.loc[:,'Gene type'] = [re.sub('.*/.*', 'Repetitive element', x) for x in df['Gene type']]
    #df['Gene type'] = [re.sub('LTR/.*', 'LTR', x) for x in df['Gene type']]
    #df['Gene type'] = [re.sub('DNA/hAT.*', 'DNA/hAT', x) for x in df['Gene type']]
    verboten_types = [
        'Repetitive element',
        'rRNA',
        'retained_intron',
        'Unknown',
        '7SK',
        'nonsense_mediated_decay'
    ]
    #df = df.loc[[(x not in verboten_types) for x in df['Gene type']], :].copy()

    #df = df[df['Gene type']!='processed_transcript']
    #df = df[df['Gene type']!='nonsense_mediated_decay']
    #df = df[[x for x in df.columns if (scheme.gene_from_fname(x) not in negative_metadata.random_proteins + [
    #    'EPB41L5', 'CCIN', 'HCT116', 'PCBP1:100P', 'PCBP1:100Q', 'PCBP1:dKH'])]]
    return df


def make_enst_gene_name_biotype_map_file(
    gtf_filename='/opt/genome/Homo_sapiens.GRCh38.83/exon_CDS_ensembl_havana_only_Homo_sapiens.GRCh38.83',
    out_filename='enst_transcript_id_name_biotype_map.txt',
                         ):
    biotypeLookupFileMaker.biotypeLookupFileMaker(
        gtf_filename=gtf_filename, out_filename=out_filename)

def make_biotype_map_file(
    to_biotype_file='/opt/genome/transcripts.gtf',
    refseq_to_ensg='/opt/genome/refseq_to_ensg.txt', 
    outfname='/opt/genome/gene_name_to_biotype.txt',
    enst_name='/opt/genome/enst_name.txt'):

    biotypeLookupFileMaker.biotypeLookupFileMaker(
        to_biotype_file=to_biotype_file, refseq_to_ensg=refseq_to_ensg,
        outfname=outfname, enst_name=enst_name)

def make_simple_mapping_file():
    biotypeLookupFileMaker.biotypeLookupFileMaker.make_simple_mapping_file()

not_used = """
def add_gene_name_and_type_from_refseq(df, refseq_name='/opt/genome/refseq/biotype_refseq_geneid_symbol.txt',
                     add_name=True, add_type=True):
    
    transl = pandas.read_csv(refseq_name, sep='\t')
    
    df['Gene'] = [str(x).split('.')[0] for x in df.Gene]
    
    if add_name:
        to_symbol = dict(zip(
            [x.split('.')[0] for x in transl['Refseq'].tolist()],
            transl['Symbol'].tolist()
        ))
        df['Gene name'] = [to_symbol.get(x, 'Unknown') for x in df['Gene'].tolist()]
        
    if add_type:
        to_type = dict(zip(
            [x.split('.')[0] for x in transl['Refseq'].tolist()],
            transl['Gene type'].tolist()
        ))    
        df['Gene type'] = [to_type.get(x, 'Unknown') for x in df['Gene'].tolist()]
        
    return df"""

