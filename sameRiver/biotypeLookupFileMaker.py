import pandas, re, collections, os

# Not used:
ignored_biotypes = set([
'antisense', 'nonsense_mediated_decay', 
'processed_transcript',
'retained_intron',
'transcribed_processed_pseudogene', 'unprocessed_pseudogene', 
'processed_pseudogene', 'sense_overlapping',
'transcribed_unprocessed_pseudogene', 
'translated_unprocessed_pseudogene',
'macro_lncRNA',
'TEC'
'IG_C_gene',
'IG_C_pseudogene',                      
'IG_D_gene',              
'IG_J_gene',                             
'IG_V_gene',                             
'IG_V_pseudogene',                       
'TR_C_gene',                             
'TR_J_gene',                             
'TR_J_pseudogene',
'TR_V_gene', 
'TR_V_pseudogene',
'IG_C_gene',
'unitary_pseudogene', 
'bidirectional_promoter_lncrna',
 'antisense_RNA', 
 'polymorphic_pseudogene',  'TEC', 
 'transcribed_unitary_pseudogene',
  '3prime_overlapping_ncrna', 'sense_intronic',
])

snRNA = ["U6", "U2", "U13_", "U7", "U3", "U1", "U5", "U13", "U4", "U8", "U17", "U14",]
rRNA = ['LSU-rRNA_Hsa', 'SSU-rRNA_Hsa', '5S']

def make_enst_gene_name_biotype_map_file(
    gtf_filename='/opt/genome/Homo_sapiens.GRCh38.83/exon_CDS_ensembl_havana_only_Homo_sapiens.GRCh38.83',
    out_filename='enst_transcript_id_name_biotype_map.txt'):

    rows = collections.defaultdict(dict)
    
    def re_result(search_result):
        if search_result is None:
            return ''
        else:
            return search_result.group(1)
        
    with open(gtf_filename) as f:
        for li in f:
            #s = li.rstrip('\n').split('\t')
            txpt = re_result(re.search('transcript_id=* *"*([^;"]+)', li))     
            name = re_result(re.search('gene_name=* *"*([^;"]+)', li))
            biotype = re_result(re.search('transcript_biotype=* *"*([^;"]+)', li))

            if (biotype == ''):
                biotype = re_result(re.search('biotype=* *"*([^;"]+)', li))
            if biotype == '':
                biotype = re_result(re.search('transcript_type=* *"*([^;"]+)', li))
            if biotype == '':
                biotype = re_result(re.search('gene_type=* *"*([^;"]+)', li))
            if biotype == '':
                biotype = re_result(re.search('type=* *"*([^;"]+)', li))

            # Transcripts from the repeats genome.
            if biotype == 'repeat' and 'tRNA-' in name:
                biotype = 'tRNA'
            if name in snRNA:
                biotype = 'snRNA'
            if name in rRNA:
                biotype = 'rRNA'

            if txpt == '':
                continue
            
            if txpt not in rows:
                rows[txpt] = {'transcript_id': txpt, 'gene_name': name, 'transcript_biotype': biotype}
    
    os.makedirs(os.path.dirname(out_filename), exist_ok=True)
    with open(out_filename, 'w') as f:
        keys = ['transcript_id', 'gene_name', 'transcript_biotype']
        f.write('\t'.join(keys) + '\n')
        for row in rows.values():
            f.write('\t'.join([row[x] for x in keys]) + '\n')


def make_biotype_map_file(to_biotype_file='/opt/genome/transcripts.gtf',
                          refseq_to_ensg='/opt/genome/refseq_to_ensg.txt', 
                          outfname='/opt/genome/gene_name_to_biotype.txt',
                          enst_name='/opt/genome/enst_name.txt',
                         ):
    
    to_ensg = pandas.read_csv(enst_name, sep='\t')
    to_ensg_dict = dict(zip(to_ensg['gene_name'].tolist(),
                 to_ensg['transcript_id'].tolist()))
    if 'SNORD3C' in to_ensg_dict:
        print(to_ensg_dict['SNORD3C'])
    print(str(to_ensg_dict)[:100])
    to_biotype = {}
    with open(to_biotype_file) as f:
        for li in f:
            if li[0] == '#':
                continue

            biotype = re.search('gene_biotype "([^"]+)"', li)
            txpt_id = re.search('transcript_id "([^"]+)"', li)
            if (biotype is not None) and (txpt_id is not None):
                to_biotype[txpt_id.group(1)] = biotype.group(1)
    #print(str(to_biotype)[-1100:])
    
    outli = "Gene name\tBiotype\n"
    for name, ensg in to_ensg_dict.items():
        outli += '{0}\t{1}\n'.format(name, to_biotype.get(ensg, 'Unkown'))
    #print(outli[:1000])
    with open(outfname, 'w') as f:
        f.write(outli)

def make_simple_mapping_file():
    mapping_fname='/opt/genome/Homo_sapiens.GRCh38.83/enst_transcript_id_name_biotype_map.txt'
    to_biotypes_df = pandas.read_csv(mapping_fname, sep='\t')

    order = ['vaultRNA','snRNA', 'snoRNA', 'scaRNA', 
     'protein_coding',
     'lincRNA',
     'non_coding']
    other_types = [
     'antisense_RNA', 
     'transcribed_unitary_pseudogene', 
     'IG_J_pseudogene', 'TR_D_gene', 'TR_J_pseudogene', 
     'bidirectional_promoter_lncrna',
     'TR_C_gene', 'IG_C_pseudogene', 'IG_C_gene', 'IG_J_gene', 'TR_V_pseudogene', 
     '3prime_overlapping_ncrna', 'IG_D_gene', 
     'TR_J_gene', 'TR_V_gene', 'IG_V_gene', 
     'IG_V_pseudogene', 'unitary_pseudogene', 
     'polymorphic_pseudogene', 
     'non_stop_decay', 
     'transcribed_unprocessed_pseudogene', 
     'sense_overlapping',  'sense_intronic', 
     'transcribed_processed_pseudogene',
     'TEC', 
     'unprocessed_pseudogene', 'processed_pseudogene', 
     'antisense',  'nonsense_mediated_decay', 
     'retained_intron', 'processed_transcript', 
     'translated_unprocessed_pseudogene', 'macro_lncRNA', 
    ]

    to_biotypes = collections.defaultdict(set)
    all_types = collections.defaultdict(int)

    for gene_name, biotype in zip(to_biotypes_df.gene_name, to_biotypes_df.transcript_biotype):
        to_biotypes[gene_name].add(biotype)
        all_types[biotype] += 1

    keys = sorted(all_types, key=lambda x: all_types[x])

    to_a_biotype = {}
    for gene_name, biotypes in to_biotypes.items():

        if len(biotypes) == 1:
            to_a_biotype[gene_name] = list(biotypes)[0]
            continue

        _types = [t for t in list(biotypes) if t in order]

        if len(_types):
            _i = sorted([order.index(t) for t in _types])
            to_a_biotype[gene_name] = order[_i[0]]
            continue

        _types = [t for t in list(biotypes) if t in other_types]
        if len(_types):
            _i = sorted([other_types.index(t) for t in _types])
            to_a_biotype[gene_name] = other_types[_i[0]]
            continue

        to_a_biotype[gene_name] = list(biotypes)[0]

    df = pandas.DataFrame.from_dict(orient='index', data=to_a_biotype)
    df['gene_name'] = df.index
    df['Biotype'] = df[0]
    del df[0]
    mapper = os.path.dirname(mapping_fname) + '/gene_name_to_biotype.txt'
    df.to_csv(mapper, sep='\t')#, index=False)