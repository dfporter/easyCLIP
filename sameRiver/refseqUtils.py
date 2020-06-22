
# coding: utf-8



def simplify_gene2refseq(refseq_name='/opt/genome/refseq/gene2refseq'):
    """Downloaded from ftp://ftp.ncbi.nih.gov/gene/DATA/"""
    transl = {}
    with open(refseq_name) as f:
        for li in f:
            s = li.split('\t')
            transl[s[3]] = s[1]
    with open('/opt/genome/refseq/simple_gene2refseq', 'w') as f:
        for k, v in transl.items():
            f.write('\t'.join([k, v]) + '\n')

            
def make_gene2refseq2symbol_human(
        gene2refseq='/opt/genome/refseq/simple_gene2refseq',
        human_info='/opt/genome/refseq/Homo_sapiens.gene_info'):
    """From:
    ftp://ftp.ncbi.nih.gov/gene/DATA/
    ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/"""
    transl = {}
    with open(gene2refseq) as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            transl[s[0].split('.')[0]] = s[1]
    id_to_symbol = {}
    with open(human_info) as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            #geneid = s[1]
            #symbol = s[2]
            id_to_symbol[s[1]] = s[2]
    outli = "Refseq\tGeneID\tSymbol\n"
    for refseq in transl:
        if transl[refseq] in id_to_symbol:
            outli += '\t'.join([refseq, transl[refseq], id_to_symbol[transl[refseq]]]) + '\n'
        #else:
            #outli += '\t'.join([refseq, transl[refseq], 'Unknown']) + '\n'
    with open('/opt/genome/refseq/refseq_geneid_symbol.txt', 'w') as f:
        f.write(outli)

def make_refseq_to_biotype(
        refseq_symbol='/opt/genome/refseq/refseq_geneid_symbol.txt',
        name_biotype='/opt/genome/gene_name_to_biotype.txt'
                ):
    df = pandas.read_csv(name_biotype, sep='\t')
    to_biotype = dict(zip(df['Gene name'].tolist(),
                         [x.strip('\n') for x in df.Biotype]))
    with open(refseq_symbol) as f:
        outli = next(f)
        outli = outli.rstrip('\n') + '\tGene type\n'
        for li in f:
            #print(li)
            s = li.rstrip('\n').split('\t')
            biotype = to_biotype.get(s[2], 'Unknown')
            outli += '\t'.join(s + [biotype + '\n'])
    with open('/opt/genome/refseq/biotype_refseq_geneid_symbol.txt', 'w') as f:
        f.write(outli)

