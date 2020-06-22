
import glob, re

def get_protein_lengths_from_tcga_lolipop_svg_files():
    protein_lengths = {}
    input_fnames = [fname for fname in glob.glob('figs/edited_lolipop/*svg')]

    for n, fname in enumerate(input_fnames):

        lines = open(fname, 'r').readlines()
        gene = re.search('lolipop/(.+)_lollipop.svg', fname).group(1)

        # Example:
        #<text text-anchor="middle" style="font-size: 10px; font-family: arial;" x="922.5999999999999" y="177" \
        # dy="1em">1922aa</text>

        for n, li in enumerate(lines):
            if '<text' not in li:
                continue
            val = re.search('>([^<]+)aa</text', li)
            if not val:
                continue
            protein_length = int(val.group(1))
            print(gene, protein_length)
            protein_lengths[gene] = protein_length
            
    #protein_lengths['PABPC4L'] = 370

    return protein_lengths