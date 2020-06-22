import glob, re, collections

def rect_dim(li):
    st = float(re.search('x="([^"]+)"', li).group(1))
    end = st + float(re.search('width="([^"]+)"', li).group(1))
    return [st, end]

def get_text(li):
#    <text x="393.9329715061058" y="153" text-anchor="middle" dy="0.3em" fill="#000000" 
# style="font-size: 12px; font-family: arial;">DEAD</text>
    try:
        x = float(re.search('x="([^"]+)"', li).group(1))
    except:
        x = None
    try:
        y = float(re.search('y="([^"]+)"', li).group(1))
    except:
        y = None
    try:
        text = re.search('>([^<]+)</text', li).group(1)
    except:
        text = None
    return x, y, text

def get_rect(fname: str, verbose: bool=False):
    lines = open(fname, 'r').readlines()
    n_rect = 0
    gene = {'domains':[], 'texts':[],
           'name':re.search('lolipop/(.+)_lollipop.svg', fname).group(1)}
    for n, li in enumerate(lines):
        # First rect is the background.
        # Second is the gene body.
        #if li[:5] == '<rect' or '<g><rect' in li:
        if '<rect' in li:
            #print(li)
            n_rect += 1
            if n_rect == 2:
                (gene_start, gene_end) = rect_dim(li)
                gene['Start'] = gene_start
                gene['End'] = gene_end
            elif n_rect == 1:
                continue  # Background.
            else:
                gene['domains'].append(rect_dim(li))
        if '<text' in li:
            x,y,text = get_text(li)
            if x is not None and y is not None and text is not None:

                gene['texts'].append([x, y, text])
    gene['texts'] = [x for x in gene['texts'] if (
        'Mutations' not in x[-1]) and (not re.match('\d+', x[-1]) and (
        not re.match('\w\d+.+', x[-1])))]
    
    if verbose:
        print(gene)
    
    # make copies
    texts = [x[:] for x in gene['texts']]
    domain_cp = [_[:] for _ in gene['domains']]
    for x,y,text in texts:
        for j,w in enumerate(domain_cp):

            (rect_a, rect_b) = w
            if rect_a < x < rect_b:
                gene['domains'][j].append(text)
    
    outlines = []
    for dom in gene['domains']:
        frac_a = (dom[0]-gene['Start'])/(gene['End']-gene['Start'])
        frac_b = (dom[1]-gene['Start'])/(gene['End']-gene['Start'])
        outlines.append(f"{gene['name']}\t{dom[-1]}\t{frac_a}\t{frac_b}\n")

    if verbose:
        print(''.join(outlines))
        
    return gene, outlines
