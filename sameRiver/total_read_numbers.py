from typing import Mapping, List
import os, glob

def for_split_bedgraph_lines(fname: str) -> List:
    with open(fname) as f:
        next(f)  # Skip header.
        for li in f:
            sp = li.rstrip('\n').split('\t')
            yield (sp[0], int(sp[1]) , int(sp[2]), sp[3])


def get_auc(fname: str) -> float:
    auc = 0
    for s in for_split_bedgraph_lines(fname):
        auc += (s[2] - s[1]) * float(s[3])
    
    return auc

def total_read_numbers(folder: str, outfile='./data/') -> Mapping[str, float]:
    """Total area under the curve for bedgraphs in a folder.
    Assuming each read is a single point with value 1, this is the total read number.
    Write the results to a file.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    
    aucs = {}
    for fname in glob.glob(f"{folder}/*_+.wig"):
        print(f"{fname}...", end=" ")
        aucs[fname] = get_auc(fname)
        minus_fname = fname.split('+.wig')[0] + '-.wig'
        aucs[fname] += get_auc(minus_fname)
        print(f"AUC={aucs[fname]:,}")
    
    aucs = {os.path.basename(k).rstrip('_+.wig'): v for k,v in aucs.items()}
    
    with open(outfile, 'w') as f:
        f.write("Dataset\tTotal read number\n")
        [f.write(f"{name}\t{val}\n") for name, val in aucs.items()]
        
    return aucs