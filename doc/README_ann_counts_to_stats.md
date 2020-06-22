# easyCLIP ann_counts.txt files to statistical results and per-protein/per-read counts


Here's a typical workflow for analyzing a set of positives against
 some negatives, both having come from experiments that genereated
 annotated counts files with exp.annotate_counts_file().

```python
import importlib
import sameRiver
import sameRiver.metadata
import sameRiver.metadata.negative_metadata as negative_metadata
import sameRiver.metadata.positive_metadata_pcbp1_190416 as positive_metadata
import sameRiver.negativeCounts
import sameRiver.positiveCounts
import sameRiver.scheme
import sameRiver.statsForCountsNB
importlib.reload(sameRiver.statsForCountsNB)
importlib.reload(sameRiver.scheme)

# If never run before:
negatives = sameRiver.negativeCounts.negativeCounts(negative_metadata)
negatives.normalize_ann_counts()
negatives.get_lowest_values_for_all_proteins()
# Optional: write_txt=True to write some txt's of the data.
negatives.save(write_object=True, write_txt=True)

# If loading:
negatives = sameRiver.negativeCounts.negativeCounts.load()

# If never run before:
positives = sameRiver.positiveCounts.positiveCounts(positive_metadata)
positives.normalize_ann_counts()
positives.save(write_object=True, write_txt=True)

# If loading:
positives = sameRiver.positiveCounts.positiveCounts.load()

#nb = sameRiver.statsForCountsNB.statsForCountsNB.load()
nb = sameRiver.statsForCountsNB.statsForCountsNB(
    negatives=negatives, positives=positives)
nb.calculate_pvalues(which='per_read')
nb.calculate_pvalues(which='per_protein')
nb.save()
nb.write_targets(which='per_read', outfname='default')
nb.write_targets(which='per_protein', outfname='default')



```

Here's a positive_metadata.py file:
```python


top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/Runs/hiseq_pcbp1_190416/'
scheme_file = top_dir + '/scheme.xlsx'
ann_counts_file = top_dir + '/ann_counts.txt'
bed_file_dir = top_dir + '/beds/'

positive_proteins = [
    'PCBP1', 'PCBP1:100P', 'PCBP1:100Q',
    'PCBP1:dKH', 'CELF1',
]
```

Here's a negative_metadata.py file:
```python
top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/all/'

scheme_file_with_random_proteins = top_dir + '/scheme.xlsx'
ann_counts_file = top_dir + '/ann_counts.txt'
bed_file_dir = top_dir + '/beds/'

random_proteins = [
    'CAPNS2', 'CCIN', 'CDK4', 'CHMP3',
    'DCTN6', 'EPB41L5', 'ETS2', 'IDE',
    'ITPA', 'TPGS2', 'UBA2',
    ]
```


Here's another way to define metadata for analysis.py scripts to create targets/stats:
```python
class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/v2all/'
gene = 'hnRNPC'
positive_metadata = Namespace(
    top_dir = top_dir,
    scheme_file = top_dir + '/scheme.xlsx',
    ann_counts_file = top_dir + '/data/signal_at_{}_peak_locations.xlsx'.format(gene),
    bed_file_dir = top_dir + '/beds/',
    positive_proteins = [
        'CAPNS2', 'CCIN', 'CDK4', 'CHMP3',
        'DCTN6', 'EPB41L5', 'ETS2', 'IDE',
        'ITPA', 'TPGS2', 'UBA2',
        'FBL', 'HCT116', 'hnRNPC',
        'SF3B1',
        'PCBP1', 'PCBP1:100P', 'PCBP1:100Q',
        'PCBP1:dKH', 'CELF1',
        'Rbfox1', 'Rbfox2', 'hnRNPD',
    ]
)

negative_metadata = Namespace(
    top_dir = top_dir,
    scheme_file_with_random_proteins = top_dir + '/scheme.xlsx',
    ann_counts_file = top_dir + '/data/signal_at_{}_peak_locations.xlsx'.format(gene),
    bed_file_dir = top_dir + '/beds/',
    random_proteins = [
        'CAPNS2', #'CCIN',
        'CDK4', 'CHMP3',
        'DCTN6', #'EPB41L5',
        'ETS2', 'IDE',
        'ITPA', 'TPGS2', 'UBA2',
    ])

```