# easyCLIP: genome build

For simply building the standard genome:

```bash
# In directory GRCh38.gencode.v28/:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz

# Or:
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
gunzip *

mkdir star_index/

# sjdbOverhang is 75 for 75 read length.
/share/software/user/open/star/2.5.4b/bin/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v28.primary_assembly.annotation.gtf --sjdbOverhang 75

# Or:
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v29.primary_assembly.annotation.gtf --sjdbOverhang 75
```


NEWER:

```bash
git clone https://github.com/dfporter/RepEnrich2.git --branch py3
git clone https://github.com/dfporter/easyCLIP-dev
```
From https://drive.google.com/drive/folders/0B8_2gE04f4QWNmdpWlhaWEYwaHM, download
 hg38_repeatmasker_clean.txt.gz.

```bash
# Get file with the locations of repeat elements.
gunzip hg38_repeatmasker_clean.txt.gz

# Download the genome.
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Mapping to the 1267 genomes outlined in hg38_repeatmasker_clean takes a long time.
# To speed things up, we concatenate some of the repeat families into a single genome.
# The output file (collapsed_by_family_hg38_repeatmasker_clean.txt)
# is in the same format as the input file.
python compress_repeat_genomes.py hg38_repeatmasker_clean.txt collapsed_by_family_hg38_repeatmasker_clean.txt

# Build the genomic bowtie2 indexes, writing to hg38re/.
python RepEnrich2_setup.py --threads 10 collapsed_by_family_hg38_repeatmasker_clean.txt hg38.fa hg38re

# Get a genomic GTF.
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
gunzip gencode.v29.primary_assembly.annotation.gtf.gz

# Subset the gtf to those lines with good transcript support levels (TSL 1 or NA).
python sameRiver/gtf.py gencode.v29.primary_assembly.annotation.gtf
# Writes gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA

```

```python
import sys, importlib
sys.path.append('./easyCLIP-dev')  # Wherever sameRiver is put.
import sameRiver
import sameRiver.mapping
importlib.reload(sameRiver.mapping)

# Get some utilities.
g = sameRiver.mapping.repeatsGenome()

# If needed:
import sameRiver.gtf
sameRiver.gtf.subset_to_only_tsl1_and_NA('gencode.v29.primary_assembly.annotation.gtf')
# Outputs gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA

g.setup_genomes(
	repeats_fasta_directory='RepEnrich2/hg38re/',
	genomic_gtf='gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA',
	igv_output_directory='RepEnrich2/for_igv/',
	repeats_gtf='repeats_as_separate_chroms.gtf',
	combined_gtf='repeats_and_genome.gtf',
	)

# The setup_genomes() function does these these calls:

# Make the .fa and .fa.fai needed to build an IGV genome.
g.make_files_for_igv_genome(repeats_fasta_directory='hg38re/', output_directory='for_igv/')

# Make a GTF for the repeats genome.
g.make_gtf_of_repeats_genome(repeats_fasta_directory='hg38re/', output_filename='repeats_as_separate_chroms.gtf')

g.combine_genomic_and_repeats_gtf(
	repeats_gtf='repeats_as_separate_chroms.gtf',
	genomic_gtf='gencode.v29.primary_assembly.annotation.gtf.exons_only_tsl1andNA',
	combined_gtf='repeats_and_genome.gtf',
)
```

```bash

```

OLDER:

Now for easyCLIP we built a custom repeats chromosome, gtf and STAR index based on
 repeatmasker.org/genomes/hg38/.
 That creation is described in the [repetitive_genome](https://github.com/dfporter/repetitive_genome) git repository.
 It ends with the creation of repeats.fa and repeats.gtf.

When building a genome for RepEnrich:
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# From RepEnrich, obtain: hg38.fa.out.
```

Once the .fa and .gtf files for both the regular genome and the repeats genome 
 were created, we did an additional couple steps.

 First, we combined the .gtf data from the regular genome and the repeats genome.
 
 Then, we filtered out from the regular genome gtf transcripts that did not have either support
 level "1" (pattern of intron/exons is well supported) or "NA" (no introns/alternative splice forms, or a pseudogene). This created the file combined_tsl1andNA.gtf.

 From that file, we used IGV to create combined_tsl1andNA.genome for visualization.

