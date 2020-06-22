import collections, os, sys
from pathlib import Path


def collapse(
    annotation_file: str='hg38_repeatmasker_clean.txt', output_filename: str='default') -> None:
    """Edit the annotation file (e.g. hg38_repeatmasker_clean.txt) to
    put repeats in the same family together (with exceptions) and output
    the edited annotation file in the same format.
    """

    # Use a default output filename if not given a filename.
    if output_filename is None:
        output_filename = Path(
            os.path.dirname(annotation_file), 'collapsed_by_family_' + os.path.basename(annotation_file))

    outf = open(output_filename, 'w')

    special_repeats = set()

    # For line describing a genomic instance of a repeat:
    with open(annotation_file) as fh:

        # The header is split into two lines for visualization, hence the weird column names.
        header = next(fh)
        outf.write(header)

        header = header.split()

        header2 = next(fh)
        outf.write(header2)

        blank = next(fh)
        outf.write(blank)

        col_number_for_repeat = header.index('matching') + 1  # It's a weird format.
        col_number_for_family = header.index('repeat') + 1

        # Count instances.

        kept_alu = ['AluJb', 'AluSx1', 'AluSx', 'AluY']

        for n, li in enumerate(fh):

            s = li.split()

            # Example column entries: DNA/hAT-Charlie, SINE/MIR.
            family = s[col_number_for_family].split('/')[-1]

            # Some repeat families we leave as they are.
            # The Alu and L1 genomes end up too large to be bai indexed,
            # so we leave out some of the largest members of these groups to be counted
            # separately:
            
            if family=='L1':
                # These are three subfamilies of L1 elements, see:
                # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030137
                # Evolutionary History of Mammalian Transposons Determined by Genome-Wide Defragmentation
                if 'L1ME' in s[col_number_for_repeat]:  # L1M stands for L1 mammalian.
                    s[col_number_for_repeat] = 'L1ME'
                elif 'L1MB' in s[col_number_for_repeat]:
                    s[col_number_for_repeat] = 'L1MB'
                elif 'L1PA' in s[col_number_for_repeat]:  # L1P stands for L1 primate.
                    s[col_number_for_repeat] = 'L1PA'
                else:
                    s[col_number_for_repeat] = 'L1'

            elif family=='Alu':
                # For Alu, we just keep the top 4 types separate; the rest are set as just 'Alu'.
                if s[col_number_for_repeat] in kept_alu:
                    pass  # Don't combine these into their family.
                else:
                    s[col_number_for_repeat] = s[col_number_for_family].split('/')[-1]   

            elif family in ['rRNA', 'scRNA', 'snRNA', 'tRNA']:
                pass  # Don't combine these into their family.

            else:  # Simplify all other repeats to be just their family.
                s[col_number_for_repeat] = s[col_number_for_family].split('/')[-1]                

            li = ' '.join(s)

            outf.write(li + '\n')

    outf.close()


def usage():
    print(collapse.__doc__)
    print("USAGE: python compress_repeat_genomes.py input_annotation_filename [output_filename]")
    print("Generally: python compress_repeat_genomes.py hg38_repeatmasker_clean.txt collapsed_by_family_hg38_repeatmasker_clean.txt")
    print("The default output filename is the input with a collapsed_by_family_ prefix.")   


if __name__ == '__main__':

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        usage()
    
    elif len(sys.argv) == 2:
        collapse(annotation_file=sys.argv[1])

    elif len(sys.argv) == 3:
        collapse(annotation_file=sys.argv[1], output_filename=sys.argv[2])
