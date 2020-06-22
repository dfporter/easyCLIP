import HTSeq, collections, pandas, random, os, sys

import sameRiver
import sameRiver.set_of_named_mRNAs
import sameRiver.bedgraphs

class signal_RNAs():
    """Holds genomic annotations (RNAs) and coverage data (signal).
    Comprises two properties:
    .beds : set of bedgraphs object
    .RNAs : set of RNAs object
    Holds a variety of methods for accessing the data.
    """
    def make_serializable(self):
        """Get rid of HTSeq objects, which have no __dict__ and can cause errors when pickling.
        Replace with a dict of {fname1: [[*iv, value], [*iv, value], ...], ..}
        """
        self.beds.make_serializable()
        self.gaos = None

    def recover_bedgraph_objects(self):
        self.beds.recover_bedgraph_objects()
        self.make_genomic_array()

    def __init__(self, set_of_bedgraphs=None, set_of_RNAs=None):
        assert(type(sameRiver.set_of_named_mRNAs.set_of_named_mRNAs()) == type(set_of_RNAs))
        #assert(type(sameRiver.bedgraphs.set_of_bedgraphs()) is type(set_of_bedgraphs))

        self.beds = set_of_bedgraphs
        self.RNAs = set_of_RNAs
    
    def signal_in_element(self, element_name):
        signal = collections.defaultdict(int)
        
        for name, rna in self.RNAs.mRNAs.items():
            for iv in rna.elements[element_name]:
                
                if iv.start == iv.end:
                    print('Zero length interval. Skipping.')
                    continue
                    
                for fname, bed in self.beds.bedgraphs.items():
                    signal[fname] += bed.total_area_under_curve(iv=iv)
                    
        return signal
    
    def exonic_and_intronic_signal(self, bedgraphs='all'):
        # Fix.
        (self.total_exonic, self.total_intronic) = (collections.defaultdict(int), collections.defaultdict(int))
        
        for gene_name, bed_fname in self.counts.items():
            
            gene_type = name.split('::')
            
            if len(gene_type) < 2:
                continue
                
            if gene_type[1] == 'exon':
                self.total_exonic[bed_fname] += value
            if gene_type[1] == 'intron':
                self.total_intronic[bed_fname] += value
                
        return self.total_exonic, self.total_intronic
    

    def make_genomic_array(self, no_clobber=False):
        """Create a GenomicArrayOfSets object holding RNA names.
        """
        
        if no_clobber and (hasattr(self, 'gaos')) and (len(self.gaos)>2):
            print("Genomic array exists and no_clobber on. {} genomic interval lists.".format(
                len(self.elements)))
            return
        
        self.gaos = HTSeq.GenomicArrayOfSets('auto', stranded=True)
        
        print("Creating a genomic array from {} RNAs...".format(len(self.RNAs.mRNAs)))
        
        if os.path.exists('error_log.txt'):
            error_file = open('error_log.txt', 'a')
        else:
            error_file = open('error_log.txt', 'w')
        error_file.write("===Error log begun===\n")

        element_lists = 0
        for name, rna in self.RNAs.mRNAs.items():
            for _type, list_of_ivs in rna.elements.items():
                
                if _type not in ['exon', 'intron']:
                    continue

                element_lists += 1
                for iv in list_of_ivs:
                    try:
                        self.gaos[iv] += '::'.join([name, _type])
                    except:
                        error_file.write("Error processing {}\n".format(iv))
                        error_file.write("Name {} Properties {}\n".format(name, rna.__dict__))
                        error_file.write(str(sys.exc_info()) + '\n')

        print("Created genomic array from {} genomic interval lists that were exons or introns.".format(
            element_lists))
        error_file.write('===Error log ended===\n')
        error_file.close()


    def write_counts_file(
        self, fname='counts.txt', report_introns_separately=False, report_unassigned=True,
        no_clobber=False):
        
        if not hasattr(self, 'counts'):
            self.assign_to_genes()
            
        #if no_clobber and os.path.exists(fname):
        #    print("Not overwriting counts file {} (no_clobber on).".format(fname))
            
        exon_counts = collections.defaultdict(int)
        intron_counts = collections.defaultdict(int)
        
        print("Creating a counts file with {} RNAs/rows.".format(len(self.counts)))
        
        for gene_type, dict_of_beds in self.counts.items():
            gene = gene_type.split('::')[0]

            if report_introns_separately and len(s := gene_type.split('::')) > 1:
                if s[1] == 'intron':
                    intron_counts[gene] = dict_of_beds
                else:
                    exon_counts[gene] = dict_of_beds
            else:
                exon_counts[gene_type] = dict_of_beds
        
        counts_df = pandas.DataFrame.from_dict(exon_counts, orient='index')
        
        print("counts file contains {} rows, {} columns, and is written to {}. Here's a slice:".format(
            len(counts_df.index), len(counts_df.columns), fname))
        
        print(counts_df.iloc[:2,:5])
        
        counts_df.to_csv(fname, sep='\t', index=True)
        
        if report_introns_separately:
            introns_df = pandas.DataFrame.from_dict(intron_counts, orient='index')
            introns_df.to_csv('intron_counts.txt', sep='\t', index=True)
            
    def assign_to_genes(self, **kwargs):
        
        if 'no_clobber' in kwargs and (kwargs['no_clobber']):
            if (hasattr(self, 'counts')) and (len(self.counts) > 2):
                print("Not assigning reads to genes because self.counts already exists.")
                return
            
        self.counts = {'_ambiguous': collections.defaultdict(int),
                      '_no_feature': collections.defaultdict(int)}
        
        if 'verbose' in kwargs and (kwargs['verbose']):
            verbose = True
        else:
            verbose = False
            
        #if verbose:
        print(f"Assigning signal to RNAs in {len(self.beds.bedgraphs)} bedgraphs.")
        
        for name in self.RNAs.mRNAs.keys():
            self.counts[name + '::exon'] = collections.defaultdict(int)
            self.counts[name + '::intron'] = collections.defaultdict(int)

        for fname, bed in self.beds.bedgraphs.items():
            
            if verbose:
                print(f"Assigning signal to RNAs in {fname}...")
            
            bed_chroms = list(bed.ga.chrom_vectors.keys())
            
            (bed_numeric, bed_chr_string) = ({}, {})
            (ga_numeric, ga_chr_string) = ({}, {})
            
            for numeric_chrom_name in [str(x) for x in range(22)]:
                if numeric_chrom_name in bed_chroms:
                    bed_numeric[numeric_chrom_name] = ''
                if numeric_chrom_name in self.RNAs.chromosomes:
                    ga_numeric[numeric_chrom_name] = ''
            
            for chr_name in ['chr' + str(x) for x in range(22)]:
                if chr_name in bed_chroms:
                    bed_chr_string[chr_name] = ''
                if chr_name in self.RNAs.chromosomes:
                    ga_chr_string[chr_name] = ''
            
            translate_bed_numeric = False
            translate_bed_chr_string = False
            
            if (len(bed_numeric) == 0 and len(bed_chr_string) == 0) or (
                len(ga_numeric) == 0 and len(ga_chr_string) == 0):
                pass  # Did not recognize chr names in one or both.
            elif (len(bed_numeric) and len(bed_chr_string)) or (
                len(ga_numeric) and len(ga_chr_string)):
                pass  # Both formats (1, 2, 3 and chr1, chr2) found in one or both.
            elif len(bed_numeric) > 0 and len(ga_numeric):
                pass  # Both numeric.
            elif len(bed_chr_string) and len(ga_chr_string):
                pass  # Both in the format chr1, chr2...
            elif len(bed_numeric) and len(ga_chr_string):
                # Translate bed file coordinates to chr1/chr2 format.
                translate_bed_numeric = True
            elif len(bed_chr_string) and len(ga_numeric):
                # Translate bed file chr1/chr2 format to numeric.
                translate_bed_chr_string = True
            
            for iv, value in bed.ga.steps():
                
                # Skip intervals with no signal:
                if value <= 0:
                    continue

                # Roll a dice for reminder print statements.
                _print = bool(random.randint(0, 1E6) == 1)

                # Translate chromosome names in the bedgraph to match chr\d if asked.
                if translate_bed_numeric and (iv.chrom in bed_numeric):
                    _print and (print(f"{iv.chrom} -> chr{iv.chrom} (this note pops up every ~1/1E6 times)"))
                    iv.chrom = 'chr{}'.format(iv.chrom)
                
                # Translate chromosomes in the bedgraph chr\d -> \d if asked.
                if translate_bed_chr_string and (iv.chrom in bed_chr_string):
                    _print and (print(f"Translated {iv.chrom} to {iv.chrom.split('chr')[-1]}"))
                    iv.chrom = iv.chrom.split('chr')[-1]
                
                # Just a reminder about this decision.
                if iv.chrom == 'repeats':
                    _print and print("Ignoring strandedness when mapping to repeat elements. (this notification pops up every ~1/1E6 times)".format(iv.chrom, iv.chrom))
                    iv.strand = '+'  # Ignore strandedness when mapping to repeat elements.

                # If this interval has a length >1, it may overlap different genes at different positions.
                # So we walk through it by its overlaps:
                for gene_iv, genes in self.gaos[iv].steps():
                    #print(iv, f' (signal {value}) > gene_iv {gene_iv}, genes overlapping {genes}')

                    if len(genes) == 1:
                        self.counts[list(genes)[0]][fname] += value * (gene_iv.end - gene_iv.start)

                    if len(genes) > 1:
                        # If one gene is an intron and the other an exon, take the exon. Otherwise it's ambiguous.

                        _type_counts = collections.defaultdict(set)
                        
                        for gene in genes:
                            
                            _type_name = gene.split('::')
                            
                            if len(_type_name) < 2:
                                continue
                            
                            _type_counts[_type_name[1]].add(gene)
                            
                        if ('exon' in _type_counts) and (len(_type_counts['exon']) == 1):
                            self.counts[list(_type_counts['exon'])[0]][fname] += value * (gene_iv.end - gene_iv.start)
                        else:
                            self.counts['_ambiguous'][fname] += value * (gene_iv.end - gene_iv.start)

                    if len(genes) == 0:
                        self.counts['_no_feature'][fname] += value * (gene_iv.end - gene_iv.start)
                
                    #for gene in genes:
                    #    print(f'after adding the above: gene {gene} ', self.counts[gene])
