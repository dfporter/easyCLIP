import re, collections, pandas, random, importlib, os, HTSeq
import sameRiver.RNA
importlib.reload(sameRiver.RNA)
from sameRiver.RNA import RNA


class set_of_named_mRNAs:
    """Holds a dict of RNA objects, and methods for reading a GTF.
    """
    
    def __init__(self, mRNAs={}):
        self.mRNAs = mRNAs
        self.log = 'set_of_named_mRNAs() initiated.'
        self.log += ' Passed a set of {} length to initiate.\n'.format(len(mRNAs))
        print(self.log)

    def get_name(self, li, translator=None, require_transcript_id=True):
        li = li.rstrip('\n')
        
        name = re.search('transcript_id=* *"*([^;"]+)', li)
        
        if name is not None:
            return name.group(1)
        
        elif require_transcript_id:
            return None
        
        else:
            name = re.search('Name=* *"*([^;"]+)', li)
            if (name is not None) and translator is not None:
                name = translator.get(name.group(1), None)
                if name is None:
                    return None
                else:
                    return name
            else:
                return None
    
    def example(self, n_examples=1, print_str=True):
        li = ''
        if n_examples > 1:
            li += "I contain {0} RNA objects. Here are {1} at random.\n".format(len(self.mRNAs), n_examples)
        else:
            li += "I contain {0} RNA objects. Here is {1} at random.\n".format(len(self.mRNAs), n_examples)
                
        for n in range(n_examples):
            name = random.choice(list(self.mRNAs.keys()))
            li += name + '\n' + str(self.mRNAs[name].__dict__) + '\n'

        if print_str:
            print(li)

        return li
            
    def find_introns(self):
        for name, rna in self.mRNAs.items():
            error = self.mRNAs[name].find_introns()
            if error:
                print(f"Error trying to find introns in {name}: {rna.__dict__}")

    def purge_NT_chrom(self):
        to_del_rna = set()

        for name, rna in self.mRNAs.items():

            new_elements = collections.defaultdict(list)

            for _type, list_of_ivs in rna.elements.items():
                for iv in list_of_ivs:
                    if ((iv.chrom[:3] != 'NT_') and (iv.chrom[:3] != 'NW_')):
                        new_elements[_type].append(iv)

            rna.elements = new_elements

            if ('exon' not in rna.elements):
                to_del_rna.add(name)
            elif (len(rna.elements['exon']) == 0):
                to_del_rna.add(name)

        print("Before NT_ chrom purge {0} RNAs.".format(len(self.mRNAs)))
        for rna in list(to_del_rna):
            del self.mRNAs[rna]
        print("After NT_ chrom purge  {0} RNAs.".format(len(self.mRNAs)))

    def housekeeping(self):
        #self.purge_NT_chrom()

        new_elements = collections.defaultdict(list)
        removed_elements = collections.defaultdict(int)
        initial_elements = collections.defaultdict(int)
        kept_elements = collections.defaultdict(int)


        # Check for errors.
        for name, rna in self.mRNAs.items():

            for _type, list_of_ivs in rna.elements.items():
                #list_of_ivs = list(set(list_of_ivs))
                _sorted = sorted(list_of_ivs, key=lambda x: x.start)
                keep = []
                skip_next = False
                any_overlaps = False

                for n, iv in enumerate(_sorted):
                    if skip_next:
                        skip_next = False
                        removed_elements[_type] += 1
                        continue

                    if n == len(_sorted) - 1:
                        keep.append(iv)
                        break

                    if _sorted[n].start == _sorted[n+1].start:
                        if _sorted[n].end >= _sorted[n+1].end:
                            keep.append(iv)
                        else:
                            keep.append(_sorted[n+1])
                        skip_next = True
                        continue
                    else:
                        keep.append(iv)

                kept_elements[_type] += len(keep)
                initial_elements[_type] += len(_sorted)
                rna.elements[_type] = keep
                
        print("""Housekeeping():
        Went through {0} RNAs and removed those elements (exons, CDS, ...) with .start values the same.
        There should not be anything removed here if the GTF parsing has done its job and the GTF has no errors.
        Initial elements {1}.
        Kept elements    {2}.
        Removed elements {3}.""".format(len(self.mRNAs), initial_elements, kept_elements, removed_elements))

    def one_transcript_per_gene(self, biotype_mapping_file='./temp/enst_transcript_id_name_biotype_map.txt'):
        print("Keeping only one transcript version per gene (by length).")
        print(f"Number of transcripts before collapsing: {len(self.mRNAs)}")
        
        to_del = set()
        collapsed = {}
        txpt_to_len = collections.defaultdict(dict)
        txpt_to_type = collections.defaultdict(dict)

        df = pandas.read_csv(biotype_mapping_file, sep='\t', index_col=False)
        transcript_to_type = dict(zip(df['transcript_id'].tolist(), df['gene_name'].tolist()))

        for txpt_id in self.mRNAs:

#            print('one_transcript_per_gene(): name in self.mRNAs =', name)
            
            #symbol = transcript_to_gene.get(name, name)
            symbol = self.mRNAs[txpt_id].gene_names
            if len(symbol) > 1:
                print(symbol, ': more than one gene id for this transcript {}'.format(txpt_id))

#            print('...symbol from transcript_to_gene=', symbol)

            if txpt_id == 'ENST00000383925':
                print("set_of_named_mRNAs.one_transcript_per_gene(): ENST00000383925 present, symbol ", symbol)
            
            if len(self.mRNAs[txpt_id].elements['exon']) == 0:
                continue
            
            _min = min([iv.start for iv in self.mRNAs[txpt_id].elements['exon']])
            _max = max([iv.end for iv in self.mRNAs[txpt_id].elements['exon']])

            txpt_to_len[list(symbol)[0]][txpt_id] = _max - _min

            txpt_to_type[list(symbol)[0]][txpt_id] = transcript_to_type[txpt_id]#self.mRNAs[txpt_id].transcript_biotypes

#            print('...txpt_to_len=', txpt_to_len)

            if txpt_id == 'ENST00000383925':
                
                print(
                    "set_of_named_mRNAs.one_transcript_per_gene(): ENST00000383925 len ",
                    _max - _min)
                
        print("Genes: {0}. Average txpts/gene: {1}.".format(len(txpt_to_len), len(self.mRNAs)/len(txpt_to_len)))
        
        biotype_priority = ['snoRNA', 'snRNA', 'protein_coding']
        
        def priority(symbol):
            for txpt_id, biotype in txpt_to_type[symbol].items():
                for biotype_to_cf in biotype_priority:
                    if biotype == biotype_to_cf:
                        return biotype_to_cf  # Return the first match by priority.
            return False  # This gene symbol has no biotype with a priority to it.


        longest_txpt = {}
        for symbol in txpt_to_len:
            txpt_by_biotype = {}            

            if not priority(symbol): # If this gene symbol has no biotype with a priority to it:
                longest_txpt[symbol] = sorted(txpt_to_len[symbol].keys(), key=lambda x: txpt_to_len[symbol][x])[-1]
                continue

            for txpt_id in txpt_to_type[symbol]:
                biotype = txpt_to_type[symbol][txpt_id]
                txpt_by_biotype.setdefault(biotype, [])
                txpt_by_biotype[biotype].append(txpt_id)

            for biotype in biotype_priority:
                if biotype not in txpt_by_biotype:
                    continue
                longest_txpt[symbol] = sorted(
                    txpt_by_biotype[biotype], key=lambda x: txpt_to_len[symbol][x])[-1]
                break
            
        for symbol, txpt in longest_txpt.items():
            #print('...collapsing symbol {} to only txpt {}.'.format(symbol, txpt))
            collapsed[symbol] = self.mRNAs[txpt]
            #print('...now ', collapsed)
            
        self.mRNAs = collapsed
        if 'RNU1-1' in self.mRNAs:
            print("""set_of_named_mRNAs.one_transcript_per_gene(): 
                  after collapsing, RNU1-1 present, symbol """, symbol)     
                  
        print("Transcripts after collapsing:  {0}".format(len(self.mRNAs)))
            
    def add_gtf_line(self, line, name=None, translator=None, require_transcript_id=True,
                    require_col_2_value=None):
        
        # Try to find a transcript id in this gtf line.
        (name is None) and (name := self.get_name(line))

        if name is None:
            return

        # Initialize an RNA() object.
        if name not in self.mRNAs:
            self.mRNAs[name] = RNA()
            self.mRNAs[name].gene_ids = set()
            self.mRNAs[name].gene_names = set()
            self.mRNAs[name].gene_biotypes = set()
            self.mRNAs[name].transcript_biotypes = set()

        # Add info from this line.
        if (m := re.search('gene_id "([^"]+)"', line)) is not None:
            self.mRNAs[name].gene_ids.add(m.group(1))
        else:
            print(f"No gene_id on line {line}")

        if (m := re.search('gene_name "([^"]+)"', line)) is not None:
            self.mRNAs[name].gene_names.add(m.group(1))
        else:
            print(f"No gene_name on line {line}")
        
        if (m := re.search('transcript_type "([^"]+)"', line)) is not None:
            self.mRNAs[name].transcript_biotypes.add(m.group(1))
        
        if (m := re.search('gene_type "([^"]+)"', line)) is not None:
            self.mRNAs[name].gene_biotypes.add(m.group(1))
        
        self.mRNAs[name].add_gtf_line(line, require_col_2_value=require_col_2_value)

    def create_set_of_named_mRNAs(
            self, mRNA_refseq_to_protein_accession=None, gtf_filename=''):
        """Calls self.read_gtf to read a gtf, then changes SNHG::introns to snoRNA biotypes.
        Can also create a refseq accession translator if mRNA_refseq_to_protein_accession is a csv filename.
        self.read_gtf() sets self.mRNAs to {transcript_id: RNA object}.
        """

        self.log += 'Asked to create a set of named RNAs with the following parameters:'
        self.log += 'gtf_filename={} mRNA_refseq_to_protein_accession={}\n'.format(
            gtf_filename, mRNA_refseq_to_protein_accession)

        os.makedirs('./temp', exist_ok=True)
        import sameRiver.biotypeLookupFileMaker as bl
        bl.make_enst_gene_name_biotype_map_file(
            gtf_filename=gtf_filename, out_filename='./temp/enst_transcript_id_name_biotype_map.txt')

        print(f"Inside create_set_of_named_mRNAs(): len(self.mRNAs)={len(self.mRNAs)}")
        if mRNA_refseq_to_protein_accession is not None:
            df = pandas.read_csv(mRNA_refseq_to_protein_accession, sep='\t')
            self.mrna_to_protein = dict(zip(df['RNA_nucleotide_accession.version'].tolist(),
                                  df['protein_accession.version'].tolist()))
            self.protein_to_mrna = dict(zip(self.mrna_to_protein.values(), self.mrna_to_protein.keys()))

            # Read GTF info into the set of mRNAs. Use the mRNA<->protein refseq map defined above.
            self.read_gtf(gtf_filename=gtf_filename, protein_to_mrna=self.protein_to_mrna)

        else:
            self.read_gtf(gtf_filename=gtf_filename, protein_to_mrna=None)

        print(f"Inside create_set_of_named_mRNAs(), have just loaded gtf: len(self.mRNAs)={len(self.mRNAs)}")
        # Expand the boundaries of snoRNA elements by 2 in both directions
        # to avoid reads that belong to the snoRNA but are 2 nt past the border
        # from being assigned to an intron.
        n_adjusted = collections.defaultdict(int)

        def _shift(iv):
            iv.start -= 2
            iv.end += 2
            return iv

        for name, rna in self.mRNAs.items():
            if 'snoRNA' in rna.transcript_biotypes:
                n_adjusted['snoRNA'] += 1
                adjusted = collections.defaultdict(list)
                if 'SNORD10' in rna.gene_names:
                    print(f"expanding {name}: {rna.__dict__}")
                rna.elements['exon'] = [_shift(_iv) for _iv in rna.elements['exon']]
                if 'SNORD10' in rna.gene_names:
                    print(f"expanded {name}: {rna.__dict__}")

        # Finished expanding snoRNA border.

        print(f"Inside create_set_of_named_mRNAs(), have just expanded snoRNA: len(self.mRNAs)={len(self.mRNAs)}")

        print("Expanded by 2 nt each side the following elements to prevent misassignments:")
        print(n_adjusted)

    def read_gtf(self, gtf_filename='', protein_to_mrna=''):
        
        print("Reading GTF file {0}".format(gtf_filename))
        
        fh = open(gtf_filename)
        
        for line_n, li in enumerate(fh):
            
            if li[0] == '#':
                continue
                
            self.add_gtf_line(
                li, translator=protein_to_mrna,
                #require_col_2_value=['havana', 'ensembl_havana', 'ensembl']
                )
            
            if not(line_n % 1E6):
                print(line_n, li)
                
        fh.close()

        self.log += f'From the gtf file, read {line_n} lines and loaded {len(self.mRNAs)} named RNAs.'
        self.log += f' Here are is a call to example():\n{self.example(n_examples=2, print_str=True)}'
    
        self.chromosomes = set()
        
        for _RNA in self.mRNAs.values():
            for element_type in _RNA.elements:
                #if len(_RNA.elements[element_type]) > 0:
                for ele in _RNA.elements[element_type]:
                    self.chromosomes.add(ele.chrom)
        
        print("Chromosomes in gtf: ")
        print(self.chromosomes)

        self.log += 'Chromosomes in the gtf: {}\n'.format(self.chromosomes)
        
    def keep_only_mRNA(self):

        self.log += 'Asked to keep only mRNA in set_of_named_mRNAs.\n'

        to_delete = set()
        for name in self.mRNAs:
            if ('CDS' not in self.mRNAs[name].elements) or (
                len(self.mRNAs[name].elements['CDS'])<1):
                to_delete.add(name)
            if 'exon' not in self.mRNAs[name].elements:
                to_delete.add(name)
        for name in to_delete:
            del self.mRNAs[name]        

    def process(self):
        # Purge non-mRNA.
        self.keep_only_mRNA()
            
        for name in self.mRNAs:
            self.mRNAs[name].process()
            self.mRNAs[name].name = name
            
    def average_region_lengths(self):
        # Determine the size of the average regions to normalize to.
        region_lens = collections.defaultdict(int)
        region_counts = collections.defaultdict(int)
        region_average_lens = {}
        for ca in self.mRNAs.values():
            _this_total_region = collections.defaultdict(int)
            for _type, region in zip(ca.region_types, ca.regions):
                if region is None:
                    continue
                _this_total_region[_type] += region[1] - region[0]
            for _type, total_len in _this_total_region.items():
                region_lens[_type] += total_len
                region_counts[_type] += 1

        for region in region_counts:
            if float(region_counts[region]) != 0:
                region_average_lens[region] = int(float(region_lens[region])/float(region_counts[region]))
            else:
                region_average_lens[region] = 0
                
        print("Average region sizes: {0}".format(region_average_lens))
        self.average_region_lens = region_average_lens
        return self.average_region_lens

    def define_a_genomic_array_of_sets(self):
        print("Running define_a_genomic_array_of_sets()...")
        self.gaos = HTSeq.GenomicArrayOfSets('auto', stranded=True)
        for name, rna in self.mRNAs.items():
            for _type, list_of_ivs in rna.elements.items():
                
                if _type not in ['exon', 'intron']:
                    continue

                for iv in list_of_ivs:
                    self.gaos[iv] += '::'.join([name, _type])
        return self.gaos

    def RNAs_overlapping_snoRNA(self):
        self.overlaps_snoRNA = set()

        for name, rna in self.mRNAs.items():
            if 'snoRNA' not in rna.gene_biotypes:
                continue

            for iv in rna.elements['exon']:
                for step, genes in self.gaos[iv].steps():
                    self.overlaps_snoRNA |= genes

        self.overlaps_snoRNA |= set([x.split('::')[0] for x in list(self.overlaps_snoRNA)])
        return self.overlaps_snoRNA




