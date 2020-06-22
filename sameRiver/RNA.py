import HTSeq, collections, re

class RNA():
    """A minimal object to hold the information from a GTF as ivs.
    GTF is 1-based. Everything else is 0-based, so we switch to 0-based
    while reading in the file. HTSeq is [a,b). It appears that
    GTF is actually [a, b) as well despite being 1-based.
    """
    def __init__(self):
        self.elements = collections.defaultdict(list)
    
    def add_gtf_line(self, line, require_col_2_value=None):
        s = line.rstrip('\n').split('\t')
        
        if (require_col_2_value is not None) and (require_col_2_value):
            if s[1] not in list(require_col_2_value):
                print(f'wrong col 2 in GTF. Got {s[1]} instead of {require_col_2_value}.')
                return
            
        #if re.match('chr.*', s[0]):
        #    s[0] = s[0][:3]
        #s[0] = s[0].split(' ')[0]

        # Adjust to 0-based. GTF already [a,b).
        self.elements[s[2]].append(
            HTSeq.GenomicInterval(s[0], int(s[3])-1, int(s[4])-1, s[6]))

    def find_introns(self):
        """Returns false if no error encountered, true otherwise.
        """

        # Don't do this again if introns are already defined.
        if 'intron' in self.elements and len(self.elements['intron']):
            return False
            
        self.elements['intron'] = []
        
        if ('exon' not in self.elements) or (len(self.elements['exon'])==1):
            return False
        
        exons = self.elements['exon']
        exons = sorted(exons, key=lambda x: x.start)
        
        for exon_n, exon in enumerate(exons):
            
            if exon_n >= len(exons) - 1:
                return False

            if exons[exon_n].end >= exons[exon_n+1].start:
                print(f"Exons not intrepretable:", exons[exon_n], exons[exon_n+1], '\n-\n')
                return True
            
            self.elements['intron'].append(HTSeq.GenomicInterval(
                exon.chrom, exons[exon_n].end, exons[exon_n+1].start, exon.strand))

        return False
