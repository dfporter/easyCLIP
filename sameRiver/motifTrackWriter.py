import re, HTSeq, os

class motifTrackWriter:

    @staticmethod
    def rc(_):
        d = {'A': 'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        return ''.join([d.get(base.upper(), 'N') for base in _[::-1]])
    
    @classmethod
    def add_motif_locations(cls, chrom, seq, strand, motif, ga, value=1):

        m = re.finditer(motif, seq, re.IGNORECASE)
        for v in m:
            iv = HTSeq.GenomicInterval(chrom, v.span()[0], v.span()[1], strand)
            ga[iv] += value

    @classmethod
    def write_motif_track(cls, genomic_fasta, motif, motif2=None):

        ga = HTSeq.GenomicArray('auto', stranded=True)
        
        for chrom, seq in genomic_fasta.items():
            print(chrom, )

            cls.add_motif_locations(chrom, seq, '+', motif, ga)

            if motif2 is not None:
                cls.add_motif_locations(chrom, seq, '+', motif2, ga, value=2)

            cls.add_motif_locations(chrom, seq, '-', cls.rc(motif), ga)

            if motif2 is not None:
                cls.add_motif_locations(chrom, seq, '-', cls.rc(motif2), ga, value=2)
        
        cls._mk('beds/')
        cls._mk('beds/motif_bedgraphs/')

        ga.write_bedgraph_file('beds/motif_bedgraphs/{}_+.wig'.format(motif), '+')
        ga.write_bedgraph_file('beds/motif_bedgraphs/{}_-.wig'.format(motif), '-')

    @staticmethod
    def _mk(_dir):
        if not os.path.exists(_dir):
            os.system('mkdir ' + _dir)