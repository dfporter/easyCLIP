import sameRiver.area
import collections, HTSeq


class gapped_area(sameRiver.area.complex_area):
    """An area that aligns to the genome with gaps, and each region of alignment
    is an area object. The area objects in an a gapped_area can be given
    a named type, such as CDS.
    After exons are read in, they may be reorganized into regions, and two sets of coordinates
    corresponding to exactly the same places are made:
    self.regions = relative coordinates
    self.areas = genomic coordinates of the areas
    """

    def __init__(self):
        super().__init__()
        self.elements = collections.defaultdict(list)
        self.exons = {}  # Genomic coordinates.
        self.regions = []  # Relative coordinates.
        self.region_types = []
        self.genomic_left_pos_of_region = []
        self.areas = []  # Genomic coordinates that correspond 1-to-1 with the regions.
        
    def add_gtf_line(self, line):
        s = line.rstrip('\n').split('\t')
        self.elements[s[2]].append(
            HTSeq.GenomicInterval(s[0], int(s[3]), int(s[4]), s[6]))

    def define_relative_regions_from_genomic_subareas(self):
        """Build a function genomic coordinate -> relative coordinate."""
        self.length_of_leftward_exons = {}
        for n, exon in self.exons.items():
            if n == 0:
                self.length_of_leftward_exons[0] = 0
                continue
            self.length_of_leftward_exons[n] = self.length_of_leftward_exons[n-1]
            self.length_of_leftward_exons[n] += self.exons[n-1].iv.end - self.exons[n-1].iv.start
            
    def relative_from_genomic(self, pos):
        if not (hasattr(self, 'length_of_leftward_exons')):
            self.define_relative_regions_from_genomic_subareas()
            
        for n, exon in self.exons.items():
            if exon.iv.start <= pos <= exon.iv.end:
                return self.length_of_leftward_exons[n] + pos - exon.iv.start
        return 0
    
    def genomic_from_relative(self, pos):
        for n, region in enumerate(self.regions):
            if region[0] <= pos <= region[1]:
                return self.genomic_left_pos_of_region[n] + pos - region[0]
    
    def make(self):
        """Creates a dict self.exons of area objects, each with the iv of the exon.
        """
        exons = self.elements['exon']
        exons = sorted(exons, key=lambda x: x.start)
        for n, exon in enumerate(exons):
            self.exons[n] = sameRiver.area.area()
            self.exons[n].set_iv(exon)
