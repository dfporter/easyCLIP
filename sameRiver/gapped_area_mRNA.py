import re, collections, HTSeq
import sameRiver
import sameRiver.gapped_area

class gapped_area_mRNA(sameRiver.gapped_area.gapped_area):
    
    def __init__(self):
        super().__init__()
        self.bedgraphs = {}
        self.raw_signal_by_type = collections.defaultdict(list)
        
    def read_coverage(self, ga, cutoff=0):
        pass
    
    def add_bedgraph_signal(self, ga, name='Signal'):
        self.bedgraphs[name] = HTSeq.GenomicArray('auto', stranded=True)
        
        for n, region in enumerate(self.regions):  # For each relative iv. Probably per exon.
            _area = self.areas[n]  # Genomic iv. This exon.
            _type = self.region_types[n]  # The type of this region.
            if _area.start == _area.end:
                continue
                
            for _iv, value in ga[_area].steps():  
                self.bedgraphs[name][_iv] = value
                self.raw_signal_by_type[_type] += [value] * (_iv.end - _iv.start)
    
    def zero(self):
        if hasattr(self, 'bedgraphs'):
            for name in self.bedgraphs:
                self.bedgraphs[name] = HTSeq.GenomicArray('auto', stranded=True)
        
        if hasattr(self, 'raw_signal_by_type'):
            for _type in self.raw_signal_by_type:
                self.raw_signal_by_type[_type] = []
                
    def add_coverage_from_mRNA(self, ga, normalizer, cutoff=0, keep_raw_signal=True):
        """Add to the passed ga object a normalized coverage distribution."""
        
        self.add_bedgraph_signal(ga, name='Signal')
                
        total = sum([sum(_arr) for _arr in self.raw_signal_by_type.values()])
        if total == 0:
            return
        
        if total < cutoff:
            return
        
        for _type, value in self.raw_signal_by_type.items():
            
            if _type not in normalizer.ga:
                continue
                
            _arr = [x/total for x in value]
            
            if self.areas[0].strand == '-':
                _arr = _arr[::-1]
                _type = {'left_UTR': 'right_UTR', 'right_UTR': 'left_UTR', 'CDS': 'CDS'}[_type]
                
            if len(_arr) < 20:
                return
            
            normalizer.resize_array_and_add(_arr, _type)
        
        norm_total = 0
        for _arr in self.raw_signal_by_type.values():
            norm_total += sum(_arr)
        
    def process(self):
        self.make()  # Defines self.exons as dict of ivs, key=exon num counting from genomic left.
        
        cds = self.elements['CDS']  # List of ivs.
        cds = sorted(cds, key=lambda x: x.start)
        
        self.define_relative_regions_from_genomic_subareas()

        cds_left = self.relative_from_genomic(cds[0].start)
        cds_right = self.relative_from_genomic(cds[-1].end)

        pos_in_rel_seq = 0
        self.region_types = []  # List of region types, which correspond in order to areas.
        self.genomic_start_pos_of_region = []
        self.areas = []  # List of genomic ivs, which correspond in order to region_types.
        
        for n, exon in self.exons.items():

            exon_left = self.relative_from_genomic(exon.iv.start)
            exon_right = self.relative_from_genomic(exon.iv.end)

            if (cds_left == exon_left) and (cds_right == exon_right):
                self.regions.append((exon_left, exon_right))
                self.region_types.append('CDS')
                self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, exon.iv.end, exon.iv.strand))
            elif exon_right <= cds_left:
                self.regions.append((exon_left, exon_right))
                self.region_types.append('left_UTR')
                self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, exon.iv.end, exon.iv.strand))
            elif exon_left >= cds_right:
                self.regions.append((exon_left, exon_right))
                self.region_types.append('right_UTR')
                self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, exon.iv.end, exon.iv.strand))
            else: # The exon overlaps with the CDS, but not exactly overlaps.
                if (exon_left < cds_left < exon_right) and (exon_left < cds_right < exon_right):
                    # CDS is within exon.
                    self.regions.append((exon_left, cds_left))
                    self.region_types.append('left_UTR')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, cds[0].start, exon.iv.strand))
                    self.regions.append((cds_left, cds_right))
                    self.region_types.append('CDS')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        cds[0].start, cds[-1].end, exon.iv.strand))
                    self.regions.append((cds_right, exon_right))
                    self.region_types.append('right_UTR')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        cds[-1].end, exon.iv.end, exon.iv.strand))
                elif (cds_left < exon_left) and (cds_right > exon_right):
                    # CDS surrounds exon.
                    self.regions.append((exon_left, exon_right))
                    self.region_types.append('CDS')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, exon.iv.end, exon.iv.strand))
                elif (exon_left < cds_left < exon_right):
                    # CDS starts within exon.
                    self.regions.append((exon_left, cds_left))
                    self.region_types.append('left_UTR')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, cds[0].start, exon.iv.strand))
                    self.regions.append((cds_left, exon_right))
                    self.region_types.append('CDS')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        cds[0].start, exon.iv.end, exon.iv.strand))                    
                elif (exon_left < cds_right < exon_right):
                    # CDS ends within exon.
                    self.regions.append((exon_left, cds_right))
                    self.region_types.append('CDS')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom,
                                        exon.iv.start, cds[-1].end, exon.iv.strand))
                    self.regions.append((cds_right, exon_right))
                    self.region_types.append('right_UTR')
                    self.areas.append(HTSeq.GenomicInterval(exon.iv.chrom, 
                                        cds[-1].end, exon.iv.end, exon.iv.strand))