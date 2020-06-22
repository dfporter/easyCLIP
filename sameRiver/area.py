# coding: utf-8

import re, pandas, HTSeq, collections
import numpy as np
import matplotlib.patches as patches


class area:
    """A single region that aligns to the genome without gaps.
    """
    
    def set_iv(self, _iv):
        assert(type(_iv) == type(HTSeq.GenomicInterval('I', 0, 10)))
        self.iv = _iv
        self.length = _iv.end - _iv.start

    def total_coverage(self, ga):
        if not(hasattr(self, 'iv')):
            return 0
        if self.iv is None:
            return 0
        if self.iv.start == self.iv.end:
            return 0
        total = 0
        for _iv, value in ga[self.iv].steps():
            total += value
        return total
        
    def normalize_index(self, index, length, normalized_length=50):
        a = index/length
        return int(a * normalized_length)

    def add_coverage_to_normalized_length(self, ga, norm_ga, norm_len=100, cutoff=0):
        self.length = self.iv.end - self.iv.start
        total = self.total_coverage(ga)
        if total == 0:
            return
        for _iv, value in ga[self.iv].steps():
            norm_pos_start = self.normalize_index(_iv.start - self.iv.start, self.length, normalized_length=norm_len)
            norm_pos_end =  self.normalize_index(_iv.end - self.iv.start, self.length, normalized_length=norm_len)
            if value >= cutoff:
                for norm_pos in range(norm_pos_start, norm_pos_end):
                    norm_ga[norm_pos] += value/total

class subarea(area):
    """A single region that aligns to the genome without gaps and is part of a larger,
    also continuous, area.
    """
    
    def set_parent_iv(self, _iv):
        self.parent_iv = _iv
    
    def set_sub_iv_from_relative_coordinates(self, start, end):
        if self.parent_iv.strand == '+':
            genomic_left = self.parent_iv.start + start
            genomic_right = self.parent_iv.start + end
        else:
            genomic_right = self.parent_iv.end - start
            genomic_left = self.parent_iv.end - end
        self.iv = HTSeq.GenomicInterval(self.parent_iv.chrom, genomic_left,
                                        genomic_right, self.parent_iv.strand)
        self.length = self.iv.end - self.iv.start
        if genomic_left > genomic_right:
            print("Error, start greater than end:")
            print(self.parent_iv, start, end, genomic_right, genomic_left)
            
            
class complex_area(area):
    """A single region that aligns to the genome without gaps, and has subdivisions.
    """
    
    def __init__(self):
        self.name = 'complex_area'
        self.iv = None  # Genomic coordinates.
        self.sub_areas = {}  # Genomic coordinates.
        self.regions = {}  # Relative coordinates.
    
    def define_genomic_subareas_from_relative_regions(self):
        """This is for a single continuous parent area, such as one exon.
        This does not deal with things that have gaps in their alignment to the genome.
        """
        if self.iv is None:
            print("complex_area.define_genomic_subareas_from_relative_regions function called, but the \
            genomic interval of this complex_area object has not been set.")
            return
        
        for name, region in self.regions.items():
            if region is None:
                self.sub_areas[name] = None
                continue
            self.sub_areas[name] = subarea()
            self.sub_areas[name].set_parent_iv(self.iv)
            self.sub_areas[name].set_sub_iv_from_relative_coordinates(*region)
            
    def set_regions(self, regions):
        """Set a dict of name=>(start, stop) relative coordinates.
        These are not genomic coordinates."""
        assert(type(regions)==type({}))
        self.regions = regions
        if self.iv is not None:
            self.define_genomic_subareas_from_relative_regions()
            
    def add_coverage(self, ga, normalizer, norm_len=10, cutoff=0):
        if self.iv is None:
            return
        self.add_coverage_to_normalized_length(
                ga, normalizer.ga['Whole area'], norm_len=normalizer.lengths['Whole area'],
                cutoff=cutoff)
        total = self.total_coverage(ga)
        all_regions = list(self.regions.keys())
        if self.iv.strand == '-':
            all_regions = all_regions[::-1]
        for n, name in enumerate(all_regions):
            if name not in normalizer.ga:
                print('name {0} not in normalizer {1}'.format(name, normalizer.ga.keys()))
                continue
            if name in self.sub_areas:
                if self.sub_areas[name] is not None:
                    self.sub_areas[name].add_coverage_to_normalized_length(
                        ga, normalizer.ga[name], norm_len=normalizer.lengths[name], cutoff=cutoff)


class complex_area_set:
    
    def __init__(self, complex_areas=None):
        if complex_areas is None:
            self.complex_areas = []
        else:
            assert(type(complex_areas) == type([]))
            self.complex_areas = complex_areas
    
    def add(self, complex_area):
        self.complex_areas.append(complex_area)
    
    def average_region_lengths(self):
        # Determine the size of the average regions to normalize to.
        region_lens = collections.defaultdict(int)
        region_counts = collections.defaultdict(int)
        region_average_lens = {}
        for ca in self.complex_areas:
            for name, region in ca.regions.items():
                if region is None:
                    continue
                region_lens[name] += region[1] - region[0]
                region_counts[name] += 1

        for region in region_counts:
            if float(region_counts[region]) != 0:
                region_average_lens[region] = int(float(region_lens[region])/float(region_counts[region]))
            else:
                region_average_lens[region] = 0
                
        print("Average region sizes: {0}".format(region_average_lens))
        self.average_region_lens = region_average_lens
        return self.average_region_lens
    
    def read_coordinates(self, coords):
        for ca in self.complex_areas:
            if ca.name in coords:
                ca.set_iv(coords[ca.name])
                ca.define_genomic_subareas_from_relative_regions()
            else:
                print("Tried to read coordinates into {0} but couldn't find the name.".format(ca.name))
        
    def add_coverage(self, ga, normalizer, norm_len=10, cutoff=0):
        for ca in self.complex_areas:
            ca.add_coverage(ga, normalizer, norm_len=norm_len, cutoff=cutoff)