import numpy as np
from sameRiver import area
import collections


class signal_normalizer:
    
    def __init__(self, region_lengths):
        self.lengths = region_lengths
        self.ga = dict([(region, np.zeros(int(ave_len))) for region, ave_len in self.lengths.items()])
        self.added_to_area = collections.defaultdict(int)
    
    def add_region(self, name, length):
        self.lengths[name] = length
        self.ga[name] = np.zeros(int(length))
    
    def put_regions_together(self, order):
        self.combine = []
        for name in order:
            if name in self.ga:
                self.combine += list(self.ga[name])
        return self.combine

    def normalize_regions_by_length(self):
        _area = area.area()
        for name in self.ga:
            k = 100/len(self.ga[name])
            self.ga[name] = [k*x for x in self.ga[name]]
            
    def normalize_index(self, index, length, normalized_length=50):
        a = index/length
        return int(a * normalized_length)
    
    def resize_array_and_add(self, _arr, name):
        
        if name not in self.ga:
            return
        
        input_len = len(_arr)
        self.added_to_area[name] += 1
        for pos, value in enumerate(_arr):
            norm_pos = self.normalize_index(pos, input_len, normalized_length=len(self.ga[name]))
            self.ga[name][norm_pos] += value
            
    def normalize_areas_by_number_of_arrays_added(self):
        for name in self.added_to_area:
            if self.added_to_area[name] == 0:
                self.ga[name] = []
            self.ga[name] = [x/self.added_to_area[name] for x in self.ga[name]]