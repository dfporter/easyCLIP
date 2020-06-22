import HTSeq, os, importlib

import sameRiver
import sameRiver.bedgraphs

importlib.reload(sameRiver.bedgraphs)
class artifactSplitter:
    
    def __init__(self):
        pass
    
    @staticmethod
    def remove_artifactual_regions_in_chrom_vector(
        ga, chrom_name, chrom_vector, strand='+', dist_cutoff=250, height_cutoff=10):
        
        #previous_read_loc_dists = {'before': 0, 'after': 0, 'at': 0, 'value': 0}
        previous_read_loc_dists = [0, 0, 0, 0]
        #this_read_loc_dists = {'before': 0, 'after': 0, 'at': 0, 'value': 0}
        this_read_loc_dists = [0, 0, 0, 0]
        to_remove = []
        for iv, val in chrom_vector[strand].steps():

            if val == 0:
                continue

            #this_read_loc_dists['at'] = iv.start
            this_read_loc_dists[2] = iv.start
            #this_read_loc_dists['value'] = val
            this_read_loc_dists[3] = val
            #this_read_loc_dists['before'] = iv.start - previous_read_loc_dists['at']
            this_read_loc_dists[0] = iv.start - previous_read_loc_dists[2]

            #previous_read_loc_dists['after'] = iv.start - previous_read_loc_dists['at']
            previous_read_loc_dists[1] = iv.start - previous_read_loc_dists[2]

            if (
                #previous_read_loc_dists['before'] >= dist_cutoff) and (
                previous_read_loc_dists[0] >= dist_cutoff) and (
                #previous_read_loc_dists['after'] >= dist_cutoff) and (
                previous_read_loc_dists[1] >= dist_cutoff) and (
                #previous_read_loc_dists['value'] >= height_cutoff ):
                previous_read_loc_dists[3] >= height_cutoff ):

                #to_remove.append(previous_read_loc_dists['at'])
                to_remove.append(previous_read_loc_dists[2])

            previous_read_loc_dists = this_read_loc_dists.copy()

        for pos in to_remove:
            ga[HTSeq.GenomicInterval(chrom_name, pos, pos+1, strand)] = 0
        
        return ga
    
    @classmethod
    def remove_artifactual_regions(
        cls, ga, strand='both', dist_cutoff=250, height_cutoff=10):
        
        for chrom_name, chrom_vector in ga.chrom_vectors.items():
            print(chrom_name, end=', ')

            for strand in ['+', '-']:
                ga = cls.remove_artifactual_regions_in_chrom_vector(
                    ga, chrom_name, chrom_vector, strand=strand, dist_cutoff=250, height_cutoff=10)
                    
        return ga
    
    @classmethod
    def remove_artifacts_from_bedgraphs_folder(
        cls, bedgraphs_folder, dist_cutoff=250, height_cutoff=10):

        out_folder = os.path.dirname(bedgraphs_folder.rstrip('/')) + '/read_start_filtered/'
        
        print('filtered .wigs to: ', out_folder)

        if not os.path.exists(out_folder):
            os.system('mkdir {0}'.format(out_folder))
            
        bg = sameRiver.bedgraphs.set_of_bedgraphs(
            bedgraphs_folder=bedgraphs_folder)

        for bedfname, bedgraph in bg.bedgraphs.items():
            bedgraph.ga = cls.remove_artifactual_regions(
                bedgraph.ga, dist_cutoff=dist_cutoff, height_cutoff=height_cutoff)
        
        bg.write_bedgraphs(top_dir=out_folder)

    @classmethod
    def remove_artifacts_from_list_of_bedgraphs_folders(
        cls, list_of_bedgraphs_folders, dist_cutoff=250, height_cutoff=10):
        
        for bedgraphs_folder in list_of_bedgraphs_folders:
            cls.remove_artifacts_from_bedgraphs_folder(
                bedgraphs_folder, dist_cutoff=dist_cutoff, height_cutoff=height_cutoff)