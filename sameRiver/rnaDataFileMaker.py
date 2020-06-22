import os, pickle, importlib

import sameRiver
import sameRiver.RNA
import sameRiver.set_of_named_mRNAs
#import sameRiver.gapped_area_mRNA

importlib.reload(sameRiver.set_of_named_mRNAs)
#importlib.reload(sameRiver.gapped_area_mRNA)

class rnaDataFileMaker():
    
    def __init__(self):
        pass
    
    def make_from_gtf_file(
        self, gtf_filename='/opt/genome/Homo_sapiens.GRCh38.83/',
        one_transcript_per_gene=True,
        output_data_filename='./data/rna.data'):

        # The mRNAs={} is actually needed to prevent some weird (namespace?) bug
        # when calling this function multiple times from Ipython.
        RNAs = sameRiver.set_of_named_mRNAs.set_of_named_mRNAs(mRNAs={})
        RNAs.create_set_of_named_mRNAs(gtf_filename=gtf_filename)

        print(f"Read in {len(RNAs.mRNAs)} RNAs.")
        
        RNAs.housekeeping()

        if one_transcript_per_gene:
            RNAs.one_transcript_per_gene()

        print('Made one transcript per gene.')
        
        RNAs.find_introns()
        os.makedirs(os.path.dirname(output_data_filename), exist_ok=True)
        pickle.dump(RNAs, file=open(output_data_filename, 'wb'))
        
        print(f"Wrote data file to {output_data_filename}")

        with open('data/rna_data_file_creation_log.txt', 'w') as f:
            f.write(RNAs.log)
        
        return RNAs