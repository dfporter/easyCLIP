
top_dir = '/oak/stanford/groups/khavari/users/dfporter/seq/all/'
top_dir = '/Users/dfporter/pma/dataAndScripts/clip/miseq/Runs/hiseq_pcbp1_190416/'
top_dir = '/Users/dp/pma/dataAndScripts/clip/meta/'

scheme_file_with_random_proteins = top_dir + '/scheme.xlsx'
ann_counts_file = top_dir + '/ann_counts.txt'
bed_file_dir = top_dir + '/beds/'

random_proteins = [
    'CAPNS2',
    #'CCIN', # Dataset too small.
    'CDK4', 'CHMP3',
    'DCTN6',
    #'EPB41L5',  # Dataset too small.
    'ETS2', 'IDE',
    'ITPA', 'TPGS2', 'UBA2',
    ]

random_protein_bed_file_names = [
    'Exp31_CDK4_05GCCATG_AAC.bed',
    'Exp31_TPGS2_05CCACTC_AAC.bed',
    'Exp31_TPGS2_05CCACTC_CAG.bed',
    'Exp31_UBA2_05TGAGTG_AAC.bed',
    'Exp32_CDK4_05GCCATG_TCA.bed',
    'Exp32_ITPA_05GCCTAA_AGT.bed',
    'Exp32_ITPA_05GCCTAA_TCA.bed',
    'Exp32_UBA2_05TGAGTG_TCA.bed',
    'Exp33_CHMP3_17TGTTGG_AGT.bed',
    'Exp33_CHMP3_17TGTTGG_TCA.bed',
    'Exp33_DCTN6_17GCTGTA_TCA.bed',
    'Exp33_ETS2_17TGCCGA_AGT.bed',
    'Exp33_ETS2_17TGCCGA_TCA.bed',
    'Exp33_IDE_17GCTCAT_TCA.bed',
    'Exp33_UBA2_17TGAGTG_AGT.bed',
    'Exp56_CCIN_tcCTGATC_AGT.bed',
    'Exp56_CCIN_tcCTGATC_TCA.bed',
    'Exp56_EPB41L5_tcAGCTAG_AGT.bed',
    'Exp56_EPB41L5_tcAGCTAG_TCA.bed',
    
    # 2020 batch.
    'Exp92_UBA2_GCCTAA_AAC.bed',
    'Exp92_ETS2_CGAAAC_AGT.bed',
    'Exp92_EPB41L5_GTAGCC_TCA.bed',

    ]

random_protein_bed_file_locations = [
    bed_file_dir + '/' + x for x in random_protein_bed_file_names]
