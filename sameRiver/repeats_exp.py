
def mapping(self, cutadapt_path: str):
    # Requires bowtie2 to be on path and 'STAR'/'STAR index' to be in self.file_paths.
    # Requires samtools to be on path.

    #not_mapping = """
    sameRiver.remove_short_reads.remove_short_reads_in_paired_fastq(
        cutadapt_path + '/concatenated_R1_reads_for_mapping.fastq',
        cutadapt_path + '/concatenated_R2_reads_for_mapping.fastq',
        cutadapt_path + '/long_concatenated_R1_reads_for_mapping.fastq',
        cutadapt_path + '/long_concatenated_R2_reads_for_mapping.fastq',
        )

    if 'STAR' not in self.file_paths:
        self.file_paths['STAR'] = 'STAR'
        print("Assuming STAR is on PATH.")

    if 'sams' not in self.file_paths:
        self.file_paths['sams'] = self.file_paths['fastq'] + '/sam/'

    for _dir in [self.file_paths['sams'], self.file_paths['beds'],
            self.file_paths['sams'] + '/split/']:

        if not os.path.exists(_dir):
            os.system('mkdir ' + _dir)

    # First map to genome.
    cmd = self.file_paths['STAR']
    cmd += ' --genomeDir {}'.format(self.file_paths['STAR index'])
    cmd += ' --runThreadN 10'
    cmd += ' --readFilesIn {} {}'.format(

        #'{}/repeats.Unmapped.out.mate1'.format(self.file_paths['sams']),
        #'{}/repeats.Unmapped.out.mate2'.format(self.file_paths['sams']),

        cutadapt_path + '/long_concatenated_R1_reads_for_mapping.fastq',
        cutadapt_path + '/long_concatenated_R2_reads_for_mapping.fastq',            
        )
    cmd += ' --outReadsUnmapped Fastx'
    # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
    # of a paired end read) reads in separate file(s).
    # output in separate fasta/fastq files, Unmapped.out.mate1/2
    cmd += ' --outFileNamePrefix {}/genome.'.format(self.file_paths['sams'])

    print(cmd, '\n')
    subprocess.check_output(cmd.split(' '))

    # Then map to repeats.
    # python RepEnrich2.py /data/mm9_repeatmasker.txt /data/Sample_Output_Folder/ \
    # Sample_Name /data/RepEnrich2_setup_mm9 /data/sample_name_multimap_R1.fastq \
    # --fastqfile2 /data/sample_name_multimap_R2.fastq /data/sample_name_unique.bam \
    # --cpus 16 --pairedend TRUE
    args = Object()

    args.annotation_file = 'hg38_repeatmasker_clean.txt'
    outputfolder = 'output_RepEnrich/'
    outputfile_prefix = 'repeats'
    setup_folder = 'hg38re/'
    repeat_bed = setup_folder + os.path.sep + 'repnames.bed' 
    unique_mapper_bam = 'unique_mapped.bam'
    fastqfile_1 = '{}/genome.Unmapped.out.mate1'.format(self.file_paths['sams'])
    fastqfile_2 = '{}/genome.Unmapped.out.mate2'.format(self.file_paths['sams'])
    cpus = 12
    b_opt = "-k 1 -p 1 --quiet"
    simple_repeat = 'Simple_repeat'
    paired_end = 'TRUE'
    allcountmethod = "FALSE"
    is_bed = "FALSE"

    hg38_repeats_path = ''
    cmd = f'bowtie2 -q -x {hg38_path}'  # ../hg38/hg38
    cmd += f' -1 {fa1} -2 {fa2} -S temp.sam'
    
    cmd = 'mv {gen} {gen}.not_mapq_filtered'.format(
        gen=self.file_paths['sams'] + '/genome.Aligned.out.sam')
    subprocess.check_output(cmd.split(' '))

    cmd = 'samtools view -q 10 -h -o {gen} {gen}.not_mapq_filtered'.format(
        gen=self.file_paths['sams'] + '/genome.Aligned.out.sam')
    subprocess.check_output(cmd.split(' '))
    #"""
    
    # Get the SQ and PG lines from the repeats header.
    seqlines = ''
    pg_line = ''
    with open(self.file_paths['sams'] + '/repeats.Aligned.out.sam') as fh:
        for li in fh:
            if li[:3] == '@SQ':
                seqlines += li
            if li[:3] == '@PG':
                pg_line += li
            if li[0] != '@':
                break

    # Get all the header lines from the genomic sam.
    hd_line = ''
    post_seqlines = ''
    with open(self.file_paths['sams'] + '/genome.Aligned.out.sam') as fh:
        for li in fh:
            if li[:3] == '@HD':
                hd_line = li
            elif li[:3] == '@SQ':
                seqlines += li
            elif li[0] == '@':
                post_seqlines += li
            else:
                break

    # Write a combined sam including the combined header and combined sam.
    with open(self.file_paths['sams'] + '/all.Aligned.out.sam', 'w') as fh:
        fh.write(hd_line)
        fh.write(seqlines)
        fh.write(post_seqlines)
        fh.write(pg_line)

        with open(self.file_paths['sams'] + '/genome.Aligned.out.sam') as genomefh:
            for li in genomefh:
                if li[0] == '@':
                    continue
                fh.write(li)

        with open(self.file_paths['sams'] + '/repeats.Aligned.out.sam') as repeatsfh:
            for li in repeatsfh:
                if li[0] == '@':
                    continue
                fh.write(li)

    # Filter out secondary alignments and read mates.
    self.filter_sam(
        self.file_paths['sams'] + '/all.Aligned.out.sam',
        self.file_paths['sams'] + '/all.filtered.sam')

    # Split the giant sam file and make bed and bedgraph files.
    split_sam_filenames = sameRiver.sam_to_bed_and_wig.split_sam(
        self.file_paths['sams'] + '/all.filtered.sam',
        self.file_paths['sams'] + '/split/')

    # Collapse duplicates.
    sam_dups_dir = self.file_paths['sams'] + '/before_duplicates_collapsed/'
    if not os.path.exists(sam_dups_dir):
        os.system('mkdir ' + sam_dups_dir)

    for sam_fname in split_sam_filenames:
        dups_file = '{}/{}'.format(
            sam_dups_dir, re.sub('.sam$', '.with_dups.sam', os.path.basename(sam_fname)))
        os.system('mv {} {}'.format(sam_fname, dups_file))
        sameRiver.collapse_duplicates.collapse_sam(dups_file, sam_fname)

    # Make bed and bedgraph files.
    #for sam_fname in split_sam_filenames:
    for sam_fname in glob.glob(self.file_paths['sams'] + '/split/*sam'):
        sameRiver.sam_to_bed_and_wig.sam_to_bed_and_wig(
            sam_fname,
            '{}/{}.bed'.format(
                self.file_paths['beds'],os.path.basename(sam_fname).split('.sam')[0]),
            self.file_paths['beds']  # Wigs directory.
            )
