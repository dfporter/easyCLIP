# Readme

Peak and motif identification is done using two notebooks:

motif_and_locations.ipynb:
 Find peaks, output locations, signal at locations, significance values for signal
 vs controls.

write_fasta_for_homer.ipynb:
 Apply cutoffs to determine which peaks to include, then write fasta file of the
 sequences under each peak for HOMER analysis.

The P value assignment process is CRITICAL, EVEN IF NO P VALUE IS USED, because
 during the calculation we apply a minimum read cutoff that effectively drops
 from any further consideration RNAs that don't reach that cutoff.
 This is a highly non-ideal way to do things.
 Specifically, in statsForCountsNB.calculate_pvalues() there is a parameter
 cutoff_for_consideration such that for each RNA, if the maximum "count" (per million
 reads or per billion/ten billion proteins) across all positive/negative replicates
 is below the cutoff, the RNA is not included in the output table.
 This was done to speed up the calculation, but a value of 5 reads per million was
 used for the original motif identification, both with and without FDR cutoffs,
 and this cutoff was necessary for the expected motif to reach the top spot
 for Rbfox1 and CELF1.

This fact is problematic, because it means the inclusion/exclusion of other
 samples affects the output for a given protein, even if it we are only
 using peak heights.
 At present, we are simply living with this decision.

```bash
export PATH=$PATH:/Users/dfporter/homer/bin/
```
```python
# Write call to HOMER.
cmd = "{fp} {d}/{g}.fa  fasta ./{g} -fasta {d}/randoms_for_{g}.fa  -rna -homer1".format(
        fp='./bin/findMotifs.pl', d=seq_dir, g=prot)
```