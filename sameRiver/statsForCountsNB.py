import sameRiver.statsForCounts
import collections, random, pandas, os, re, importlib
import numpy as np
import scipy.stats as stats
import statsmodels
import statsmodels.api as sm
from sameRiver.biotypeAdder import biotypeAdder
from typing import Union, Mapping, List

class statsForCountsNB(sameRiver.statsForCounts.statsForCounts):
    """
    Input is an ann_counts.txt and a percentCrosslinked.xlsx file.
    The outputs from statsForCountsNB generally used are:
    tables/pvals_per_read.xlsx, and tables/pvals_per_protein.xlsx
    """
    def make_df_averaged_across_replicates(self, df, scheme):
        """Much faster. Not used."""

        by_prot = collections.defaultdict(list)
        for col in self.positives.numeric_columns(df):
            prot = scheme.gene_from_fname(col)
            by_prot[prot].append(df[col].values)

        for prot in by_prot:
            arrs = by_prot[prot]
            n = len(by_prot[prot])
            by_prot[prot] = np.sum(arrs, axis=0)/n

        by_prot['gene_name'] = df['gene_name']
        by_prot['Gene type'] = df['Gene type']
        _df = pandas.DataFrame(by_prot)
        _df.index = _df['gene_name']
        return _df

    def calculate_pvalues(
        self, which='per_read', log_values=False, 
        log_if_values_above=1E9, apply_bh_adjust=True,
        # This cutof_for_consideration value is important. It has to be
        # 5 for the Figure 1 HOMER motifs to be correct, or CELF1 and Rbfox1
        # have their correct motif shift down the list while an incorrect motif
        # takes the top spot. In general, however, it shouldn't be used.
        cutoff_for_consideration=-5.,
        test_mode=False, slow=False):
        """ For all RNAs and proteins, fit a gaussian KDE to randoms and get a p value.
        Sets self.pvals_kde as {gene_type (RNA) => {protein=>p value}}

        gene_type: "SMAD2::intron" (gene::type)
        counts_by_protein: dict of {gene_type: {protein: [replicate 1, replicate 2]}}
            That is to say a dict by RNA, then by protein (CDK4).
        """
        self.stats_log = collections.defaultdict(int)

        print("dict {gene_type (RNA) => {protein=>count}")
        if which == 'per_read':
            pos_rpg = self.positives.reads_per_million
            neg_rpg = self.negatives.reads_per_million
        elif which == 'raw':
            pos_rpg = self.positives.raw_reads_per_gene
            neg_rpg = self.negatives.raw_reads_per_gene
        elif which == 'per_protein':
            pos_rpg = self.positives.reads_per_protein
            neg_rpg = self.negatives.reads_per_protein
        else:
            raise ValueError(
                f"which={which}. Possible values are 'per_read', 'raw', and 'per_protein'.")
        
        pos_counts_by_protein = pos_rpg.df
        neg_counts_by_protein = neg_rpg.df
        
        missing_rnas_in_negatives = set(pos_counts_by_protein.index) - set(neg_counts_by_protein.index)
        for rna in missing_rnas_in_negatives:
            neg_counts_by_protein.loc[rna] = 0
        #neg_counts_by_protein = neg_counts_by_protein.append(
        #    pandas.Series(0, index=neg_counts_by_protein.columns), ignore_index=True)

        print(f"Determining p values by negative binomial for {pos_counts_by_protein.shape[0]} RNAs...")

        self.pvals[which] = collections.defaultdict(dict)

        self.neg_sub_df = {}
        self.pos_sub_df = {}

        # Create dicts of {protein -> dataframe} by subsetting columns.
        # For the negative proteins:

        for prot in neg_rpg.proteins():
            self.neg_sub_df[prot] = neg_counts_by_protein.loc[:,
                neg_rpg.columns_for_a_protein(prot)]
            self.neg_sub_df[prot].fillna(value=0, inplace=True)
            lower_bound = min([self.negatives.lowest_positive_vals[which][prot]/10, 1])
            self.neg_sub_df[prot].clip(lower_bound, None, inplace=True)
            self.neg_sub_df[prot] = self.neg_sub_df[prot].mean(axis=1)

        # For the positive proteins:
        for prot in pos_rpg.proteins():
            self.pos_sub_df[prot] = pos_counts_by_protein.loc[:,
                pos_rpg.columns_for_a_protein(prot)]
            self.pos_sub_df[prot].fillna(value=0, inplace=True)
            self.pos_sub_df[prot] = self.pos_sub_df[prot].mean(axis=1)
        
        if len(self.neg_sub_df) < 2 or len(self.pos_sub_df) < 1:
            return False  # Abort and return false if there aren't enough proteins to compare.

        all_proteins = neg_rpg.proteins() | pos_rpg.proteins()

        print(f"Columns for hnRNPC: {pos_rpg.columns_for_a_protein('hnRNPC')}")
        if test_mode:
            #pos_counts_by_protein['sum'] = pos_counts_by_protein.sum(axis=1)
            #pos_counts_by_protein.sort_values(by=['sum'], inplace=True, ascending=False)
            pos_counts_by_protein = pos_counts_by_protein.head(20)

        # For every RNA:
        for gene_type in pos_counts_by_protein.index:

            verbose = bool(random.randint(0, 10000) == 1)

            neg_vals = {
                prot: self.neg_sub_df[prot].loc[gene_type] for prot in self.neg_sub_df \
                    if (gene_type.split('::')[0] != prot)}

            pos_counts_by_protein = {
                prot: self.pos_sub_df[prot].loc[gene_type] for prot in self.pos_sub_df}

            #print("Neg vals {}\n pos vals {}\n".format(neg_vals, pos_counts_by_protein))

            verbose and print(f"{gene_type}:")

            verbose and print(f"Calculating p values from neg. binom. fit to randoms: on RNA {len(self.pvals[which])}/{len(neg_counts_by_protein)}.")

            if (np.max([x for x in neg_vals.values()]) < cutoff_for_consideration) and (
                np.max([x for x in pos_counts_by_protein.values()]) < cutoff_for_consideration):
                #print("Setting to 1")
                #self.pvals[which][gene_type] = {prot:1 for prot in all_proteins}
                continue

            neg_vals_at_rna = np.array([x for x in neg_vals.values()])
            #print(neg_vals_at_rna)
            #if slow:
            #    self.pvals[which][gene_type] = self.pval_at_rna_by_nbinom_slow(  
            #        pos_counts_by_protein,
            #        neg_vals_at_rna,
            #        gene_type, log_values=log_values,
            #        log_if_values_above=log_if_values_above, which=which)
            #else:
            self.pvals[which][gene_type] = self.pval_at_rna_by_nbinom(  
                pos_counts_by_protein,
                neg_vals_at_rna,
                gene_type, log_values=log_values,
                log_if_values_above=log_if_values_above, which=which,
                verbose=verbose)

            verbose and print(
                f"pos={pos_counts_by_protein} neg={neg_vals_at_rna}. pvals={self.pvals[which][gene_type]}")

            #print(self.pvals[which][gene_type])

        if not hasattr(self, 'pvals_dfs'):
            self.pvals_dfs = {}

        self.pvals_dfs[which] = pandas.DataFrame.from_dict(self.pvals[which], orient='index')

        if apply_bh_adjust:
            
            fdr_cols = []
            for prot in self.pvals_dfs[which].columns[:]:
                arr = statsmodels.stats.multitest.fdrcorrection(
                    self.pvals_dfs[which][prot], alpha=0.05,
                    method='indep', is_sorted=False)

                self.pvals_dfs[which]['FDR' + prot] = arr[1]
                fdr_cols.append('FDR' + prot)
            for col in fdr_cols:
                self.pvals_dfs[which][re.sub('FDR', '', col)] = self.pvals_dfs[which][col]
                del self.pvals_dfs[which][col]

            print("Warning: P values in pvals_df were BH ajusted, but self.pvals were not.")

        print(self.stats_log)

        return True  # Return True if finished successfully.

    def mask_low_absolute_counts(self, which='per_read', raw_cutoff=5):

        for protein in self.positives.raw_reads_per_gene.proteins():

            cols = self.positives.raw_reads_per_gene.columns_for_a_protein(protein)

            # Total reads per RNA.
            _l = self.positives.raw_reads_per_gene.df[cols].sum(axis=1)

            # RNAs to set their P values to 1 for having too low an absolute read count.
            too_low = [name for (name, v) in zip(_l.index, _l) \
                if (v<raw_cutoff and (name in self.pvals_dfs[which].index))]

            self.pvals_dfs[which].loc[too_low, protein] = 1.
                #is_too_low(gene, sig) for (gene, sig) in zip(_p.index, _p)]
    """
    def pval_at_rna_by_nbinom_slow(
        self, pos_dict_of_counts, neg_vals_at_rna, gene_and_type,
        log_if_values_above=1E9,
        log_values=False, which='per_read'):

        wo_one_protein = {}
        pvals = {}
        for prot in neg_dict_of_counts:
            wo_one_protein = dict(zip(
                neg_dict_of_counts.keys(),
                neg_dict_of_counts.values()))
            del wo_one_protein[prot]

            pval_this_negative = self.pval_at_rna_by_nbinom(
                {prot: neg_dict_of_counts[prot]}, wo_one_protein, gene_and_type,
                log_if_values_above=log_if_values_above,
                log_values=log_values, which=which)
            pvals.update(pval_this_negative)

        only_positives = {prot:val for prot, val in zip(
            pos_dict_of_counts.keys(), pos_dict_of_counts.values()
            ) if (prot not in neg_dict_of_counts)}

        pvals_positives = self.pval_at_rna_by_nbinom(
                only_positives, neg_vals_at_rna, gene_and_type,
                log_if_values_above=log_if_values_above,
                log_values=log_values, which=which)

        pvals.update(pvals_positives)

        return pvals
    """
        
    def pval_at_rna_by_nbinom(
        self, pos_dict_of_counts: Mapping[str, List], neg_vals_at_rna: np.array, gene_and_type,
        log_if_values_above=1E9,
        log_values=False, which='per_read',
        verbose=False):
        """For a given RNA, get the p values for all proteins by negative binomial.
        gene_and_type: "SMAD2::exon"
        dict_of_counts: dict of {protein: [replicate 1, replicate 2]}
        """

        if len(neg_vals_at_rna) == 0:
            return None

        log_scale_high_value = (np.mean(neg_vals_at_rna) > log_if_values_above)

        if log_values or log_scale_high_value:
            log_this_gene = True
            neg_vals_at_rna = np.log10(neg_vals_at_rna)
        else:
            log_this_gene = False
        
        #if not np.any(neg_vals_at_rna):
            #print("No positive values in negatives.")
        #    neg_vals_at_rna = np.array([
        #        self.negatives.lowest_positive_vals[which][x]/10 for x in \
        #            self.negatives.metadata.random_proteins])
            #print(f"negatives now {neg_vals_at_rna}")
        mean_negative = np.average(neg_vals_at_rna)
        std_negative = np.std(neg_vals_at_rna)

        vmr = (std_negative**2)/mean_negative

        verbose and print(f'vmr for negatives={vmr}')
        # Use a poisson if the var/mean is low enough:
        if vmr < 2:
            verbose and print("Using poisson.")
            self.stats_log['vmr<2'] += 1
            pois = stats.poisson(mean_negative)
            return self.use_dist(pos_dict_of_counts, log_this_gene, pois)

        verbose and print("Wil try to use NB.")
        self.stats_log['vmr>=2'] += 1

        # Try to fit a NB useing statsmodels.
        q = sm.NegativeBinomial(
            neg_vals_at_rna, np.array([1] * len(neg_vals_at_rna)), loglike_method='nb2')
        try:
            res = q.fit(disp=0)
        except:  # If a NB can't be fit, revert to a poisson.
            print(f"Could not run q.fit(disp=0) on neg_vals_at_rna= {neg_vals_at_rna}. Using poisson.")
            pois = stats.poisson(mean_negative)
            return self.use_dist(pos_dict_of_counts, log_this_gene, pois)

        # Create a scipy.stats.nbinom object to use its cdf, based on the statsmodels fit parameters.
        # There is no cdf function for the statsmodels object.
        mu = res.predict()[0]  # alpha = res.params[1]
        size = 1. / res.params[1]  # prob = size / (size + mu)

        verbose and print(f"Fit NB mu={mu}")
    
        pvals = self.use_dist(
            pos_dict_of_counts, log_this_gene, stats.nbinom(size, size/(size + mu)))

        return pvals

    def use_dist(self, pos_dict_of_counts, log_this_gene, dist_est):
        _pvals = {}

        for protein, average in pos_dict_of_counts.items():

            if np.isnan(average) or average == 0:
                _pvals[protein] = 1.
                continue
                
            if log_this_gene:
                average = np.log10(average)

            _pvals[protein] = 1 - dist_est.cdf(average)
        
        return _pvals

    def targets(self, protein, cutoff=0.01, which='per_read'):
        pvals = self.pvals_dfs[which]

        targets = pvals[pvals[protein]<cutoff].copy()
        targets = set(targets.index)
        return targets


