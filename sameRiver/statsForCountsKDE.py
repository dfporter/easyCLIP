class statsForCountsKDE(sameRiver.statsForCounts.statsForCounts):
    
   def all_pvals_from_kde_of_randoms(
        self, which='per_million', log_values=False, bandwidths=None,
        log_if_values_above=1E9, scale_all_numbers_by=None, use_this_counts_by_protein=False):
        """ For all RNAs and proteins, fit a gaussian KDE to randoms and get a p value.
        Sets self.pvals_kde as {gene_type (RNA) => {protein=>p value}}
        """
        
        if use_this_counts_by_protein:
            _counts_by_protein = use_this_counts_by_protein
            print("Using given argument as the dict {gene_type (RNA) => {protein=>count}")
        if (which == 'per_million') or (which == 'per_read'):
            _counts_by_protein = self.per_million_by_protein
        elif which == 'raw':
            _counts_by_protein = self.raw_counts_by_protein
        elif which == 'per_protein':
            _counts_by_protein = self.per_protein_by_protein
        else:
            raise ValueError(
                "which={}. Possible values are 'per_million' (reads), 'raw', and 'per_protein'.".format(which))

        if scale_all_numbers_by is not None:
            counts_by_protein = self.scale_all_numbers(_counts_by_protein, scale_all_numbers_by)
        else:
            counts_by_protein = _counts_by_protein
            
        print("Determining p values by KDE for {} RNAs...".format(
            len(counts_by_protein)))

        self.get_lowest_values_for_all_proteins(counts_by_protein=counts_by_protein)

        pvals = collections.defaultdict(dict)

        for gene_type, dict_of_proteins in counts_by_protein.items():

            if random.randint(0, 1000) == 1:
                print("all_pvals_from_kde_of_randoms() on RNA {}/{}.".format(
                    len(pvals), len(counts_by_protein)))

            pvals[gene_type] = self.p_value_from_kde_of_randoms_at_gene(  
                gene_type, counts_by_protein, log_values=log_values,
                bandwidths=bandwidths, log_if_values_above=log_if_values_above)

        self.pvals_kde = pvals
        return pvals


    def p_value_from_kde_of_randoms_at_gene(
        self, gene_type, counts_by_protein, log_values=False,
        bandwidths=None, log_if_values_above=1E9, scaling_factor=None):

        if gene_type not in counts_by_protein:
            return None

        _dict_of_counts = counts_by_protein[gene_type]

        undetected, detected = (0, [])
        for random_protein in self.list_randoms():

            if gene_type.split('::')[0] == random_protein:
                continue

            replicate_vals_at_protein = counts_by_protein[gene_type].get(random_protein, [np.nan])

            #if np.all(np.isnan(replicate_vals_at_protein)):
            #    undetected += 1
            #    continue

            cleaned = np.nan_to_num(replicate_vals_at_protein)
            cleaned = np.clip(cleaned, self.lowest_positive_value_for_a_protein(
                random_protein, counts_by_protein=counts_by_protein), None)

            #cleaned = cleaned[cleaned>0]

            detected.append(np.mean(cleaned))

        if len(detected) == 0:
            return None

        detected = np.array(detected)
        #detected = np.clip(detected, 0.00001, None)

        if bandwidths is None:
            bandwidth_guess = len(detected)**(-0.2)
            bandwidths = [bandwidth_guess * 1, bandwidth_guess * 2]

        _x = {}
        _cdf = {}

        log_scale_high_value = (np.mean(detected) > log_if_values_above)

        if log_values or log_scale_high_value:
            print(log_scale_high_value)
            print(np.mean(detected))
            print(log_if_values_above)
            log_this_gene = True
            detected = np.log10(detected)
        else:
            log_this_gene = False

        if scaling_factor:
            detected = scaling_factor * detected
            
        for bandwidth in bandwidths:
            _x[bandwidth], _cdf[bandwidth] = self.cdf(
                KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(detected[:, np.newaxis]),
                cdf_min=-10, cdf_max=max(detected)*2)

        fraction_detected = len(detected)/(len(detected) + undetected)

        _pvals = {}

        for protein, replicate_vals in counts_by_protein[gene_type].items():

            #if protein != 'TPGS2':
            #    continue

            #print(detected)
            #print(_cdf)
            
            if np.any(np.isnan(replicate_vals)) or (np.mean(replicate_vals) == 0):
                _pvals[protein] = 1.
                continue
            
            if scaling_factor:
                replicate_vals = scaling_factor * replicate_vals
                
            if log_this_gene:
                average = np.log10(np.mean(replicate_vals))
            else:
                average = np.mean(replicate_vals)

            #print('>>', average)
            _possible_pvals = []
            for bandwidth in _x:
                _index = bisect.bisect(_x[bandwidth], average)
                _possible_pvals.append(
                    (1 - _cdf[bandwidth][_index-1]) * fraction_detected
                )

            if log_this_gene and random.randint(0, 1E5) == 1:
            #if 1:
                print('(Random example) *** {} {}'.format(gene_type, protein))
                print('log scale ? {}'.format(log_this_gene))
                print('Bandwidths: ', _cdf.keys())
                print('P vals', _possible_pvals)
                print('Randos: ', detected)
                print('This protein: ', replicate_vals)

            _pvals[protein] = max(np.nan_to_num(_possible_pvals))

            if _pvals[protein] < 1E-4 and random.randint(0, 1E3) == 1:
            #if 1:
                print('(Random example) *** {} {}'.format(gene_type, protein))
                print('log scale ? {}'.format(log_this_gene))
                print('Bandwidths: ', _cdf.keys())
                print('P vals', _possible_pvals)
                print('Prot mean {} randos mean {}'.format(average, np.mean(detected)))

        return _pvals

            @staticmethod
    def cdf(kde, cdf_min=-7, cdf_max=50):

        _x = np.linspace(cdf_min, cdf_max, 1000)

        Y = [10**x for x in kde.score_samples(_x[:, np.newaxis])] 

        cumulative = [Y[0]]    
        for val in Y[1:]:
            cumulative.append(cumulative[-1] + val)

        total = cumulative[-1]
        cumulative = [x/total for x in cumulative]

        return _x, cumulative
    