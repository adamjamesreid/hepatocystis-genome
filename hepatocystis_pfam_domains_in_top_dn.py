#!/usr/bin/env python3
"""
Script for checking for enrichment for PFAM domains in the top 1% of Hepatocystis genes ranked by dN.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License
"""

import hep_project_shared_functions as hpf
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact
import statsmodels.stats.multitest

def load_pfam_domains(pfam_domains_path):
    """
    Loads a tab separated table where column 1 contains Hepatocystis gene names and column 2 contains the PFAM domains in the corresponding genes
    """
    domains_list = list()
    domains_data = hpf.l(pfam_domains_path)
    domains_data = domains_data[1: len(domains_data)]
    pfam_domains_dict = dict()
    for line in domains_data:
        split_line = line.split()
        pfam_domain_entry = split_line[1]
        if pfam_domain_entry != "NA":
            pfam_domains = pfam_domain_entry.split(",")
            pfam_domains_dict[split_line[0]] = pfam_domains
            for pfam_domain in pfam_domains:
                if pfam_domain not in domains_list:
                    domains_list.append(pfam_domain)
    return domains_list, pfam_domains_dict


def load_dn_data(dn_path):
    """
    Loads the dN results table
    """
    dn_data = hpf.l(dn_path)
    dn_data = dn_data[1: len(dn_data)]
    dn_data = [n for n in dn_data if n.split()[1] != "NA"]
    return dn_data


def run_fisher_tests(domains_list, dn_data, pfam_domains_dict, top_gene_cutoff):
    """
    Function for running Fisher exact test to look for enrichment of PFAM domains in the top x% of Hepatocystis genes ranked by dN.
    Input: 1) list of all PFAM domains that occur in Hepatocystis proteins, 2) dN values of Hepatocystis genes, 3) a dictionary that matches Hepatocystis genes with their PFAM domains, 
        4) ranked list cutoff value for sorting Hepatocystis genes into genes with top dN and the rest of the genes.
    Output: 1) list of Fisher exact test p-values for each PFAM domain, uncorrected for multiple hypothesis testing, 2) list of Fisher exact test odds ratios for each PFAM domains
    """
    pvalues_list = list()
    oddsratios_list = list()

    for selected_domain in domains_list:
        topgene_selecteddomain_count = 0
        topgene_not_selecteddomain_count = 0
        not_topgene_selecteddomain_count = 0
        not_topgene_not_selecteddomain_count = 0
        for counter, line in enumerate(dn_data):
            split_line = line.split()
            hep_gene = split_line[0]
            if hep_gene in pfam_domains_dict:
                top_entry = False
                if counter < top_gene_cutoff:
                    top_entry = True
                pfam_domains = pfam_domains_dict[hep_gene]
                for pfam_domain in pfam_domains:
                    if pfam_domain == selected_domain:
                        if top_entry == True:
                            topgene_selecteddomain_count += 1
                        else:
                            not_topgene_selecteddomain_count += 1
                    else:
                        if top_entry == True:
                            topgene_not_selecteddomain_count += 1
                        else:
                            not_topgene_not_selecteddomain_count += 1
        oddsratio, pvalue = fisher_exact([[topgene_selecteddomain_count, topgene_not_selecteddomain_count], [not_topgene_selecteddomain_count, not_topgene_not_selecteddomain_count]])
        pvalues_list.append(pvalue)
        oddsratios_list.append(oddsratio)
    return pvalues_list, oddsratios_list


def main():
    TOP_PERCENTAGE = 1
    pfam_domains_path = "hep_proteins_pfam_domains.txt"
    domains_list, pfam_domains_dict = load_pfam_domains(pfam_domains_path)

    dn_path = "hep_pberghei_povale_3-way_codeml_dn_results.txt"
    dn_data = load_dn_data(dn_path)
    
    top_gene_cutoff = int(round((TOP_PERCENTAGE / 100) * len(dn_data), 0))
    pvalues_list, oddsratios_list = run_fisher_tests(domains_list, dn_data, pfam_domains_dict, top_gene_cutoff)
    qvals = statsmodels.stats.multitest.multipletests(pvalues_list, alpha=0.05, method="fdr_bh")

    print("domain,pvalue,corrected_pvalue,odds_ratio")
    for i in range(0, len(domains_list)):
        if qvals[0][i] == True:
            print(domains_list[i] + "," + str(pvalues_list[i]) + "," + str(qvals[1][i]) + "," + str(oddsratios_list[i]))

if __name__ == "__main__":
    main()



    

