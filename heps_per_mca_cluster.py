# Given orthoMCL data for Pb and Hep, count orthologue groups shared by both and unique to Pb

import sys
import re
import scipy, scipy.stats
import statsmodels.stats.multitest

ortho_file = sys.argv[1]
mca_cluster_file = sys.argv[2]
gene_desc_file = sys.argv[3]

species_of_interest = 'Hepatocystis_DNA'

# FDR cutoff for clusters with unexpectedly low numbers of Hepatocystis genes
fdr_cutoff = 0.05

# For a given Pb id, we want to know if it shares an orthologous group with Hep


class Orthomcl_result:
  ''' parse orthomcl results and provide access to them '''
  def __init__ (self, orthomcl_file):
    self.orthomcl_file = orthomcl_file
    self.ngenes = dict()
    self.ntaxa = dict()
    self.taxa = dict()
    self.genes = dict()
    self.groups = list()
    self.group_for_gene = dict()
    self.genes_for_group = dict()
    self.taxa_for_group = dict()
    self.parse()

  def parse (self):
    with open(self.orthomcl_file) as f:
      for x in f:
        x = x.rstrip()
        head, members = x.split('\t ')
        name, numbers = head.split('(')
        self.groups.append(name)
        m = re.findall(r'(\d+) genes', numbers)
        ngenes = m[0]
        m = re.findall(r'(\d+) taxa', numbers)
        ntaxa = m[0]

        self.ngenes[name] = ngenes
        self.ntaxa[name] = ntaxa
        
        mem_list = members.split(' ')
        for y in mem_list:
          gene_id, species_plus = y.split('(')
          species, _ = species_plus.split(')')
          if species not in self.taxa:
            self.taxa[species] = list()
          self.taxa[species].append(gene_id)
          if gene_id not in self.genes:
            self.genes[gene_id] = list()
          self.genes[gene_id] = species
          self.group_for_gene[gene_id] = name
          if name not in self.genes_for_group:
            self.genes_for_group[name] = list()
          self.genes_for_group[name].append(gene_id)
          if name not in self.taxa_for_group:
            self.taxa_for_group[name] = list()
          self.taxa_for_group[name].append(species)

  def shares_with (self, gid, taxaid):
    ''' Does gene gid share an ortho group with a gene from taxa taxaid '''
    group = self.group_for_gene[gid]
    #print(group, gid)
    if group in self.taxa_for_group and taxaid in self.taxa_for_group[group]:
      #print(f'return True for', group, gid, taxaid)
      return True
    else:
      #print(f'return False for', group, gid, taxaid)
      return False
      
  def gene_in_group (self, gid):
    ''' Does gene gid exist in a group '''
    if gid in self.group_for_gene:
      #return self.group_for_gene
      return True
    else:
      return False

  def taxon_in_group (self, taxon, group):
    ''' Does a taxon exist in a group? '''
    if taxon in self.taxa_for_group[group]:
      return True
    else:
      return False

  def total_shared_between_taxa (self, taxon_a, taxon_b):
    ''' For a pair of taxa, determine how many ortho groups are shared between them '''
    count = 0
    #Loop through groups and ask if taxon_a and taxon_b are in there and keep count
    for x in self.groups:
      if self.taxon_in_group(taxon_a, x) and self.taxon_in_group(taxon_b, x):
        count = count + 1

    return count

  def total_shared_between_taxa_excluding (self, taxon_a, taxon_exc, *taxa_b):
    ''' number of groups shared between taxon_a and taxa_b (list), but not taxon_exc '''
    count = 0
    # Loop through groups and ask if taxon_a and taxa_b are present, but not taxon_exc
    for x in self.groups:
      if self.taxon_in_group(taxon_a, x) and not self.taxon_in_group(taxon_exc, x):
        exc_results = [self.taxon_in_group(b, x) for b in taxa_b]
        if exc_results.count(True) == len(exc_results):
          count = count + 1
    return count

# Instantiate Orthomcl_result class
orh = Orthomcl_result(ortho_file)

# Get Pb gene descriptions
gene_desc = dict()

with open(gene_desc_file) as gf:
  for x in gf:
    x = x.rstrip()
    gid, desc = x.split('\t')
    gene_desc[gid] = desc

clusters = dict()

with open(mca_cluster_file) as f:
  for x in f:
    x = x.rstrip()
    if x.startswith('feature_symbol'):
      continue
    gene_id, cluster = x.split('\t')
    if cluster not in clusters:
      clusters[cluster] = list()
    clusters[cluster].append(gene_id)

# Get totals for shared ortho groups (for chisq calculations)
# a: between Pb and Hep
# b: between Pb and Po/Pv but not Hep
groups_shared_Pb_Hep = orh.total_shared_between_taxa("PbergheiANKA", species_of_interest)
groups_shared_Pb_Po_Pv_not_Hep = orh.total_shared_between_taxa_excluding("PbergheiANKA", species_of_interest, "Povalewallikeri", "PvivaxP01")
#print(groups_shared_Pb_Hep, groups_shared_Pb_Po_Pv_not_Hep)

pvals = list()
results = list()

for c in clusters:
  genes_in_cluster = len(clusters[c])
  genes_in_groups = 0
  shared_with_hep = 0
  shared_with_Pow = 0
  shared_with_Pv = 0
  shared_with_ovale_vivax_not_hep = 0
  for g in clusters[c]:
    if orh.gene_in_group(g):
      genes_in_groups = genes_in_groups + 1
      if orh.shares_with(g, species_of_interest):
        shared_with_hep = shared_with_hep + 1
      if orh.shares_with(g, 'Povalewallikeri'):
        shared_with_Pow = shared_with_Pow + 1
      if orh.shares_with(g, 'PvivaxP01'):
        shared_with_Pv = shared_with_Pv + 1
      if orh.shares_with(g, 'Povalewallikeri') and orh.shares_with(g, 'PvivaxP01') and not orh.shares_with(g, species_of_interest):
        shared_with_ovale_vivax_not_hep = shared_with_ovale_vivax_not_hep + 1
    
  hep_share = shared_with_hep * 100 / genes_in_cluster
  ovale_share = shared_with_Pow * 100 / genes_in_cluster
  vivax_share = shared_with_Pv * 100 / genes_in_cluster

  # Calculate difference between ovale/vivax average and hep
  midpoint_diff = hep_share - ((ovale_share + vivax_share) / 2)

  #print(c, genes_in_groups, genes_in_cluster, hep_share, ovale_share, vivax_share, midpoint_diff, sep='\t')

  # Calculate chi-sq/Fisher's p-value
  expected_shared_with_hep = (shared_with_hep + shared_with_ovale_vivax_not_hep) * (groups_shared_Pb_Hep / (groups_shared_Pb_Hep + groups_shared_Pb_Po_Pv_not_Hep))
  expected_shared_with_ovale_vivax_not_hep = (shared_with_hep + shared_with_ovale_vivax_not_hep) * (groups_shared_Pb_Po_Pv_not_Hep / (groups_shared_Pb_Hep + groups_shared_Pb_Po_Pv_not_Hep))

  chisq, p = (None, None)
  fisher_oddsratio, fisher_p = (None, None)
  if(shared_with_hep > 1 and shared_with_ovale_vivax_not_hep > 1):
    chisq, chisq_p = scipy.stats.chisquare([shared_with_hep, shared_with_ovale_vivax_not_hep], f_exp=[expected_shared_with_hep, expected_shared_with_ovale_vivax_not_hep])
    #print([[shared_with_hep, shared_with_ovale_vivax_not_hep], [expected_shared_with_hep, expected_shared_with_ovale_vivax_not_hep]])
    fisher_oddsratio, fisher_p = scipy.stats.fisher_exact([[shared_with_hep, shared_with_ovale_vivax_not_hep], [expected_shared_with_hep, expected_shared_with_ovale_vivax_not_hep]])
    pvals.append(fisher_p)
    results.append([c, genes_in_groups, genes_in_cluster, shared_with_hep, shared_with_ovale_vivax_not_hep, expected_shared_with_hep, expected_shared_with_ovale_vivax_not_hep, chisq, chisq_p, fisher_oddsratio, fisher_p])

qvals = statsmodels.stats.multitest.multipletests(pvals, alpha = 0.05, method='fdr_bh')
#print(pvals, qvals)
for x in range(0, len(results)):
  res_str = '\t'.join([str(s) for s in results[x]])
  fdr = qvals[1][x]
  if (fdr <= fdr_cutoff):
    print(res_str, fdr, sep='\t')
    #print(clusters[results[x][0]])
    # for a particular cluster, determine genes which are in a group with Pow and Pv, but not Hep
    for g in clusters[results[x][0]]:
      #print(g)
      if orh.gene_in_group(g):
        if orh.shares_with(g, 'Povalewallikeri') and orh.shares_with(g, 'PvivaxP01') and not orh.shares_with(g, species_of_interest):
          #print(g)
          if g not in gene_desc:
            gene_desc[g] = ''
          print(results[x][0], g, gene_desc[g], sep='\t')