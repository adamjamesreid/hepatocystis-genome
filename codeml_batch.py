#!/usr/bin/env python3
"""
Wrapper script for running codeml (http://envgen.nox.ac.uk/bioinformatics/docs/codeml.html) as batch.
Requires codeml to be installed and in path.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd. 
This program is distributed under the terms of the GNU General Public License
"""

import sys
import os
import hep_project_shared_functions as hpf
import argparse

codeml_ctl_template = """seqfile = {seq_file_path} * sequence data filename
treefile = {tree_file_path}      * tree structure file name
outfile = {out_file_path}           * main result file name
noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
*        ndata = 10
clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
* dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
model = 1
* models for codons:
* 0:one, 1:b, 2:2 or more dN/dS ratios for branches
* models for AAs or codon-translated AAs:
* 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
* 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
* 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
* 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
* 13:3normal>0
icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0
* codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
* AA: 0:rates, 1:separate
fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
kappa = 2  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
omega = .4 * initial or fixed omega, for codons or codon-based AAs
fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0  * different alphas for genes
ncatG = 8  * # of categories in dG of NSsites models
getSE = 1  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .5e-6
cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
method = 0  * Optimization method 0: simultaneous; 1: one branch a time
"""

def write_codeml_ctl_file(seq_file, tree_file, out_file, temp_folder):
    """
    Writes a codeml control file (codeml.ctl) to the temporary files folder
    """
    ctl_file_contents = codeml_ctl_template.format(seq_file_path=seq_file, tree_file_path=tree_file, out_file_path=out_file)
    out_path = temp_folder + "/codeml.ctl"
    with open(out_path, "w") as out_file:
        out_file.write(ctl_file_contents)


def write_temp_alignment_file(alignment_file_path, temp_folder):
    """
    Writes an alignment file to the temporary files folder
    """
    fasta_data = hpf.l(alignment_file_path)
    out_file_path = temp_folder + "/temp_seq.fa"
    with open(out_file_path, "w") as out_file:
        for line in fasta_data:
            if line.startswith(">"):
                line = line.split("_")[0]
            out_file.write(line + "\n")


def write_temp_treefile(treefile_content, temp_folder):
    """
    Writes a Newick tree file to the temporary files folder
    """
    out_file_path = temp_folder + "/temp_tree.treefile"
    with open(out_file_path, "w") as out_file:
        out_file.write(treefile_content[0] + "\n")


def main(alignments_folder, output_folder, temp_folder, treefile_path):
    treefile_content = hpf.l(treefile_path)

    os.system("mkdir -p " + temp_folder)
    os.system("mkdir -p " + output_folder)
    os.chdir(temp_folder)

    alignment_files = os.listdir(alignments_folder)
    for alignment_file in alignment_files:
        alignment_file_path = alignments_folder + "/" + alignment_file
        write_temp_alignment_file(alignment_file_path, temp_folder)
        write_codeml_ctl_file("temp_seq.fa", "temp_tree.treefile", "temp_out.txt", temp_folder)
        write_temp_treefile(treefile_content, temp_folder)
        os.system("codeml codeml.ctl")
        results_file_name = alignment_file.split(".fa")[0] + "_codeml.txt"
        results_file_path = output_folder + "/" + results_file_name
        os.system("cp " + temp_folder + "/temp_out.txt" + " " + results_file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alignments_folder", type=str, help="Path to folder with TranslatorX alignments of transcripts")
    parser.add_argument("output_folder", type=str, help="Path for output folder")
    parser.add_argument("temp_folder", type=str, help="Path for folder with temporary files")
    parser.add_argument("treefile_path", type=str, help="Path to a Newick phylogeny tree of the species that are compared")
    args = parser.parse_args()
    main(args.alignments_folder, args.output_folder, args.temp_folder, args.treefile_path)







