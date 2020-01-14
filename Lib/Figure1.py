
import numpy as np
import pandas as pd
import re
from Bio import SeqIO
import argparse
import Figure2

#Python - R interface
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

import matplotlib.pyplot as plt


def define_samples_cutoff(infile_counts, ei_cutoff=None, counts_cutoff=10):
    """
    Loads up raw counts to define the samples to use
    :param infile_counts: infile of the EI counts
    :param ei_cutoff: E or I to use for cutoff
    :param counts_cutoff: minimum counts to use
    :return:
    """

    # Find samples where we have a minimum number of counts
    df_counts = pd.read_csv(infile_counts, sep = ",",index_col=0)
    if ei_cutoff is not None:
        assert ei_cutoff in ["I","E"]
        # Subset to either the Effector or Immunity genes
        total_use = df_counts.loc[:,[x for x in df_counts.columns if re.search("^{}_".format(ei_cutoff),x)]].sum(1)
    else:
        total_use = df_counts.sum(1)
    return df_counts.index[np.where((total_use > counts_cutoff))[0]]


def hash_seq_lengths(infile_seqs):
    """
    Hash the lengths of all the sequences for normalization
    :param infile_seqs: fasta file of the sequences of interest
    :return:
    """
    lengths = {}
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            name = seq.id
            assert name not in lengths, "%s appears twice in your sequence file" %(name)
            lengths[name] = float(len(seq))
    return lengths


def load_data(infile_counts, infile_seqs, infile_total_reads):
    """
    Load up and normalize the data
    :param infile_counts: infile for a dataframe of counts
    :param infile_seqs: sequence fasta file corresponding to that dataframe
    :param infile_total_reads: infile for a dataframe of the total number of reads
    :return:
    """

    # Load up the counts mapping to each gene
    df_counts = pd.read_csv(infile_counts, sep = ",", index_col=0)
    df_counts["name"] = df_counts.index
    # Load up the total library counts
    df_total_reads = pd.read_csv(infile_total_reads, sep=",", index_col=0)
    # Load up the sequence lengths
    lengths = hash_seq_lengths(infile_seqs)

    # Merge on total library size
    df_counts = df_counts.merge(df_total_reads, "inner", on=["name"])

    assert df_counts.shape[0] > 0, "No lines after merging your gene counts and library counts"

    # Normalize the data
    df_norm = pd.DataFrame()
    for j in range(df_counts.shape[1]):
        # Not relevant columns
        if (df_counts.columns[j] == "name") | (df_counts.columns[j] == "reads"):
            pass
        else:
            # Normalized by sequence length and by total library reads
            df_norm[df_counts.columns[j]] = (df_counts.iloc[:, j] / lengths[df_counts.columns[j]]) / df_counts.loc[:, "reads"]
    df_norm.index = df_counts["name"]
    return df_norm



def define_groups(df_scatter):
    """
    Define groups where the immunity gene is likely to be orphan based on abundance level
    :param df_scatter:
    :return:
    """
    orphan_status = []
    for (srs, dat) in df_scatter.iterrows():
        # This is our pseudo count
        if dat["Bfr"] == -14.0:
            orphan_status.append("Absent")
        # If it is 10X more
        elif dat["EI"] - dat["Bfr"] >= 1:
            orphan_status.append("High")
        else:
            orphan_status.append("Normal")

    orphan_status = np.array(orphan_status)
    df_scatter["OrphanStatus"] = orphan_status
    return df_scatter


def figure1_plot(df_norm_ei, df_norm_bfr, outfile, plot_groups = True):
    """
    Plot the scatter plot for Figure1A
    :param df_norm_ei: immunity or effector data
    :param df_norm_bfr: bfr abundance data
    :param outfile: outfile to write the plot
    :param plot_groups: include the groups on the plot
    :return:
    """

    # Define log10 transformed data with pseudo count
    df_plot = pd.DataFrame(
            {"samples": df_norm_ei.index, "EI": np.log10(1e-14 + df_norm_ei.apply(lambda x: np.mean(x[x > 0]), 1)),
         "Bfr": np.log10(1e-14 + df_norm_bfr)})

    if plot_groups:
        df_plot = define_groups(df_plot)
        plot_R = "TRUE"
    else:
        plot_R = "FALSE"

    df_plot.to_csv("temp2.csv",sep=",")

    # Plot the figure
    ro.globalenv['df_plot'] = df_plot
    ro.r('''
        source('Lib/Figure1.R')
        size <- 1.5
        size <- 2.0
        GA_Bfragilis_Scatter(df_plot, "%s", size, size, pseudo = -14.0, plot_groups = %s)
    ''' % (outfile, plot_R))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots the relationship between GA3 gene abundance and the Bacteroides fragilis abundance data')
    parser.add_argument('--infile_counts_ei', help='A dataframe of counts mapping to GA3 EI genes')
    parser.add_argument('--infile_seqs_ei', help='A fasta file of the sequences represented by infile_counts_ei')
    parser.add_argument('--infile_total_reads', help='A dataframe of the total number of reads in each library')
    parser.add_argument('--infile_bacteroides_counts', help='A dataframe of counts mapping to Bacteroides marker genes')
    parser.add_argument('--infile_bacteroides_seqs', help='A fasta file of the sequences represented by infile_bacteroides_counts')
    parser.add_argument('--indir_marker_genes', help='A directory containing markers_s__Bacteroides_XXX.fasta files')
    args = parser.parse_args()

    # Figure 1A
    df_norm = load_data(args.infile_counts_ei, args.infile_seqs_ei, args.infile_total_reads)
    df_norm_E = df_norm.loc[:,[x for x in df_norm.columns if re.search("^E_", x)]]
    df_norm_I = df_norm.loc[:,[x for x in df_norm.columns if re.search("^I_", x)]]
    df_data = Figure2.create_bacteroides_abundance_table(args.infile_bacteroides_counts, args.infile_total_reads, args.infile_bacteroides_seqs, args.indir_marker_genes)
    df_norm_bfr = df_data["markers_s__Bacteroides_fragilis.fasta"]

    # Threshold to samples with sufficient evidence for Immunity genes
    samples_use = define_samples_cutoff(args.infile_counts_ei, ei_cutoff="I", counts_cutoff=10)
    df_norm_I_sub = df_norm_I.loc[samples_use,:]
    df_norm_bfr_sub = df_norm_bfr.loc[samples_use]

    figure1_plot(df_norm_I_sub, df_norm_bfr_sub, "Figure1A.pdf", plot_groups=True)

    #Figure 1B
    samples_use = define_samples_cutoff(args.infile_counts_ei, ei_cutoff="E", counts_cutoff=10)
    df_norm_E_sub = df_norm_E.loc[samples_use, :]
    df_norm_bfr_sub = df_norm_bfr.loc[samples_use]

    figure1_plot(df_norm_E_sub, df_norm_bfr_sub, "Figure1B.pdf", plot_groups=True)
