
import matplotlib
matplotlib.use('agg')

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import re
import os
import tempfile
import subprocess
from collections import defaultdict

import Figure1
import numpy as np
import pandas as pd
import argparse

#Python - R interface
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

import matplotlib.pyplot as plt


def get_genome_sequences(gene_oi, len_oi, infile_seqs, infile_db):
    """
    Runs a BLAST search to get the genes from infile_seqs
    :param gene_oi: gene that we want the sequence of
    :param len_oi: Length of the gene of interest
    :param infile_seqs: fasta file containing gene_oi
    :param infile_db: BLAST database of the genomes
    :return:
    """
    Df = do_BLAST(infile_seqs, infile_db)
    Df["Query"] = [x.split(" ")[0] for x in Df["Query"]]
    Df = Df.query('Query == "%s"' %(gene_oi) )
    Df = Df.query('(query_len == %i.0) & (align_len == %i.0)' %(len_oi, len_oi))

    assert Df.shape[0] > 0
    genes_oi = defaultdict(list)

    for (i, dat) in Df.iterrows():
        genes_oi[dat["Query"]].append(dat["Hit"])

    seqs = SubsetFastaWithList("/home/averster/nexus/vol1/T6S/Data/Strains/GCA_AllStrains_mRNA_Curated_Alt_StopFix.fasta", [x.split(" ")[0] for x in genes_oi[gene_oi]])
    return seqs


def SubsetFastaWithList(infile, list_oi):
    """
    Loads up a fasta infile and subsets the sequences to using list_oi. Only sequences whose id appears in list_oi
    :param infile:
    :param list_oi:
    :return:
    """
    f = ReadFasta(infile)
    items = []
    for record in f:
        if record.id in list_oi:
            items.append(record)
    return items


def ReadFasta(infile):
    '''
    Reads a gziped genome fasta file downloaded from NCBI
    Takes the input from NCBI_GetGenomeFromTaxID()
    '''
    f = open_maybe_gzip(infile)
    for record in SeqIO.parse(f, "fasta"):
        yield record
    f.close()


def open_maybe_gzip(infile):
    #Check to see if this is gzipped
    if infile[-3:] == ".gz":
        f = gzip.open(infile,'rb')
    else:
        f = open(infile, 'rb')
    return f


def do_BLAST(infile_seqs, infile_db, E_VALUE_THRESH = 10e-10, FRAC_ID_THRESH = 0.95, FRAC_LEN_THRESH = 0.50, blast_prog = "blastn", return_table = True):
    """
    Runs BLAST of a fasta file against a database
    :param infile_seqs:
    :param infile_db:
    :param E_VALUE_THRESH:
    :param FRAC_ID_THRESH:
    :param FRAC_LEN_THRESH:
    :param blast_prog:
    :param return_table:
    :return:
    """
    #If seqs is a seq object, write it to a tempfile
    if isinstance(infile_seqs,list):
        if all([isinstance(x,Seq) for x in infile_seqs]) | all([isinstance(x,SeqRecord) for x in infile_seqs]):
            seq_file_use = tempfile.NamedTemporaryFile().name
            with open(seq_file_use,"w") as f:
                SeqIO.write(infile_seqs, f, "fasta")
    elif isinstance(infile_seqs,Seq) | isinstance(infile_seqs,SeqRecord):
        seq_file_use = tempfile.NamedTemporaryFile().name
        with open(seq_file_use,"w") as f:
            SeqIO.write(infile_seqs, f, "fasta")
    else:
        seq_file_use = infile_seqs

    #Create a blastdb if we have a sequence file
    if isinstance(infile_db,list):
        if all([isinstance(x,Seq) for x in infile_db]) | all([isinstance(x,SeqRecord) for x in infile_db]):
            db_file_use = tempfile.NamedTemporaryFile().name
            infile_db_use = db_file_use + "BLASTDB"
            with open(db_file_use,"w") as f:
                SeqIO.write(infile_db, f, "fasta")
            if blast_prog == "blastn":
                t = "nucl"
            elif blast_prog == "blastp":
                t = "prot"
            else:
                sys.exit("Can't decide on nucl or prot from your blast_prog")
            cmd = ["makeblastdb", "-in", db_file_use, "-dbtype", t, "-out", infile_db_use]
            subprocess.call(cmd)
    elif isinstance(infile_db,Seq) | isinstance(infile_db,SeqRecord):
        db_file_use = tempfile.NamedTemporaryFile().name
        infile_db_use = db_file_use + "BLASTDB"
        with open(db_file_use,"w") as f:
            SeqIO.write(infile_db, f, "fasta")
        if blast_prog == "blastn":
            t = "nucl"
        elif blast_prog == "blastp":
            t = "prot"
        else:
            sys.exit("Can't decide on nucl or prot from your blast_prog")
        cmd = ["makeblastdb", "-in", db_file_use, "-dbtype", t, "-out",infile_db_use]
        subprocess.call(cmd)
    else:
        infile_db_use = infile_db
    tf_out = tempfile.NamedTemporaryFile()
    assert os.path.isfile(infile_db_use + ".nhr") | os.path.isfile(infile_db_use + ".pal") | os.path.isfile(infile_db_use + ".phr")
    cmd = [blast_prog, "-query", seq_file_use, "-db", infile_db_use, "-out", tf_out.name, "-outfmt", "5"]
    subprocess.call(cmd)

    r = parse_BLAST(tf_out.name, E_VALUE_THRESH = E_VALUE_THRESH, FRAC_ID_THRESH = FRAC_ID_THRESH, FRAC_LEN_THRESH = FRAC_LEN_THRESH, verbose = False)
    if return_table:
        Df = make_summary_table(r)
        return Df
    else:
        return r


def make_summary_table(blast_hits):
    """
    Makes a summary table from the BLAST results from parse_BLAST()
    :param blast_hits:
    :return:
    """
    hits = {}
    #Start by finding all the hits
    for query in blast_hits:
        for hit in blast_hits[query]:
            if hit["hit"] not in hits:
                hits[hit["hit"]] = []
            hits[hit["hit"]].append(query + "[Eval]" + str(hit["eval"]) + "[ID]" + "%f" %(hit["id"]))
    #Now make a table
    i = 0
    Df = pd.DataFrame(columns = ["Query","Hit","AltQueries","ID","Eval","query_len", "align_len", "sbjct_start"])
    for query in blast_hits:
        for hit in blast_hits[query]:
            #Remove alt queries which are the current query
            alt_queries = []
            for query_alt in hits[hit["hit"]]:
                if not query in query_alt:
                    alt_queries.append(query_alt)
            #Now sort
            if len(alt_queries) > 0:
                alt_queries = sorted(alt_queries, key = lambda x: float(re.search("\[ID\](.*)",x).group(1)))
            Df.loc[i] = [query,hit["hit"],"_".join(alt_queries),hit["id"], hit["eval"], hit["query_len"], hit["align_len"], hit["sbjct_start"]]
            i += 1
    return Df


def parse_BLAST(infile, E_VALUE_THRESH = 1e-30, FRAC_ID_THRESH = None, FRAC_LEN_THRESH = None, verbose = False, mode = "hit_def"):
    """
    Parses a BLAST xml file into a dictionary with associated cutoffs
    :param infile:
    :param E_VALUE_THRESH:
    :param FRAC_ID_THRESH:
    :param FRAC_LEN_THRESH:
    :param verbose:
    :param mode:
    :return:
    """
    r = {}
    with open(infile) as result_handle:
        blast_record = NCBIXML.parse(result_handle)
        for record in blast_record:
            query_name = record.query
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        #frac_id = float(hsp.identities) / float(min(record.query_length, alignment.length))
                        frac_id = float(hsp.identities) / float(hsp.align_length)
                        if FRAC_ID_THRESH is not None:
                            if frac_id < FRAC_ID_THRESH:
                                continue
                        if FRAC_LEN_THRESH is not None:
                            if float(hsp.align_length) / float(record.query_length) < FRAC_LEN_THRESH:
                                continue
                        if verbose:
                            print("QUERY:", query_name)
                            print("HIT:", alignment.title)
                            #Need to divide by 3.0 since it thinks its the DNA sequence
                            print("ID:%f, EVAL:%f, QUERY_LEN: %f, ALIGN_LEN: %f" %(frac_id, hsp.expect, record.query_length, alignment.length))
                            print("SEQ_QUERY: " + hsp.query[0:75] + '...')
                            print("SEQ_MATCH: " + hsp.match[0:75] + '...')
                            print("SEQ_HIT:   " + hsp.sbjct[0:75] + '...')
                            print("")
                        if query_name not in r:
                            r[query_name] = []
                        #old version had "hit": alignment.title which generates some extra test
                        if mode == "hit_def":
                            hit_use = alignment.hit_def
                        elif mode == "title":
                            hit_use = alignment.title
                        r[query_name].append({"hit": hit_use, "eval": hsp.expect, "id": frac_id, "query_len": record.query_length, "align_len": alignment.length, "hsp_len": hsp.align_length, "query_start": hsp.query_start, "sbjct_start": hsp.sbjct_start })
    #Finally sort the results
    for query in r:
        r[query] = sorted(r[query], key = lambda x: x["id"], reverse = True)
    return r


def get_seqs(infile_pileup_list, gene_oi, infile_seqs, use_minor=False, minor_cutoff=5, min_covg = 10, cutoff = 0.90, getref=False):
    """
    Given a list of pileups, this reconstructs the sequence of a gene of interest
    :param infile_pileup_list:
    :param gene_oi:
    :param use_minor:
    :param minor_cutoff:
    :param min_covg:
    :param cutoff:
    :param getref:
    :param infile_seqs:
    :return:
    """

    # Hash the seqs
    all_seqs = []
    seq_dict = {}
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_dict[seq.id] = seq

    for infile in infile_pileup_list:

        bases = ["A", "C", "G", "T"]

        if not os.path.isfile(infile_pileup):
            print("Can't find %s" % (infile_pileup))
            continue

        Counts = get_table(infile, infile_seqs=infile_seqs, bases=bases, return_dict=True)
        array_matrix_oi = np.array(Counts[gene_oi])
        seq = ""
        for i in xrange(array_matrix_oi.shape[0]):

            # Do something if we have 0,0,0,0
            if np.sum(array_matrix_oi[i, :]) == 0:
                base_use = "N"
            base_use = bases[np.argmax(array_matrix_oi[i, :])]
            # If we want the minor allele
            if use_minor:
                minor_counts = sorted(array_matrix_oi[i, :], reverse=True)[1]
                if minor_counts > minor_cutoff:
                    base_use = bases[np.argsort(array_matrix_oi[i, :])[::-1][1]]
            major_counts = sorted(array_matrix_oi[i, :], reverse=True)[0]
            minor_counts = sorted(array_matrix_oi[i, :], reverse=True)[1]
            # If we don't have enough data to resolve this clearly
            if (major_counts - minor_counts) < min_covg:
                base_use = "N"
            seq = seq + base_use

        seqobj = SeqRecord(Seq(seq), id=srs, description="")

        # Only include sequences where we can resolve enough of the sequence
        if (1.0 - seqobj.count("N") / float(len(seqobj))) > cutoff:
            all_seqs.append(seqobj)

    if getref:
        with open(infile_seqs) as f:
            for seq in SeqIO.parse(f, "fasta"):
                if seq.id == gene_oi:
                    all_seqs.append(seq)
    return all_seqs



def get_table(infile_allele, infile_seqs,
               bases_trim = 50, bases = ["A","C","G","T"], return_dict = False):
    """
    Creates a table of allele counts from the pileup file
    :param infile_allele: A pileup file
    :param infile_seqs: reference fasta file of the sequences in infile_allele
    :param bases_trim: NUmber of bases to trim from the end
    :param bases:
    :param return_dict: return a dictionary instead of dataframe
    :return:
    """

    #First get the sequence lengths
    lengths = {}
    counts = {}
    seq_order = []
    with open(infile_seqs) as f:
        # Initialize a matrix of zeros the length of the sequence
        for seq in SeqIO.parse(f, "fasta"):
            lengths[seq.id] = len(seq)
            counts[seq.id] = [[0,0,0,0] for i in xrange(lengths[seq.id])]
            seq_order.append(seq)

    seq_i = defaultdict(int)
    with open(infile_allele) as f:
        for line in csv.reader(f, delimiter = "\t"):
            seq_id = line[0]
            pos = int(line[1]) - 1 #Need to subtract by 1 to make this consistent with the 0-start sequence
            seq = line[4]
            quality = line[5]
            depth = int(line[3])
            reference = line[2]
            if (pos >= lengths[seq_id]):
                continue
            if (pos < bases_trim) | (pos > (lengths[seq_id] - bases_trim)):
                c = [0,0,0,0]
                if reference in bases:
                    c[bases.index(reference)] = depth
                else:
                    c = [depth,0,0,0]
                counts[seq_id][pos] = c
            else:
                parsed = parse_pileup_seq(seq, quality, depth, reference)
                counts[seq_id][pos] = [parsed[bases[0]],parsed[bases[1]],parsed[bases[2]],parsed[bases[3]]]
                seq_i[seq_id] += 1
    if return_dict:
        return counts
    Counts_Data = np.vstack([counts[x.id] for x in seq_order])

    return Counts_Data, seq_order



def get_gene_length(gene_id, infile_seqs):
    """
    Calculates the length of a gene of interest
    :param gene_id: Name of the gene
    :param infile_seqs: Sequence file
    :return:
    """
    len_oi = None
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f,"fasta"):
            if seq.id == gene_id:
                len_oi = len(seq)
    return len_oi


def create_tree(infile_pileup_list, infile_seqs_ei, infile_db, infile_genome_names, cutoff_use = 20, gene_oi = "I_6_9343__bti1"):
    """
    Creates the dendrogram from Figure2A and defines clusters
    :return:
    """

    # Get the reconstructed sequences from the metagenomes
    seqs_metagenomes = get_seqs(infile_pileup_list, gene_oi, infile_seqs_ei)
    
    len_oi = get_gene_length(gene_id, infile_seqs)
    #Get the sequences from the genomes
    seqs_genomes = get_genome_sequences(gene_oi, len_oi, infile_seqs, infile_db)

    #Will Plot the cluster gram

    # There are some near identical Bfr strains that are being removed to make the tree more legible
    seqs_remove = ["B. fr 321 BFRA",
    "B. fr str. I345",
    "B. fr A7 (UDC12-2)",
    "B. fr str. 1009-4-F #7",
    "B. fr str. 1009-4-F #10",
    "B. fr 320_BFRA",
    "B. fr 885_BFRA",
    "B. fr J-143-4",
    "B. fr S24L34",
    "B. fr CL03T12C07",
    "B. fr J38-1",
    "B. fr CL03T00C08",
    "B. fr str. B1 (UDC16-1)",
    "B. fr 894_BFRA",
    "B. fr S24L26"]

    tree_samples = plot_tree(seqs_genomes + seqs_metagenomes, infile_genome_names, cutoff_use, gene_oi, seqs_remove)
    return tree_samples


def rename_genomes(DfResults, infile):
    genomes = [re.sub("_[0-9]*$", "", x.split("|")[-1]) for x in DfResults.index if "ti|" in x]

    DfGenomes = pd.read_csv(infile, sep=",", index_col=False)
    DfGenomes["GenBank FTP"] = [os.path.basename(x) for x in DfGenomes["GenBank FTP"].values]
    DfGenomes.columns = [x.replace(" ", "_") for x in DfGenomes.columns]

    name_dict = {}
    for g in genomes:
        DfGenomesSub = DfGenomes.query('GenBank_FTP == "%s"' % (g))

        organism = str(DfGenomesSub.iloc[0, 0]).replace("Bacteroides fragilis", "B. fr")
        strain = str(DfGenomesSub.iloc[0, 1])
        if (strain in organism) | (strain == "nan"):
            name = organism
        else:
            name = organism + " " + strain
        name_dict[g] = name
    index = DfResults.index.tolist()
    for i in xrange(len(index)):
        ind = index[i]
        if "ti|" in ind:
            index[i] = name_dict[re.sub("_[0-9]*$", "", ind.split("|")[-1])]
    DfResults.index = index
    DfResults.columns = index

    return DfResults


def number_diff(s1,s2):
    """
    Calculate the number of differences between 2 sequences, ignoring N
    :param s1:
    :param s2:
    :return:
    """
    n_Ns = 0
    n_diff = 0
    assert len(s1) == len(s2)
    for i in xrange(len(s1)):
        if (s1[i] == "N") | (s2[i] == "N"):
            n_Ns += 1
        elif (s1[i] != s2[i]):
            n_diff += 1
    return n_diff, n_Ns


def plot_tree(seqs_full, infile_genome_names, h_use=0, seqs_ignore=[]):
    """
    Plots a clustergram of sequences
    :param seqs_full: a list of sequences
    :param h_use: position on the tree to cut off
    :param seqs_ignore: seqs to not include in the graph
    :return:
    """
    results = defaultdict(dict)

    for (s1, s2) in itertools.combinations(seqs_full, 2):
        r = number_diff(s1.seq, s2.seq)
        results[s1.id][s2.id] = r[0]
        results[s2.id][s1.id] = r[0]

    df_results = pd.DataFrame(results)
    df_results.fillna(0, inplace=True)
    df_results = rename_genomes(df_results, infile_genome_names)
    df_results = df_results.ix[
        [x for x in df_results.index if x not in seqs_ignore], [x for x in df_results.columns if x not in seqs_ignore]]
    
    robjects.globalenv['DfResults'] = df_results
    robjects.r('''
        library('extrafont')
        DfResults.Dist <- as.dist(DfResults)
        clusters <- hclust(DfResults.Dist)

        pdf("Figure2_ClusterGram.pdf", width = 6.0, height = 4.5, family = "Arial")
        plot(clusters,  hang = -1, main = "", ylab = "N bases different", horiz = TRUE, cex = 0.63, sub = "", xlab = "")
        dev.off()

        clusters <- cutree(clusters, h = %i)
        cluster_names <- names(clusters)
    ''' % (h_use))

    cluster_names = robjects.globalenv['cluster_names']
    clusters = robjects.globalenv['clusters']
    clusters = np.array(clusters)
    cluster_names = np.array(cluster_names)

    tree_groups = []
    for i in set(clusters):
        samples = [x for x in cluster_names[clusters == i] if ("MH" in x) | ("SRS" in x) | ("V1." in x)]
        # Only accept groups with at least 3 members
        if len(samples) >= 3:
            tree_groups.append(samples)

    return tree_groups


def load_markers(indir, reduced_markers=False):
    """
    Takes Metaphlan2 markers from the bacteroides species
    Marker genes are identified from the markers_info.txt.bz2 file associated with metaphlan2
    Needs to match the genes in your counts file
    :param indir: A directory with marker genes fasta files (named as markers_s__Bacteroides_fragilis.fasta)
    :return:
    """
    infiles = [x for x in os.listdir(indir) if re.search("markers_s__Bacteroides_[a-z]+\.fasta",x)]
    markers_dict = defaultdict(list)
    markers_oi = []
    for inf in infiles:
        with open(indir + inf) as f:
            for seq in SeqIO.parse(f, "fasta"):
                if reduced_markers:
                    infsub = inf.replace("markers_s__","").replace(".fasta","")
                    if infsub not in all_seqs:
                        continue
                    if seq.id not in [x.id for x in all_seqs[infsub]]:
                        continue
                markers_dict[inf].append(seq.id)
                markers_oi.append(seq.id)
    return markers_dict, markers_oi



def create_bacteroides_abundance_table(infile_bacteroides_counts, infile_reads,
                                       infile_bacteroides_seqs, indir_marker_genes):
    """
    Creates an abundance table with data for each Bacteroides species and with the immunity gene
    :param infile_bacteroides_counts:
    :param infile_ei_counts:
    :param infile_reads:
    :param infile_bacteroides_seqs:
    :param infile_ei_seqs:
    :param orphan_gene:
    :return:
    """

    # identify which marker genes correspond to which species
    markers_dict, markers_oi = load_markers(indir_marker_genes)

    df_bacteroides = Figure1.load_data(infile_bacteroides_counts, infile_bacteroides_seqs, infile_reads)

    # Average the marker genes for each species
    df_species = pd.DataFrame()
    for sp in markers_dict:
        df_species[sp] = df_bacteroides.loc[:, markers_dict[sp]].mean(1)

    return df_species


def add_ei_data(df_species, infile_ei_counts, infile_ei_seqs, orphan_gene="I_6_9343__bti1"):
    """
    Adds on abundance data for the orphan immunity
    :param df_species:
    :param infile_ei_counts:
    :param infile_ei_seqs:
    :param orphan_gene:
    :return:
    """
    # Add on the orphanI
    df_ei = Figure1.load_data(infile_ei_counts, infile_ei_seqs, infile_reads)
    df_species["OrphanI"] = df_ei.loc[:, orphan_gene]
    return df_species


def single_species_model(Df, orphanI_col = "OrphanI", transform = ""):
    """
    Calculates the mean square error that assigns the orphan immunity gene to each species
    Because the species abundance is the average of marker genes, the orphan gene abundance should be equal to the encoding species
    :param Df:
    :param orphanI_col:
    :param transform:
    :return:
    """
    assert orphanI_col in Df.columns
    pseudo = 1e-15
    results = []
    # Transform the data is requested
    if transform == "log10":
        for col in Df.columns:
            Df.loc[:,col] = np.log10(Df.loc[:,col] + pseudo)
    # Go through each species
    for col in Df.columns:
        if col == orphanI_col:
            continue
        # Calculate the MSE assuming equal abundance to the species
        diff = Df.loc[:,col] - Df.loc[:,orphanI_col]
        mse = np.sum(diff**2) / float(Df.shape[0])
        relative_error = np.mean(np.absolute(diff) / Df.loc[:,orphanI_col] * 100)
        results.append([col, mse, relative_error])
    # Put the MSE into a DataFrame, sort it and return the data
    Df_R = pd.DataFrame(results)
    Df_R.columns = ["Species","MSE","Relative_error"]
    Df_R["Species"] = [x.replace("markers_s__Bacteroides_","").replace(".fasta","").replace("Bacteroides_","") for x in Df_R["Species"].values]
    Df_R.sort("MSE", inplace = True)
    return Df_R



def assign_species(df_bacteroides_data, clusters, gene_oi="I_6"):
    """
    This creates the graphs for Figure 2B-E
    Uses a least squares model to assign species based on Bacteroides abundance data
    :param df_bacteroides_data: from create_bacteroides_abundance_table()
    :param clusters: from plot_tree()
    :param gene_oi: Which immunity gene to assign
    :return:
    """
    pseudo = 1e-15
    color_dict = {
        "thetaiotaomicron": "#a6d854ff",
        "ovatus": "#8da0cbff",
        "fragilis": "#ffd92fff",
        "oleiciplenus": "blue",
        "propionicifaciens": "red",
        "faecis": "yellow",
        "xylanisolvens": "#a6d854ff",
        "coprophilus": "red"}
    plt.rcParams["font.family"] = "Arial"
    group = 1
    # For each of the clusters on the tree
    for cluster in clusters:
        srs_use = [x for x in cluster if x in df_bacteroides_data.index]
        # Minimum 3 sequences in the cluster to consider
        if len(srs_use) >= 3:
            # Calculate the MSE
            df_error = single_species_model(df_bacteroides_data.loc[srs_use, :], transform="")
            df_error.reset_index(inplace=True)

            plt.rc('font', size=8)
            plt.rc('xtick', labelsize=8)

            if (df_error.ix[0, "Species"] == "xylanisolvens") & (group == 3):
                color_use = "#ffd92fff"
            else:
                color_use = color_dict[df_error.ix[0, "Species"]]

            # Plot a scatter plot for the abundance of the best fitting species against the orphan immunity abundance
            # Figure 2B, D
            f, ax = plt.subplots(figsize=(1.6, 1.6))
            plt.xlim([-9.0, -7.0])
            plt.ylim([-9.0, -7.0])
            ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
            ax.scatter(
                np.log10(df_bacteroides_data.loc[srs_use, "markers_s__Bacteroides_" + DfError.ix[0, "Species"] + ".fasta"] + pseudo),
                np.log10(df_bacteroides_data.loc[srs_use, "OrphanI"] + pseudo), color=color_use, edgecolors="black")
            plt.xlabel("Bacteroides %s" % (DfError.ix[0, "Species"]), style="italic")
            plt.ylabel("Orphan I")
            plt.xticks([-9.0, -8.5, -8.0, -7.5, -7.0])
            plt.savefig("Figure2_Species_%s_Group_%i_%s_Scatter.pdf" % (DfError.ix[0, "Species"], group, gene_oi))
            plt.savefig(
                "Figure2_Species_%s_Group_%i_%s_Scatter_SecondHit.pdf" % (DfError.ix[1, "Species"], group, gene_oi))
            plt.clf()
            f.clf()

            # Plot the MSE of the top 5 species as a barplot
            # Figure 2C, E
            df_error_sub = df_error.ix[:5, :]
            df_error_sub["Species"] = ["B." + x[:2] for x in df_error_sub["Species"]]
            plt.rc('font', size=8)
            plt.rc('xtick', labelsize=6)
            f, ax = plt.subplots(figsize=(0.8, 1.6))
            ind = np.arange(df_error_sub.shape[0])
            ax.bar(ind, df_error_sub["MSE"], edgecolor="black", color="#d3d3d3", linewidth=0.25)
            plt.xticks(ind, df_error_sub["Species"], rotation=90, style="italic")
            plt.title("MSE")
            plt.savefig("Figure2_Species_%s_Group_%i_%s_MSE.pdf" % (DfError.ix[0, "Species"], group, gene_oi))
            plt.clf()
            f.clf()

            plt.rc('font', size=8)
            plt.rc('xtick', labelsize=8)
            group += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assings immunity sequences in metagenomes to their most likely species that encodes it')
    parser.add_argument('--infile_counts_ei', help='A dataframe of counts mapping to GA3 EI genes')
    parser.add_argument('--infile_seqs_ei', help='A fasta file of the sequences represented by infile_counts_ei')
    parser.add_argument('--infile_total_reads', help='A dataframe of the total number of reads in each library')
    parser.add_argument('--infile_bacteroides_counts', help='A dataframe of counts mapping to Bacteroides marker genes')
    parser.add_argument('--infile_bacteroides_seqs', help='A fasta file of the sequences represented by infile_bacteroides_counts')
    parser.add_argument('--indir_marker_genes', help='A directory containing markers_s__Bacteroides_XXX.fasta files')
    parser.add_argument('--indir_pileups', help='A directory containing pileups that you want to build a tree from')
    parser.add_argument('--infile_db', help='A BLAST database of different bacteroides strains genes')
    parser.add_argument('--infile_genome_names', help='A table with the genome names of each B. fragilis strain')
    args = parser.parse_args()

    # Load the data
    df_data = create_bacteroides_abundance_table(args.infile_bacteroides_counts, args.infile_total_reads, args.infile_bacteroides_seqs, args.indir_marker_genes)
    df_data = add_ei_data(df_data, args.infile_ei_counts, args.infile_ei_seqs)

    # Figure 2A
    tree_groups = create_tree(args.indir_pileups, args.infile_seqs_ei, args.infile_db, args.infile_genome_names)
    # Figure 2B-E
    assign_species(df_data, tree_groups)

