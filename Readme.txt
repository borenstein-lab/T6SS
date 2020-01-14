
Source code for the computational analyses in Ross BD, Verster AJ et al. 2019.

Most relevant data is included in this repository, excluding the file:

GCA_AllStrains_mRNA_Curated_Alt_StopFix.fasta

which was too large to host in a standard GitHub repository. To obtain this file, please contact elbo@tauex.tau.ac.il

Lib/Figure1.py

Plots the relationship between GA3 gene abundance and the Bacteroides fragilis
abundance data

  --infile_counts_ei INFILE_COUNTS_EI
                        A dataframe of counts mapping to GA3 EI genes
  --infile_seqs_ei INFILE_SEQS_EI
                        A fasta file of the sequences represented by
                        infile_counts_ei
  --infile_total_reads INFILE_TOTAL_READS
                        A dataframe of the total number of reads in each
                        library
  --infile_bacteroides_counts INFILE_BACTEROIDES_COUNTS
                        A dataframe of counts mapping to Bacteroides marker
                        genes
  --infile_bacteroides_seqs INFILE_BACTEROIDES_SEQS
                        A fasta file of the sequences represented by
                        infile_bacteroides_counts
  --indir_marker_genes INDIR_MARKER_GENES
                        A directory containing
                        markers_s__Bacteroides_XXX.fasta files

Lib/Figure2.py

Assings immunity sequences in metagenomes to their most likely species that
encodes it

  --infile_counts_ei INFILE_COUNTS_EI
                        A dataframe of counts mapping to GA3 EI genes
  --infile_seqs_ei INFILE_SEQS_EI
                        A fasta file of the sequences represented by
                        infile_counts_ei
  --infile_total_reads INFILE_TOTAL_READS
                        A dataframe of the total number of reads in each
                        library
  --infile_bacteroides_counts INFILE_BACTEROIDES_COUNTS
                        A dataframe of counts mapping to Bacteroides marker
                        genes
  --infile_bacteroides_seqs INFILE_BACTEROIDES_SEQS
                        A fasta file of the sequences represented by
                        infile_bacteroides_counts
  --indir_marker_genes INDIR_MARKER_GENES
                        A directory containing
                        markers_s__Bacteroides_XXX.fasta files
  --indir_pileups INDIR_PILEUPS
                        A directory containing pileups that you want to build
                        a tree from
  --infile_db INFILE_DB
                        A BLAST database of different bacteroides strains
                        genes
  --infile_genome_names INFILE_GENOME_NAMES
                        A table with the genome names of each B. fragilis
                        strain


