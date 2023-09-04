import pybedtools
import SeqIO from Bio
import numpy as np
import pickle

def get_fasta_from_bed(bed_file, fasta_file):
    """
    Extract sequences from fasta file using bed file
    """
    bed = pybedtools.BedTool(bed_file)
    fasta = pybedtools.BedTool(fasta_file)
    fasta_bed = fasta.sequence(fi=fasta_file, bed=bed_file)
    return fasta_bed

def truncate_fasta(fasta):
    """
    Truncate FASTA sequences to middle 200 bp regardless of length
    """
    fasta_trunc = []
    for record in SeqIO.parse(fasta, "fasta"):
        seq_len = len(record.seq)
        seq_start = int((seq_len - 200) / 2)
        seq_end = seq_start + 200
        record.seq = record.seq[seq_start:seq_end]
        fasta_trunc.append(record)
    return fasta_trunc

def preprocess_motif(motifs, bed_files, neg_file):
    """
    Preprocess motif and negative control sequences
    """
    ## Import BED File with Negative Control Locations
    neg_file = "neg.bed"
    ## Import FASTA File for D. Melanogaster
    fasta_file = "dmel-all-chromosome-r6.17.fasta"

    motif_dict = {}

    for motif, bed_file in zip(motifs, bed_files):
        ## Get FASTA Sequences for Motif Locations
        motif_fasta = get_fasta_from_bed(bed_file, fasta_file)

        ## Get FASTA Sequences for Negative Control Locations
        neg_fasta = get_fasta_from_bed(neg_file, fasta_file)

        ## Truncate FASTA Sequences to middle 200 bp Regardless of Length
        motif_fasta_trunc = truncate_fasta(motif_fasta)

        ## Concatenate Both FASTA Sequences
        fasta = motif_fasta_trunc + neg_fasta

        ## Convert FASTA Sequences to List of Strings and Convert to NP Array
        fasta_list = []
        for record in fasta:
            fasta_list.append(str(record.seq))
        fasta_np = np.array(fasta_list)

        ## Create Dictionary with Motif as Key and FASTA Sequences as Values
        motif_dict[motif] = fasta_np

    # Dump motif_dict to pickle file
    with open("motif_dict.pkl", "wb") as f:
        pickle.dump(motif_dict, f)