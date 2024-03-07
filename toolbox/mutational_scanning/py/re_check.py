# == Native Modules
import os.path
# == Installed Modules
import pandas as pd


def reverser_complement(input_sequence):
    pairs = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    reverse_complement_seq = ""

    for index in range(len(input_sequence) - 1, 0, -1):
        base = input_sequence[index]
        complement = pairs[base]
        reverse_complement_seq += complement

    return reverse_complement_seq


def avoidREs_nodelchar(filename: str, seqrecord_list: list):
    """
    Adapted from Brittney Thornton's original code (03/2024)
    == Restriction Enzyme check ==
    This version of avoidREs was built to work with the assign_variant_name fxn.
    It takes the filename of the .csv file (filename -- string) that stores all
    the REs and the recognition sites to avoid, as well as an array of arrays,
    where the first column in each array is filename, and the second column is
    the sequence in AGCT.
    The output is two arrays or arrays: one containing passed
    filenames and seqs, another containing failed filenames and seqs.
    Takes a filename of a CSV with a restriction enzyme list and an input list of
    strings to sort through. If any RE recognition sites are found in each string,
    the string gets added to faillist output. If no RE recognition sites are found
    in the string, it gets added to the passlist output.
    :param filename:
    :param seqrecord_list:
    :return:
    """
    passlist = []
    faillist = []
    refaillist = []
    determinator = []
    # Check paths
    if not os.path.exists(filename):
        raise "Check your restriction enzyme path as it couldn't be found"

    # Import Restriction Enzyme dataframe
    re_df = pd.read_csv(filename, names=["RE", "Motif"])
    res = list(re_df.Motif) + [reverser_complement(x) for x in list(re_df.Motif)]
    # Import Fasta sequences
    for record in seqrecord_list:
        nt_sequence = record.seq
        for re in res:
            if re in nt_sequence:
                determinator.append(True)
                refaillist.append(re)
                break
            if re not in nt_sequence:
                determinator.append(False)
        if any(determinator):
            faillist.append(record)
        else:
            passlist.append(record)
        determinator = []

    return passlist, faillist, refaillist
