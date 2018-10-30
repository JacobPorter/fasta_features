#!/usr/bin/env python
"""
Extract features useful for machine learning from a fasta file.

:Authors:
    Jacob Porter <jsporter@vt.edu>
"""

import argparse
import sys
import datetime
from collections import defaultdict
from itertools import product
from SeqIterator.SeqIterator import SeqReader

# The nucleotide bases.
nucleic_acids_string = "ACGTN"
# The amino acid bases.
amino_acids_string = "ARNDBCEQZGHILKMFPSTWYVXJUO*"


def get_features(file_input, k_length=[3], stride=[1],
                 output=sys.stdout, amino_acids=False,
                 file_type='fasta', freq=False):
    """
    Extract and output the features.

    Parameters
    ----------
    file_input: str
        The location of the input fasta or fastq file.
    k_length: list
        A list of lengths of kmers to sample from.
    stride: list
        The lengths of the stride to take kmers from.  The order
        corresponds to the order in k_length.
    output: writable
        The writable to write the output to.
    amino_acids: bool
        Use this switch to determine if the input is in nucelic acids or
        amino acids.
    file_type: str
        Either 'fasta' or 'fastq'.  This must correspond to the
        file input.

    Returns
    -------
    int
        A count of the number of processed records.

    """
    reader = SeqReader(file_input, file_type=file_type)
    records_processed = 0
    order_dict = {}
    for k in k_length:
        if amino_acids:
            acid_order = product(amino_acids_string, repeat=k)
        else:
            acid_order = product(nucleic_acids_string, repeat=k)
        acid_order = list(acid_order)
        acid_order = ["".join(t1) for t1 in acid_order]
        acid_order.sort()
        order_dict[k] = acid_order
    for record in reader:
        seq = record[1].upper()
        for k, s in zip(k_length, stride):
            kmer_counts, kmer_amount = _count_kmers(seq, k, stride=s)
            for acid in order_dict[k]:
                denominator = kmer_amount
                if not freq:
                    denominator = 1.0
                output.write(str(kmer_counts[acid] / denominator) + "\t")
        output.write("\n")
        records_processed += 1
    return records_processed


def _count_kmers(seq, k, stride=1):
    kmer_counts = defaultdict(int)
    kmer_amount = 0
    for i in range(0, len(seq), stride):
        if len(seq[i:i+k]) != k:
            continue
        kmer_counts[seq[i:i+k]] += 1
        kmer_amount += 1
    return kmer_counts, kmer_amount


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description=__doc__)
    parser.add_argument("file_input", type=str,
                        help=('The input fasta or fastq '
                              'file with the sequences.'))
    parser.add_argument("--k_length", "-k", type=int,
                        action='append',
                        help=('The length of the kmer.'),
                        required=True)
    parser.add_argument("--stride", "-s", type=int,
                        action='append',
                        help=("The strides corresponding to each kmer "
                              "length."),
                        required=True)
    parser.add_argument("--output", "-o", type=str,
                        help=("The file location to write the output to."),
                        default=sys.stdout)
    parser.add_argument("--fastq", '-q', action="store_true",
                        help=("Turn this on if the input file is a "
                              "fastq file."),
                        default=False)
    parser.add_argument("--amino_acids", "-a", action="store_true",
                        help=("Turn this on if amino acids are being used."),
                        default=False)
    parser.add_argument("--freq", "-f", action="store_true",
                        help=("Use a frequency vector instead of a count."),
                        default=False)
    args = parser.parse_args()
    file_type = 'fasta'
    if args.fastq:
        file_type = 'fastq'
    count = get_features(args.file_input, args.k_length, args.stride,
                         args.output, args.amino_acids, file_type, args.freq)
    print("There were {} records processed.".format(count), file=sys.stderr)
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
