#!/usr/bin/env python

import argparse
import logging

import edlib
import pysam
import sys
from tqdm import tqdm

import os


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def str2bool(v):
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )
    parser.add_argument(
        "--max-error",
        dest="MAX_ERROR",
        type=int,
        default=2,
        help="Max edit distance for UMI",
    )
    parser.add_argument(
        "--adapter-length",
        dest="ADAPTER_LENGTH",
        type=int,
        default=250,
        help="Length of adapter",
    )
    parser.add_argument(
        "-t", "--theads", dest="THREADS", type=int, default=1, help="Number of threads."
    )
    parser.add_argument(
        "--tsv", dest="TSV", type=str, required=False, help="TSV output file"
    )
    parser.add_argument(
        "-o", "--output", dest="OUT", type=str, required=False, help="FASTA output file"
    )
    parser.add_argument(
        "-ifwd", dest="INPUT_FAFWD", type=str, default="/dev/stdin", help="Fasta input with reads mapped on + strand LTR"
    )
    parser.add_argument(
        "-irev", dest="INPUT_FAREV", type=str, default="/dev/stdin", help="Fasta input with reads mapped on - strand LTR"
    )

    args = parser.parse_args(argv)

    return args


def align(query, pattern_info, max_ed, normalise=False):
    pattern, wildcard, equalities, forward = pattern_info

    result = edlib.align(
        pattern,
        query,
        task="path",
        mode="HW",
        k=max_ed,
        additionalEqualities=equalities,
    )
    if result["editDistance"] == -1:
        return None, None

    ed = result["editDistance"]
    if not normalise:
        locs = result["locations"][0]
        umi = query[locs[0]:locs[1]+1]
        return ed, umi

    # Extract and normalise UMI
    umi = ""
    align = edlib.getNiceAlignment(result, pattern, query)
    for q, t in zip(align["query_aligned"], align["target_aligned"]):
        if q not in wildcard:
            continue
        if t == "-":
            umi += "N"
        else:
            umi += t

    if len(umi) != 16:
        raise RuntimeError("UMI length incorrect: {}".format(umi))

    return ed, umi


def count_reads_fastx(fasta_filename):
    n_read = 0
    logging.info("Counting reads in {}".format(fasta_filename))
    with pysam.FastxFile(fasta_filename) as fh:
        for entry in fh:
            n_read += 1
    return n_read


def extract_adapters(entry, max_adapter_length):
    read_5p_seq = None
    read_3p_seq = None
    if len(entry.sequence) > max_adapter_length:
        read_5p_seq = entry.sequence[:max_adapter_length]
        read_3p_seq = entry.sequence[-max_adapter_length:]
    return read_5p_seq, read_3p_seq


def print_seq(
    entry,
    strand,
    result_5p_fwd_dist,
    result_3p_rev_dist,
    result_5p_fwd_seq,
    result_3p_rev_seq,
    out
):
    if (result_5p_fwd_seq == None) and (result_3p_rev_seq != None):
        seq = result_3p_rev_seq
        print(
            ">{};strand={};umi_fwd={};umi_rev={};umi_fwd_seq={};umi_rev_seq={};seq={}".format(
                entry.name,
                strand,
                None,
                result_3p_rev_dist,
                None,
                result_3p_rev_seq,
                entry.sequence,
            ),
            file=out,
        )
        print(seq, file=out)

    elif (result_3p_rev_seq == None) and (result_5p_fwd_seq != None):
        seq = result_5p_fwd_seq
        print(
            ">{};strand={};umi_fwd={};umi_rev={};umi_fwd_seq={};umi_rev_seq={};seq={}".format(
                entry.name,
                strand,
                result_5p_fwd_dist,
                None,
                result_5p_fwd_seq,
                None,
                entry.sequence,
            ),
            file=out,
        )
        print(seq, file=out)


def extract_umis(
    input_file_fwd,
    input_file_rev,
    max_adapter_length,
    max_pattern_dist,
    pattern_fwd,
    pattern_rev,
    output_file,
    out,
    tsv,
    mode
):

    n_umi_fwd = 0
    n_umi_rev = 0
    n_umi = 0
    strand_stats = {"+": 0, "-": 0}
    n_read_fwd = count_reads_fastx(input_file_fwd)
    n_read_rev = count_reads_fastx(input_file_rev)
    n_read =  n_read_fwd + n_read_rev
    with tqdm(total=n_read) as pbar:
            with pysam.FastxFile(input_file_fwd) as fh:
                
                for entry in fh:
                    pbar.update(1)

                    strand = "+"
                    strand_stats[strand] += 1

                    read_5p_seq,read_3p_seq = extract_adapters(entry, max_adapter_length)

                    if mode == "LTR3" and read_3p_seq:
                        result_3p_fwd_dist, result_3p_fwd_seq = align(read_3p_seq, pattern_rev, max_pattern_dist)
                        if result_3p_fwd_seq:
                            n_umi += 1
                            n_umi_fwd += 1
                            print_seq(entry, strand, None, result_3p_fwd_dist, None, result_3p_fwd_seq, out)
                    elif mode == "LTR5" and read_5p_seq:
                        result_5p_fwd_dist, result_5p_fwd_seq = align(read_5p_seq, pattern_fwd, max_pattern_dist)
                        if result_5p_fwd_seq:
                            n_umi += 1
                            n_umi_fwd += 1
                            print_seq(entry, strand, result_5p_fwd_dist, None, result_5p_fwd_seq, None, out)


            with pysam.FastxFile(input_file_rev) as fh:
                
                for entry in fh:
                    pbar.update(1)

                    strand = "-"
                    strand_stats[strand] += 1

                    read_5p_seq, read_3p_seq = extract_adapters(entry, max_adapter_length)
                    
                    if mode == "LTR3" and read_3p_seq:
                        result_3p_rev_dist, result_3p_rev_seq = align(read_3p_seq, pattern_rev, max_pattern_dist)
                        if result_3p_rev_seq:
                            n_umi += 1
                            n_umi_rev += 1
                            print_seq(entry, strand, None, result_3p_rev_dist, None, result_3p_rev_seq, out)
                    elif mode == "LTR5" and read_5p_seq:
                        result_5p_rev_dist, result_5p_rev_seq = align(read_5p_seq, pattern_fwd, max_pattern_dist)
                        if result_5p_rev_seq:
                            n_umi += 1
                            n_umi_rev += 1
                            print_seq(entry, strand, result_5p_rev_dist, None, result_5p_rev_seq, None, out)

    fwd_rev_ratio = -1
    perc = -1.0
    fwd_rev_ratio = strand_stats["+"] / strand_stats["-"]
    logging.info(
        "Found {} fwd and {} rev reads (ratio: {})".format(
            strand_stats["+"], strand_stats["-"], fwd_rev_ratio
        )
    )
    if n_read:
        perc = 100.0 * (n_umi / n_read)
        perc_fwd = 100.0 * (n_umi_fwd / n_read_fwd)
        perc_rev = 100.0 * (n_umi_rev / n_read_rev)
        logging.info(
            "{}% of reads contained an UMI with max {} mismatches and {} ({}%) of fwd reads and {} ({}%) rev reads.".format(
                perc, max_pattern_dist, n_umi_fwd, perc_fwd, n_umi_rev, perc_rev
            )
        )
    if tsv:
        print(
            output_file,
            max_pattern_dist,
            strand_stats["+"],
            strand_stats["-"],
            fwd_rev_ratio,
            perc,
            file=tsv,
        )

    out.close()


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    max_adapter_length = args.ADAPTER_LENGTH
    max_pattern_dist = args.MAX_ERROR
    output_file = args.OUT
    tsv_file = args.TSV
    input_file_fwd = args.INPUT_FAFWD
    input_file_rev = args.INPUT_FAREV

    if "LTR3" in os.path.basename(input_file_fwd) and "LTR3" in os.path.basename(input_file_rev):
        mode = "LTR3"
    elif "LTR5" in os.path.basename(input_file_fwd) and "LTR5" in os.path.basename(input_file_rev):
        mode = "LTR5"
    else:
        raise ValueError("Error: Both input fasta files must contain either 'LTR5' or 'LTR3' in their names.")


    pattern_fwd = (
        "TTTVVVVTTVVVVTTVVVVTTVVVVTTT",
        "V",
        [("V", "A"), ("V", "G"), ("V", "C")],
        True,
    )
    pattern_rev = (
        "AAABBBBAABBBBAABBBBAABBBBAAA",
        "B",
        [("B", "T"), ("B", "G"), ("B", "C")],
        False,
    )

    tsv = None
    if tsv_file:
        tsv = open(tsv_file, "w")

    with open(output_file, "w") as out:
        extract_umis(
            input_file_fwd,
            input_file_rev,
            max_adapter_length,
            max_pattern_dist,
            pattern_fwd,
            pattern_rev,
            output_file,
            out,
            tsv,
            mode
        )
    if tsv_file:
        tsv.close()

if __name__ == "__main__":
    main()
