#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import gzip
from collections import defaultdict


def main():

    parser = argparse.ArgumentParser(
        description="strips polyA followed by rest of trailing sequence from end of reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--fq_file", required=True, type=str, help="fastq file")

    parser.add_argument(
        "--min_polyA_len",
        default=10,
        help="min length of polyA stretch required for pruning",
    )

    parser.add_argument(
        "--terminal_frac_seq_search",
        type=float,
        default="0.5",
        help="search only in this tail fraction of the read sequence for polyA"
    )

    
    parser.add_argument(
        "--print_trimmed_only",
        action='store_true',
        default=False,
        help="only print those fastq records that got trimmed."
    )
    
    
    args = parser.parse_args()
    
    fq_file = args.fq_file
    min_polyA_len = args.min_polyA_len
    print_trimmed_only_flag = args.print_trimmed_only    
    
    
    fq_iterator = fastq_iterator(fq_file)
    terminal_frac = args.terminal_frac_seq_search

    record_counter = 0
    trim_counter = 0
    
    for read_tuple in fq_iterator:
        readname, readseq, L3, quals = read_tuple

        record_counter += 1
        
        readseq = readseq.upper()
        
        readseq_len = len(readseq)

        A_string = "A" * min_polyA_len

        read_search_start = readseq_len - int(terminal_frac * readseq_len)

        # search polyA
        polyA_pos = readseq.find(A_string, read_search_start)
        found_polyA = False
        if polyA_pos >= read_search_start:
            found_polyA = True
            readseq = readseq[0:polyA_pos]
            quals = quals[0:polyA_pos]
            trim_len = readseq_len - polyA_pos
            polyA_pos += 1 # for 1-based coordinates
            L3 += f"polyA:{polyA_pos};trimLen:{trim_len}"
            trim_counter += 1

        if  print_trimmed_only_flag and not found_polyA:
            continue

        # print record
        print(
            "\n".join([readname, readseq, L3, quals]))


    pct_reads_trimmed = trim_counter / record_counter * 100
    print(f"# {fq_file} : trimmed {trim_counter} / {record_counter} = {pct_reads_trimmed:.1f}%", file=sys.stderr)
    
    sys.exit(0)



def fastq_iterator(fastq_filename):

    if re.search(".gz$", fastq_filename):
        fh = gzip.open(fastq_filename, "rt", encoding="utf-8")
    else:
        fh = open(fastq_filename, "rt", encoding="utf-8")

    have_records = True
    while have_records:
        try:
            readname = next(fh).rstrip()
            readseq = next(fh).rstrip()
            L3 = next(fh).rstrip()
            quals = next(fh).rstrip()

            yield (readname, readseq, L3, quals)

        except StopIteration:
            have_records = False

    return




if __name__ == "__main__":
    main()
