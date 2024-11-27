#!/usr/bin/env python3

import sys, os, re
import pysam
import argparse
import logging
import csv


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="add sqanti read classifications into the bam file under the XS tag",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam",
        type=str,
        required=True,
        help="input bam file",
    )

    parser.add_argument(
        "--classification_tsv",
        type=str,
        required=True,
        help="input sqanti3 classification file",
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output bam file containing XS tags",
    )

    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    classification_tsv_filename = args.classification_tsv

    ## get sqanti3 classifications for read names

    read_classifications = dict()
    logger.info(
        "-parsing sqanti3 classifications from {}".format(classification_tsv_filename)
    )
    with open(classification_tsv_filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            read_name = row["isoform"]
            category = row["structural_category"]
            read_classifications[read_name] = category

    # read/write bam, add categories

    bamreader = pysam.AlignmentFile(input_bam_filename, "rb")

    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    logger.info("-reading bam: {}".format(input_bam_filename))
    logger.info("-writing bam: {}".format(output_bam_filename))

    for read in bamreader.fetch():
        read_name = read.query_name
        category = "None"
        if read_name in read_classifications:
            category = read_classifications[read_name]

        read.set_tag("XS", category, "Z")
        bamwriter.write(read)

    bamreader.close()
    bamwriter.close()

    logger.info("-done. See output bam: {}".format(output_bam_filename))

    sys.exit(0)


if __name__ == "__main__":
    main()
