#!/usr/bin/env python3
"""
Take two input csv files, one for BAMs and one for FASTQs, with paths to the
sequencing data to process.  This script will set up symlinks for the project to
those files.

The BAM listing must have two columns: sample and bam.
There should be one record per sample.
The FASTQ listing must have seven columns: sample, platform, sequencing_center,
date, library, read_one_fastq, and read_two_fastq.
The platform is the sequencing platform (typically Illumina),
the sequencing_center is where the samples were sequenced, e.g. RCIGM,
the date should be the approximate date of the sequencing, e.g. 20211007,
and the library is an identifier for the library that was used to prepare the
sample for sequencing.  The library can simply be something like "lib1" for all
samples if this is not useful information to retain.

There should be one record per read group, so each sample will have one or more
records.

Note that the first line for each must have a header with these field names and
that one of the input files can have no records (if, say, your sample set only
has FASTQs as input, you will have no BAM records but must still pass in a BAM
listing with no records).

The script will check for various errors, e.g. missing files, samples being
listed in both files, invalid sample names, etc.
"""
import argparse
import sys
import os
import re
from collections import defaultdict, namedtuple
from glob import glob
from yaml import safe_load

BAM_LISTING_FIELDS = ["sample", "bam"]
FASTQ_LISTING_FIELDS = [
    "sample",
    "platform",
    "sequencing_center",
    "date",
    "library",
    "read_one_fastq",
    "read_two_fastq",
]


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def setup_inputs(bam_listing_fh, fastq_listing_fh):
    try:
        bam_listing_fields = next(bam_listing_fh).strip().split(",")
    except StopIteration:
        raise ValueError(
            f"The BAM listing file ({bam_listing_fh.name}) "
            "must have at least a header."
        )
    if bam_listing_fields != BAM_LISTING_FIELDS:
        raise ValueError(
            f"The BAM listing file ({bam_listing_fh.name}) "
            "does not have the proper header "
            f"({','.join(BAM_LISTING_FIELDS)}."
        )
    try:
        fastq_listing_fields = next(fastq_listing_fh).strip().split(",")
    except StopIteration:
        raise ValueError(
            f"The FASTQ listing file ({fastq_listing_fh.name}) "
            "must have at least a header."
        )
    if fastq_listing_fields != FASTQ_LISTING_FIELDS:
        raise ValueError(
            f"The FASTQ listing file ({fastq_listing_fh.name}) "
            "does not have the proper header:\n"
            "{}".format(",".join(FASTQ_LISTING_FIELDS))
        )
    bams_in_input_directory = bool(glob("../input/bams/*.bam"))
    fastqs_in_input_directory = bool(glob("../input/fastqs/*.fastq.gz"))
    fastq_rgs_in_input_directory = bool(glob("../input/fastqs/*.txt"))
    if (
        bams_in_input_directory
        or fastqs_in_input_directory
        or fastq_rgs_in_input_directory
    ):
        raise ValueError(
            "The input directories (../input/bams and "
            "../input/fastqs) should not have data already."
        )
    config_path = "../config/config.yaml"
    if os.path.isfile(config_path):
        with open(config_path) as config_fh:
            config = safe_load(config_fh)
        sample_regex = re.compile(config["sample_regex"])
    else:
        raise ValueError("Couldn't find project config file ({config_path}).")

    bam_samples, fastq_samples = set(), set()
    bam_mapping = {}
    fastq_mapping = defaultdict(list)
    bam_samples_seen_multiple_times = set()
    missing_input_files = []
    invalid_sample_names = set()
    for line in bam_listing_fh:
        sample, bam = line.strip().split(",")
        if not sample_regex.match(sample):
            invalid_sample_names.add(sample)
        if sample in bam_samples:
            bam_samples_seen_multiple_times.add(sample)
        else:
            bam_samples.add(sample)
            bam_mapping[sample] = bam
            if not os.path.isfile(bam):
                missing_input_files.append(bam)
    FASTQRecord = namedtuple(
        "FASTQRecord",
        (
            "sample",
            "platform",
            "sequencing_center",
            "date",
            "library",
            "read_one",
            "read_two",
        ),
    )
    for line in fastq_listing_fh:
        fastq_record = FASTQRecord(*line.strip().split(","))
        if not sample_regex.match(fastq_record.sample):
            invalid_sample_names.add(fastq_record.sample)
        fastq_samples.add(fastq_record.sample)
        fastq_mapping[fastq_record.sample].append(fastq_record)
        if not os.path.isfile(fastq_record.read_one):
            missing_input_files.append(fastq_record.read_one)
        if not os.path.isfile(fastq_record.read_two):
            missing_input_files.append(fastq_record.read_two)
    all_samples = bam_samples | fastq_samples
    if not all_samples:
        raise ValueError("Both the BAM and FASTQ listings have no records?")
    overlap = bam_samples & fastq_samples
    if invalid_sample_names:
        raise ValueError(
            "The following samples have invalid names "
            "(compared against the regular expression "
            f"{config['sample_regex']}):\n"
            "{{}}".format("\n".join(sorted(invalid_sample_names)))
        )
    if overlap:
        raise ValueError(
            "The following samples are in the BAM and FASTQ "
            f"listings:\n"
            "{}".format("\n".join(sorted(overlap)))
        )
    if bam_samples_seen_multiple_times:
        raise ValueError(
            "The following samples are in the BAM listing "
            "multiple times:\n"
            "{}".format("\n".join(sorted(bam_samples_seen_multiple_times)))
        )
    if missing_input_files:
        raise ValueError(
            "The following input files are missing:\n"
            "{}".format("\n".join(missing_input_files))
        )
    print(
        f"Found {len(all_samples)} samples and no errors detected.  Creating symlinks now."
    )
    symlinks_created = []
    rg_lines_created = []
    try:
        for sample, bam in bam_mapping.items():
            link = f"../input/bams/{sample}.bam"
            os.symlink(bam, link)
            symlinks_created.append(link)
        for sample_fastq_records in fastq_mapping.values():
            for x, fastq_record in enumerate(sample_fastq_records):
                link = f"../input/fastqs/{fastq_record.sample}_{x}.R1.fastq.gz"
                os.symlink(fastq_record.read_one, link)
                symlinks_created.append(link)
                link = f"../input/fastqs/{fastq_record.sample}_{x}.R2.fastq.gz"
                os.symlink(fastq_record.read_two, link)
                symlinks_created.append(link)
                rg_line = f"../input/fastqs/{fastq_record.sample}_{x}.txt"
                with open(rg_line, "w") as rg_fh:
                    rg_fh.write(
                        f"@RG\tID:{x}\tLB:{fastq_record.library}\tSM:"
                        f"{fastq_record.sample}\tPL:"
                        f"{fastq_record.platform}\tCN:"
                        f"{fastq_record.sequencing_center}\tDT:"
                        f"{fastq_record.date}\n"
                    )
                rg_lines_created.append(rg_line)
    except:
        # remove any already created symlinks if there is any error
        for link in symlinks_created:
            os.remove(link)
        for rg_line in rg_lines_created:
            os.remove(rg_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter
    )
    parser.add_argument(
        "-b",
        "--bam-listing",
        default="bam_listing.csv",
        type=argparse.FileType("r"),
        help="the BAM input listing",
    )
    parser.add_argument(
        "-f",
        "--fastq-listing",
        default="fastq_listing.csv",
        type=argparse.FileType("r"),
        help="the FASTQ input listing",
    )
    args = parser.parse_args()
    setup_inputs(args.bam_listing, args.fastq_listing)
