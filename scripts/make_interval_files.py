#!/usr/bin/env python3
"""
Take a dict file for a reference genome and create a series of interval files
containing at most the specified number of bases per file.
Only takes primary chromosomes, mitochondria, and alt chromosome contigs
(not decoy).
"""
import argparse
import sys
import re
from collections import defaultdict

DEFAULT_DICT_FILE = "/projects/ps-gleesonlab3/resources/hg38/Homo_sapiens_assembly38.dict"
ORDER = [f"chr{chromosome}" for chromosome in list(range(1,23)) + ["X", "Y", "M"]]

def make_interval_files(dict_fh, nbases):
    alt_contigs_by_chromosome = defaultdict(dict)
    chromosome_lengths = dict()
    for line in dict_fh:
        if line.startswith("@SQ"):
            try:
                d = dict([entry.split(":", 1) for entry in line.strip().split("\t")[1:]])
            except:
                print(line)
                raise
            if d["SN"] in ORDER:
                chromosome_lengths[d["SN"]] = int(d["LN"])
            elif d["SN"].endswith("_alt"):
                alt_contigs_by_chromosome[d["SN"].split("_")[0]][d["SN"]] = int(d["LN"])

    current_length = 0
    current_index = 0
    sequences = ORDER[:]
    for chromosome in ORDER:
        for key in sorted(alt_contigs_by_chromosome[chromosome].keys()):
            sequences.append(key)
            chromosome_lengths[key] = alt_contigs_by_chromosome[chromosome][key]
    order = iter(sequences)
    current_chromosome = next(order)
    current_position = 0
    length_remaining_current_intervalfile = nbases
    fh = open(f"intervalfile_{current_index}.list", "w")
    while True:
        remaining_in_chromosome = chromosome_lengths[current_chromosome] - current_position
        current_position += 1
        if remaining_in_chromosome >= length_remaining_current_intervalfile:
            ending_position = current_position + length_remaining_current_intervalfile - 1
            fh.write(f"{current_chromosome}:{current_position}-{ending_position}\n")
            fh.close()
            if remaining_in_chromosome == length_remaining_current_intervalfile:
                try:
                    current_chromosome = next(order)
                    current_position = 0
                    remaining_in_chromosome = chromosome_lengths[current_chromosome]
                except StopIteration:
                    done = True
                    break
            else:
                current_position = ending_position
            current_index += 1
            fh = open(f"intervalfile_{current_index}.list", "w")
            length_remaining_current_intervalfile = nbases
        else:
            ending_position = current_position + remaining_in_chromosome - 1
            fh.write(f"{current_chromosome}:{current_position}-{ending_position}\n")
            length_remaining_current_intervalfile -= remaining_in_chromosome
            try:
                current_chromosome = next(order)
                current_position = 0
                remaining_in_chromosome = chromosome_lengths[current_chromosome]
            except StopIteration:
                break
    if not fh.closed:
        fh.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--dict", default=DEFAULT_DICT_FILE,
                        type=argparse.FileType("r"),
                        help="the .dict file to read")
    parser.add_argument("--megabases", default=30, type=int,
                        help="the number of megabases to include in each interval")
    args = parser.parse_args()
    make_interval_files(args.dict, args.megabases * 10**6)
