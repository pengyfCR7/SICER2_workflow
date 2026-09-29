#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import numpy as np
import scipy
from sicer.lib import GenomeData
from sicer.main import run_SICER, run_SICER_df


# SICER 2.1.0 still uses scipy.array, which recent SciPy versions removed.
if not hasattr(scipy, "array"):
    scipy.array = np.array


def register_genome(path):
    chroms = []
    lengths = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                chrom, length = line.rstrip("\n").split("\t")
                chroms.append(chrom)
                lengths[chrom] = int(length)
    GenomeData.species_chroms["custom"] = chroms
    GenomeData.species_chrom_lengths["custom"] = lengths


def add_common_args(parser):
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output_directory", required=True)
    parser.add_argument("--redundancy_threshold", type=int, required=True)
    parser.add_argument("--window_size", type=int, required=True)
    parser.add_argument("--fragment_size", type=int, required=True)
    parser.add_argument("--effective_genome_fraction", type=float, required=True)
    parser.add_argument("--gap_size", type=int, required=True)
    parser.add_argument("--false_discovery_rate", type=float, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--e_value", type=int, default=1000)
    parser.add_argument("--significant_reads", action="store_true")
    parser.add_argument("--verbose", action="store_true")


def parse_args():
    parser = argparse.ArgumentParser(description="Run SICER2 with custom chromosome sizes")
    commands = parser.add_subparsers(dest="mode", required=True)

    call = commands.add_parser("call")
    add_common_args(call)
    call.add_argument("--treatment_file", required=True)
    call.add_argument("--control_file")

    diff = commands.add_parser("diff")
    add_common_args(diff)
    diff.add_argument("--treatment_file", nargs=2, required=True)
    diff.add_argument("--control_file", nargs=2)
    diff.add_argument("--false_discovery_rate_df", type=float, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    start_dir = Path.cwd()
    args.output_directory = str(Path(args.output_directory).resolve())
    Path(args.output_directory).mkdir(parents=True, exist_ok=True)

    if args.mode == "call":
        args.treatment_file = str(Path(args.treatment_file).resolve())
        args.control_file = str(Path(args.control_file).resolve()) if args.control_file else None
    else:
        args.treatment_file = [str(Path(path).resolve()) for path in args.treatment_file]
        args.control_file = (
            [str(Path(path).resolve()) for path in args.control_file]
            if args.control_file else None
        )

    register_genome(args.chrom_sizes)
    args.species = "custom"
    args.input_type = "SE"
    args.subcommand = "SICER"
    args.df = args.mode == "diff"

    try:
        if args.mode == "call":
            run_SICER.main(args)
        else:
            run_SICER_df.main(args)
    finally:
        # SICER changes into a temporary directory and does not restore cwd.
        os.chdir(start_dir)


if __name__ == "__main__":
    main()
