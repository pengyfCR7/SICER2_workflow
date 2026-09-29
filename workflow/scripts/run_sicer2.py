#!/usr/bin/env python3
import argparse
import math
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import scipy
from sicer.lib import GenomeData
from sicer.main import run_SICER, run_SICER_df


# SICER 2.1.0 still calls scipy.array, removed from recent SciPy releases.
# Keep the compatibility fix local to this process; do not modify /opt packages.
if not hasattr(scipy, "array"):
    scipy.array = np.array


def load_genome(path, name="workflow_custom"):
    chroms = []
    lengths = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Invalid chrom sizes line {line_number}: {line.rstrip()}")
            chrom, length_text = fields
            length = int(length_text)
            if length <= 0 or chrom in lengths:
                raise ValueError(f"Invalid or duplicate chromosome at line {line_number}")
            chroms.append(chrom)
            lengths[chrom] = length
    if not chroms:
        raise ValueError("No chromosomes found")
    GenomeData.species_chroms[name] = chroms
    GenomeData.species_chrom_lengths[name] = lengths
    return name


def common_namespace(args, species, is_diff):
    return argparse.Namespace(
        treatment_file=args.treatment_file if is_diff else str(Path(args.treatment_file).resolve()),
        control_file=args.control_file,
        input_type="SE",
        species=species,
        redundancy_threshold=args.redundancy_threshold,
        window_size=args.window_size,
        fragment_size=args.fragment_size,
        effective_genome_fraction=args.effective_genome_fraction,
        false_discovery_rate=args.false_discovery_rate,
        false_discovery_rate_df=getattr(args, "false_discovery_rate_df", None),
        output_directory=str(Path(args.output_directory).resolve()),
        gap_size=args.gap_size,
        e_value=args.e_value,
        cpu=args.cpu,
        significant_reads=args.significant_reads,
        verbose=args.verbose,
        subcommand="SICER",
        df=is_diff,
    )


def find_one(directory, pattern):
    matches = list(Path(directory).glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one file matching {pattern!r}, found {matches}")
    return matches[0]


def write_call_outputs(args):
    raw = Path(args.output_directory)
    if args.control_file:
        source_islands = find_one(
            raw,
            f"*-W{args.window_size}-G{args.gap_size}-FDR{args.false_discovery_rate}-island.bed",
        )
        source_summary = find_one(raw, f"*-W{args.window_size}-G{args.gap_size}-islands-summary")
        shutil.copyfile(source_islands, args.stable_islands)
        header = "chrom\tstart\tend\ttreatment_reads\tcontrol_reads\tp_value\tfold_enrichment\tfdr\n"
        with open(args.stable_summary, "w", encoding="utf-8") as output:
            output.write(header)
            with source_summary.open(encoding="utf-8") as source:
                shutil.copyfileobj(source, output)
    else:
        source = find_one(raw, f"*-W{args.window_size}-G{args.gap_size}.scoreisland")
        shutil.copyfile(source, args.stable_islands)
        with open(args.stable_summary, "w", encoding="utf-8") as output:
            output.write("chrom\tstart\tend\tisland_score\n")
            with source.open(encoding="utf-8") as handle:
                shutil.copyfileobj(handle, output)


def read_diff_summary(path):
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        rows = []
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, values)))
    return header, rows


def write_table(path, header, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in header) + "\n")


def write_diff_outputs(args):
    raw = Path(args.output_directory)
    source = find_one(raw, f"*-and-*-W{args.window_size}-G{args.gap_size}-summary")
    original_header, rows = read_diff_summary(source)
    ratio_column = "Fc_A_vs_B"
    fdr_ab_column = "FDR_A_vs_B"
    fdr_ba_column = "FDR_B_vs_A"
    required = {ratio_column, fdr_ab_column, fdr_ba_column}
    if not required.issubset(original_header):
        raise RuntimeError(f"Differential summary is missing columns: {sorted(required - set(original_header))}")

    header = original_header + [
        "fold_change_test_vs_reference",
        "log2_fold_change_test_vs_reference",
    ]
    for row in rows:
        ratio = float(row[ratio_column])
        row["fold_change_test_vs_reference"] = f"{ratio:.12g}"
        row["log2_fold_change_test_vs_reference"] = f"{math.log2(ratio):.12g}"

    increased_fdr = [row for row in rows if float(row[fdr_ab_column]) <= args.false_discovery_rate_df]
    decreased_fdr = [row for row in rows if float(row[fdr_ba_column]) <= args.false_discovery_rate_df]
    increased_filtered = [
        row for row in increased_fdr
        if float(row["fold_change_test_vs_reference"]) >= args.min_fold_change
    ]
    decreased_filtered = [
        row for row in decreased_fdr
        if float(row["fold_change_test_vs_reference"]) <= 1.0 / args.min_fold_change
    ]
    write_table(args.all_islands, header, rows)
    write_table(args.increased_fdr, header, increased_fdr)
    write_table(args.decreased_fdr, header, decreased_fdr)
    write_table(args.increased_filtered, header, increased_filtered)
    write_table(args.decreased_filtered, header, decreased_filtered)


def add_common(parser):
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output-directory", required=True)
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
    parser = argparse.ArgumentParser(description="Run SICER2 with a custom genome")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    call = subparsers.add_parser("call")
    add_common(call)
    call.add_argument("--treatment_file", required=True)
    call.add_argument("--control_file")
    call.add_argument("--stable-islands", required=True)
    call.add_argument("--stable-summary", required=True)

    diff = subparsers.add_parser("diff")
    add_common(diff)
    diff.add_argument("--treatment_file", nargs=2, required=True)
    diff.add_argument("--control_file", nargs=2)
    diff.add_argument("--false_discovery_rate_df", type=float, required=True)
    diff.add_argument("--min-fold-change", type=float, required=True)
    diff.add_argument("--all-islands", required=True)
    diff.add_argument("--increased-fdr", required=True)
    diff.add_argument("--decreased-fdr", required=True)
    diff.add_argument("--increased-filtered", required=True)
    diff.add_argument("--decreased-filtered", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    start_directory = Path.cwd()
    path_outputs = (
        "stable_islands", "stable_summary", "all_islands", "increased_fdr",
        "decreased_fdr", "increased_filtered", "decreased_filtered",
    )
    for attr in path_outputs:
        value = getattr(args, attr, None)
        if value:
            setattr(args, attr, str(Path(value).resolve()))
    args.output_directory = str(Path(args.output_directory).resolve())
    Path(args.output_directory).mkdir(parents=True, exist_ok=True)
    for attr in path_outputs:
        value = getattr(args, attr, None)
        if value:
            Path(value).parent.mkdir(parents=True, exist_ok=True)
    species = load_genome(args.chrom_sizes)

    if args.mode == "call":
        if args.control_file:
            args.control_file = str(Path(args.control_file).resolve())
        namespace = common_namespace(args, species, False)
        try:
            run_SICER.main(namespace)
        finally:
            # SICER changes into its temporary directory and does not restore cwd.
            os.chdir(start_directory)
        write_call_outputs(args)
    else:
        args.treatment_file = [str(Path(path).resolve()) for path in args.treatment_file]
        if args.control_file:
            args.control_file = [str(Path(path).resolve()) for path in args.control_file]
        namespace = common_namespace(args, species, True)
        try:
            run_SICER_df.main(namespace)
        finally:
            os.chdir(start_directory)
        write_diff_outputs(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise
