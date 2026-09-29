#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


def read_chromosomes(path):
    chroms = []
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[1].isdigit() or int(fields[1]) <= 0:
                raise ValueError(f"Invalid chrom sizes line {line_number}: {line.rstrip()}")
            chroms.append(fields[0])
    if not chroms or len(chroms) != len(set(chroms)):
        raise ValueError("Chromosome list is empty or contains duplicates")
    return chroms


def run_checked(command, **kwargs):
    result = subprocess.run(command, **kwargs)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    return result


def validate_bam(samtools, bam):
    run_checked([samtools, "quickcheck", "-v", bam])
    result = run_checked(
        [samtools, "view", "-c", "-f", "1", bam],
        text=True,
        stdout=subprocess.PIPE,
    )
    paired = int(result.stdout.strip() or "0")
    if paired:
        raise ValueError(
            f"Paired reads detected in {bam} ({paired} alignments with flag 0x1); "
            "SICER2_workflow v1 accepts SE BAM only"
        )


def stream_bed(samtools, bedtools, bams):
    if len(bams) == 1:
        bed = subprocess.Popen(
            [bedtools, "bamtobed", "-i", bams[0]],
            stdout=subprocess.PIPE,
            text=True,
        )
        return None, bed
    merge = subprocess.Popen(
        [samtools, "merge", "-u", "-", *bams],
        stdout=subprocess.PIPE,
    )
    bed = subprocess.Popen(
        [bedtools, "bamtobed", "-i", "stdin"],
        stdin=merge.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
    merge.stdout.close()
    return merge, bed


def main():
    parser = argparse.ArgumentParser(description="Pool SE BAMs into a nuclear BED6 file")
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--bedtools", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--bam", nargs="+", required=True)
    args = parser.parse_args()

    chromosomes = set(read_chromosomes(args.chrom_sizes))
    for bam in args.bam:
        validate_bam(args.samtools, bam)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.stats).parent.mkdir(parents=True, exist_ok=True)
    merge, bed = stream_bed(args.samtools, args.bedtools, args.bam)
    kept = 0
    excluded = 0
    with open(args.output, "w", encoding="utf-8") as output:
        assert bed.stdout is not None
        for line in bed.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                bed.kill()
                raise ValueError(f"bedtools produced fewer than six BED columns: {line.rstrip()}")
            if fields[0] in chromosomes:
                output.write("\t".join(fields[:6]) + "\n")
                kept += 1
            else:
                excluded += 1
    bed_return = bed.wait()
    merge_return = merge.wait() if merge is not None else 0
    if bed_return or merge_return:
        raise RuntimeError(
            f"BAM pooling failed: samtools merge={merge_return}, bedtools bamtobed={bed_return}"
        )
    if kept == 0:
        raise ValueError("No reads remain on chromosomes listed in --chrom-sizes")

    with open(args.stats, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_file_count\t{len(args.bam)}\n")
        handle.write(f"nuclear_bed_reads\t{kept}\n")
        handle.write(f"excluded_contig_reads\t{excluded}\n")
        handle.write(f"input_files\t{';'.join(args.bam)}\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise
