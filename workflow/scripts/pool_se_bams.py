#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Pool SE BAMs into a nuclear BED6 file")
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--bedtools", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bam", nargs="+", required=True)
    args = parser.parse_args()

    chromosomes = {
        line.split("\t", 1)[0]
        for line in Path(args.chrom_sizes).read_text().splitlines()
        if line.strip()
    }

    for bam in args.bam:
        subprocess.run([args.samtools, "quickcheck", "-v", bam], check=True)
        paired = subprocess.run(
            [args.samtools, "view", "-c", "-f", "1", bam],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
        if int(paired or 0):
            raise ValueError(f"Paired reads detected in {bam}; only SE BAM is supported")

    if len(args.bam) == 1:
        merge = None
        bed = subprocess.Popen(
            [args.bedtools, "bamtobed", "-i", args.bam[0]],
            stdout=subprocess.PIPE,
            text=True,
        )
    else:
        merge = subprocess.Popen(
            [args.samtools, "merge", "-u", "-", *args.bam],
            stdout=subprocess.PIPE,
        )
        bed = subprocess.Popen(
            [args.bedtools, "bamtobed", "-i", "stdin"],
            stdin=merge.stdout,
            stdout=subprocess.PIPE,
            text=True,
        )
        merge.stdout.close()

    kept = 0
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as output:
        for line in bed.stdout:
            fields = line.rstrip().split("\t")
            if fields[0] in chromosomes:
                output.write("\t".join(fields[:6]) + "\n")
                kept += 1

    bed_return = bed.wait()
    merge_return = merge.wait() if merge else 0
    if bed_return or merge_return:
        raise RuntimeError("BAM to BED conversion failed")
    if not kept:
        raise ValueError("No reads remain on chromosomes listed in --chrom-sizes")
    print(f"Nuclear BED reads: {kept}")


if __name__ == "__main__":
    main()
