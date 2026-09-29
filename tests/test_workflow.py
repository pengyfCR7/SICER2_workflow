#!/usr/bin/env python3
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from sicer.lib.associate_tags_with_regions import tag_position
from sicer.src.remove_redundant_reads import remove_redundant_1chrom_single_strand_sorted

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/lib"))
sys.path.insert(0, str(ROOT / "workflow/scripts"))

from workflow_io import (
    load_samples,
    validate_parameters,
    validate_tair10_nuclear_chrom_sizes,
)
from run_sicer2 import write_diff_outputs


class WorkflowInputTests(unittest.TestCase):
    def test_control_identity_dedup_uses_files_not_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("ip1.bam", "ip2.bam", "igg_a.bam", "igg_b.bam"):
                (root / name).touch()
            os.symlink(root / "igg_a.bam", root / "igg_a_link.bam")
            os.link(root / "igg_a.bam", root / "igg_a_hardlink.bam")
            table = root / "samples.tsv"
            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts1\t{root/'ip1.bam'}\t{root/'igg_a.bam'}\tSE\n"
                f"g\ts2\t{root/'ip2.bam'}\t{root/'igg_a_link.bam'}\tSE\n",
                encoding="utf-8",
            )
            info = load_samples(table)
            self.assertEqual(len(info["group_controls"]["g"]), 1)

            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts1\t{root/'ip1.bam'}\t{root/'igg_a.bam'}\tSE\n"
                f"g\ts2\t{root/'ip2.bam'}\t{root/'igg_a_hardlink.bam'}\tSE\n",
                encoding="utf-8",
            )
            info = load_samples(table)
            self.assertEqual(len(info["group_controls"]["g"]), 1)

            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts1\t{root/'ip1.bam'}\t{root/'igg_a.bam'}\tSE\n"
                f"g\ts2\t{root/'ip2.bam'}\t{root/'igg_b.bam'}\tSE\n",
                encoding="utf-8",
            )
            info = load_samples(table)
            self.assertEqual(len(info["group_controls"]["g"]), 2)

    def test_pe_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "ip.bam").touch()
            table = root / "samples.tsv"
            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts\t{root/'ip.bam'}\t\tPE\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "only layout=SE"):
                load_samples(table)

    def test_fdr_names_and_parameter_ranges(self):
        call = {
            "redundancy_threshold": 1,
            "window_size": 200,
            "fragment_size": 150,
            "effective_genome_fraction": 0.8,
            "gap_size": 600,
            "false_discovery_rate": 0.01,
        }
        diff = {"false_discovery_rate_df": 0.01, "min_fold_change": 1.3}
        validate_parameters(call, diff)
        broken = dict(diff, false_discovery_rate_df=0)
        with self.assertRaisesRegex(ValueError, "false_discovery_rate_df"):
            validate_parameters(call, broken)

    def test_tair10_nuclear_chrom_sizes(self):
        validate_tair10_nuclear_chrom_sizes(ROOT / "config/TAIR10.nuclear.chrom.sizes")

    def test_sicer_se_shift_and_native_coordinate_redundancy(self):
        self.assertEqual(tag_position(("Chr1", 100, 150, "r1", 60, "+"), 150), 175)
        self.assertEqual(tag_position(("Chr1", 100, 150, "r2", 60, "-"), 150), 74)
        reads = np.array(
            [("Chr1", 100, 150), ("Chr1", 100, 150)],
            dtype=[("chrom", "U6"), ("start", np.int32), ("end", np.int32)],
        )
        total, retained, _ = remove_redundant_1chrom_single_strand_sorted(reads, 1)
        self.assertEqual((total, retained), (2, 1))

    def test_differential_direction_and_fold_change_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            source = raw / "test-and-reference-W200-G600-summary"
            source.write_text(
                "#chrom\tstart\tend\tReadcount_A\tNormalized_Readcount_A\tReadcountB\t"
                "Normalized_Readcount_B\tFc_A_vs_B\tpvalue_A_vs_B\tFDR_A_vs_B\t"
                "Fc_B_vs_A\tpvalue_B_vs_A\tFDR_B_vs_A\n"
                "Chr1\t0\t200\t20\t2\t10\t1\t1.5\t0.001\t0.005\t0.6667\t1\t1\n"
                "Chr1\t400\t600\t5\t0.5\t20\t2\t0.7\t1\t1\t1.4286\t0.001\t0.005\n",
                encoding="utf-8",
            )
            args = SimpleNamespace(
                output_directory=str(raw), window_size=200, gap_size=600,
                false_discovery_rate_df=0.01, min_fold_change=1.3,
                all_islands=str(root / "all.tsv"),
                increased_fdr=str(root / "increased.fdr.tsv"),
                decreased_fdr=str(root / "decreased.fdr.tsv"),
                increased_filtered=str(root / "increased.filtered.tsv"),
                decreased_filtered=str(root / "decreased.filtered.tsv"),
            )
            write_diff_outputs(args)
            data_rows = lambda path: len(Path(path).read_text(encoding="utf-8").splitlines()) - 1
            self.assertEqual(data_rows(args.increased_fdr), 1)
            self.assertEqual(data_rows(args.decreased_fdr), 1)
            self.assertEqual(data_rows(args.increased_filtered), 1)
            self.assertEqual(data_rows(args.decreased_filtered), 1)


if __name__ == "__main__":
    unittest.main()
