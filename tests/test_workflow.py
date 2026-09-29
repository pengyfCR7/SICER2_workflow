#!/usr/bin/env python3
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from sicer.lib.associate_tags_with_regions import tag_position
from sicer.src.remove_redundant_reads import remove_redundant_1chrom_single_strand_sorted

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/lib"))

from workflow_io import load_samples


class WorkflowInputTests(unittest.TestCase):
    def test_control_identity_dedup_uses_files_not_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("ip1.bam", "ip2.bam", "igg_a.bam", "igg_b.bam"):
                (root / name).touch()
            os.symlink(root / "igg_a.bam", root / "igg_a_link.bam")
            os.link(root / "igg_a.bam", root / "igg_a_hardlink.bam")
            table = root / "samples.txt"
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
            table = root / "samples.txt"
            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts\t{root/'ip.bam'}\t\tPE\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "only layout=SE"):
                load_samples(table)

    def test_sicer_se_shift_and_native_coordinate_redundancy(self):
        self.assertEqual(tag_position(("Chr1", 100, 150, "r1", 60, "+"), 150), 175)
        self.assertEqual(tag_position(("Chr1", 100, 150, "r2", 60, "-"), 150), 74)
        reads = np.array(
            [("Chr1", 100, 150), ("Chr1", 100, 150)],
            dtype=[("chrom", "U6"), ("start", np.int32), ("end", np.int32)],
        )
        total, retained, _ = remove_redundant_1chrom_single_strand_sorted(reads, 1)
        self.assertEqual((total, retained), (2, 1))

if __name__ == "__main__":
    unittest.main()
