#!/usr/bin/env python3
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/utils"))

from samples import load_samples


class WorkflowInputTests(unittest.TestCase):
    def test_shared_control_paths_can_be_deduplicated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("ip1.bam", "ip2.bam", "igg_a.bam", "igg_b.bam"):
                (root / name).touch()
            table = root / "samples.txt"
            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts1\t{root/'ip1.bam'}\t{root/'igg_a.bam'}\tSE\n"
                f"g\ts2\t{root/'ip2.bam'}\t{root/'igg_a.bam'}\tSE\n",
                encoding="utf-8",
            )
            samples = load_samples(table)
            self.assertEqual(samples["control_bam"].drop_duplicates().size, 1)

            table.write_text(
                "group\tsample\ttreatment_bam\tcontrol_bam\tlayout\n"
                f"g\ts1\t{root/'ip1.bam'}\t{root/'igg_a.bam'}\tSE\n"
                f"g\ts2\t{root/'ip2.bam'}\t{root/'igg_b.bam'}\tSE\n",
                encoding="utf-8",
            )
            samples = load_samples(table)
            self.assertEqual(samples["control_bam"].drop_duplicates().size, 2)

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

if __name__ == "__main__":
    unittest.main()
