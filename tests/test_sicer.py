import unittest

import numpy as np
from sicer.lib.associate_tags_with_regions import tag_position
from sicer.src.remove_redundant_reads import remove_redundant_1chrom_single_strand_sorted


class SicerBehaviorTests(unittest.TestCase):
    def test_se_shift_and_native_coordinate_redundancy(self):
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
