"""Analytic ambient propagation tests independent of telescope geometry."""

import math
import pathlib
import unittest

from obdeect.ray_tracing.propagation import ambient_group_index


class TestPropagation(unittest.TestCase):
    def test_exponential_refractivity_matches_reference_lookup(self):
        contents = "altitude refractive_index\n" + "\n".join(
            f"{height} {0.0003 * math.exp(-height / 8):.17g}" for height in (0, 1, 2, 3)
        )
        for altitude in (0, 123, 1000, 1789):
            self.assertAlmostEqual(
                ambient_group_index(contents, altitude),
                1 + 0.0003 * math.exp(-altitude / 8000),
                places=15,
            )
        self.assertEqual(ambient_group_index(contents, 3000), 1)
        with self.assertRaisesRegex(ValueError, "upper layer"):
            ambient_group_index(contents, 2500)

    def test_invalid_profile_rejects(self):
        for contents in (
            "altitude index\n0 1\n1 2",
            "altitude refractive_index\n0 0\n1 0.1",
            "altitude refractive_index\n0 0.1\n0 0.2",
        ):
            with self.assertRaises(ValueError):
                ambient_group_index(contents, 100)


# The C++ implementation in propagation (timing.hpp, ambient_group_index) is tested against the
# same golden values in src/propagation/cpp/tests/test_timing.cpp, on the same table. Ray
# tracing does not depend on propagation, so a change on one side fails the test of the other.
PROFILE = (
    pathlib.Path(__file__).resolve().parents[3]
    / "propagation/cpp/tests/data/atmprof_ecmwf_north_winter_fixed.dat"
)
GOLDEN = {
    1000.0: 1.000255944992766,
    2180.0: 1.000225191893913,
    4500.0: 1.000177406919158,
    60000.0: 1.0000000583180002,
}


@unittest.skipUnless(PROFILE.exists(), "propagation test data not in this source tree")
class TestAgreesWithPropagation(unittest.TestCase):
    def test_golden_values_shared_with_propagation_cpp_test(self):
        rows = [
            line.split()
            for line in PROFILE.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        contents = "altitude refractive_index\n" + "\n".join(f"{r[0]} {r[3]}" for r in rows)
        for altitude, expected in GOLDEN.items():
            self.assertAlmostEqual(ambient_group_index(contents, altitude), expected, delta=1e-12)
