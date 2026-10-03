import unittest

from digest_map.core import compute_cut_sites, compute_fragments, run_digest
from digest_map.validate import Dna, Enzyme


class MatchingTests(unittest.TestCase):
    def test_linear_basic_cut(self):
        dna = Dna("d", "AAAGATATCTTT", "linear")
        sites = compute_cut_sites(dna, [Enzyme("e", "GATATC")])
        self.assertEqual([s.position for s in sites], [6])
        self.assertEqual(sites[0].match_starts, (3,))

    def test_overlapping_matches_not_skipped(self):
        dna = Dna("d", "AAAAAA", "linear")
        sites = compute_cut_sites(dna, [Enzyme("e", "AAAA")])
        self.assertEqual([s.position for s in sites], [2, 3, 4])
        self.assertEqual(sites[0].match_starts, (0,))
        self.assertEqual(sites[2].match_starts, (2,))

    def test_circular_match_across_origin(self):
        dna = Dna("d", "TCGA", "circular")
        sites = compute_cut_sites(dna, [Enzyme("e", "GATC")])
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].position, (2 + 2) % 4)
        self.assertEqual(sites[0].match_starts, (2,))

    def test_circular_cut_position_modulo(self):
        dna = Dna("d", "ATCG", "circular")
        sites = compute_cut_sites(dna, [Enzyme("e", "GATC")])
        self.assertEqual(sites[0].position, (3 + 2) % 4)
        self.assertEqual(sites[0].match_starts, (3,))

    def test_recognition_longer_than_circular_dna_no_match(self):
        dna = Dna("d", "ACGT", "circular")
        sites = compute_cut_sites(dna, [Enzyme("e", "ACGTAC")])
        self.assertEqual(sites, [])

    def test_same_position_multiple_enzymes_cut_once(self):
        dna = Dna("d", "TTGATATCAA", "linear")
        enzymes = [Enzyme("e1", "GATATC"), Enzyme("e2", "GATATC")]
        sites = compute_cut_sites(dna, enzymes)
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].enzymes, ("e1", "e2"))
        self.assertEqual(sites[0].match_starts, (2, 2))


class FragmentTests(unittest.TestCase):
    def test_linear_fragments_include_ends(self):
        dna = Dna("d", "AAAGATATCTTT", "linear")
        frags = compute_fragments(dna, [6])
        self.assertEqual([(f.start, f.length, f.sequence) for f in frags],
                         [(0, 6, "AAAGAT"), (6, 6, "ATCTTT")])
        self.assertFalse(any(f.wraps_origin for f in frags))

    def test_linear_no_cuts(self):
        dna = Dna("d", "ACGT", "linear")
        frags = compute_fragments(dna, [])
        self.assertEqual(len(frags), 1)
        self.assertEqual((frags[0].start, frags[0].length), (0, 4))

    def test_circular_no_cuts_keeps_circle(self):
        dna = Dna("d", "ACGTAC", "circular")
        frags = compute_fragments(dna, [])
        self.assertEqual(len(frags), 1)
        self.assertEqual((frags[0].start, frags[0].length, frags[0].sequence),
                         (0, 6, "ACGTAC"))
        self.assertFalse(frags[0].wraps_origin)

    def test_circular_single_cut_full_length_linear(self):
        dna = Dna("d", "ACGTAC", "circular")
        frags = compute_fragments(dna, [2])
        self.assertEqual(len(frags), 1)
        self.assertEqual((frags[0].start, frags[0].length), (2, 6))
        self.assertEqual(frags[0].sequence, "GTACAC")
        self.assertTrue(frags[0].wraps_origin)

    def test_circular_multiple_cuts_wrap(self):
        dna = Dna("d", "AACCGGTT", "circular")
        frags = compute_fragments(dna, [2, 6])
        self.assertEqual([(f.start, f.length) for f in frags], [(2, 4), (6, 4)])
        self.assertEqual(frags[0].sequence, "CCGG")
        self.assertEqual(frags[1].sequence, "TTAA")
        self.assertFalse(frags[0].wraps_origin)
        self.assertTrue(frags[1].wraps_origin)

    def test_circular_cut_at_zero_does_not_wrap(self):
        dna = Dna("d", "AACCGGTT", "circular")
        frags = compute_fragments(dna, [0])
        self.assertFalse(frags[0].wraps_origin)
        self.assertEqual(frags[0].sequence, "AACCGGTT")

    def test_total_length_preserved(self):
        dna = Dna("d", "GATATCGATATC", "linear")
        result = run_digest("x", dna, [Enzyme("e", "GATATC")])
        self.assertEqual(sum(f.length for f in result.fragments), 12)
        self.assertEqual(len(result.fragments), 3)


if __name__ == "__main__":
    unittest.main()
