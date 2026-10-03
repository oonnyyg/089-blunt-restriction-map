import json
import unittest

from digest_map.core import (
    Digest,
    DigestError,
    Dna,
    Enzyme,
    Task,
    build_fragments,
    cut_sites,
    digest_one,
    find_matches,
    load_task,
)


def linear(seq):
    return Dna(id="d", topology="linear", sequence=seq)


def circular(seq):
    return Dna(id="d", topology="circular", sequence=seq)


ECORV = Enzyme(id="EcoRV", recognition="GATATC")  # cuts GAT|ATC


class TestFindMatches(unittest.TestCase):
    def test_linear_basic_and_no_end_wrap(self):
        dna = linear("AAGATATCGG")
        self.assertEqual(find_matches(dna, ECORV), [2])

    def test_overlapping_matches(self):
        # ATAT contains overlapping ATAT at 0 and 2
        enz = Enzyme(id="X", recognition="ATAT")
        self.assertEqual(find_matches(linear("ATATAT"), enz), [0, 2])

    def test_circular_wraps_origin(self):
        # site GATATC split across origin: seq[3:6]="GAT" + seq[0:3]="ATC"
        dna = circular("ATCGAT")
        self.assertEqual(find_matches(dna, ECORV), [3])

    def test_recognition_longer_than_dna_no_match(self):
        self.assertEqual(find_matches(circular("GAT"), ECORV), [])
        self.assertEqual(find_matches(linear("GAT"), ECORV), [])

    def test_case_insensitive_input_normalized(self):
        dna = linear("aagatatcgg".upper())
        self.assertEqual(find_matches(dna, ECORV), [2])


class TestCutSites(unittest.TestCase):
    def test_cut_position_is_midpoint(self):
        dna = linear("AAGATATCGG")
        sites = cut_sites(dna, [ECORV])
        self.assertEqual([s.position for s in sites], [5])
        self.assertEqual(sites[0].hits, [("EcoRV", 2)])

    def test_same_position_cut_once_keeps_all_hits(self):
        enz2 = Enzyme(id="Other", recognition="GATATC")
        dna = linear("AAGATATCGG")
        sites = cut_sites(dna, [ECORV, enz2])
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].hits, [("EcoRV", 2), ("Other", 2)])

    def test_circular_cut_position_modulo(self):
        # match start 3, half 3 -> cut at 6 % 6 = 0
        dna = circular("ATCGAT")
        sites = cut_sites(dna, [ECORV])
        self.assertEqual([s.position for s in sites], [0])


class TestFragmentsLinear(unittest.TestCase):
    def test_linear_includes_end_stubs(self):
        dna = linear("AAGATATCGG")
        sites = cut_sites(dna, [ECORV])
        frags = build_fragments(dna, sites)
        self.assertEqual([(f.start, f.length, f.sequence) for f in frags],
                         [(0, 5, "AAGAT"), (5, 5, "ATCGG")])
        self.assertFalse(any(f.wraps_origin for f in frags))

    def test_no_cuts_keeps_full_molecule(self):
        dna = linear("AAAA")
        frags = build_fragments(dna, [])
        self.assertEqual(len(frags), 1)
        self.assertEqual(frags[0].length, 4)


class TestFragmentsCircular(unittest.TestCase):
    def test_circular_two_cuts(self):
        # cuts at 2 and 6 on length-8 circle
        dna = circular("AACCGGTT")
        sites = cut_sites(dna, [Enzyme(id="A", recognition="ACGT")])
        # no matches; build sites manually instead
        from digest_map.core import CutSite
        sites = [CutSite(2, []), CutSite(6, [])]
        frags = build_fragments(dna, sites)
        self.assertEqual([(f.start, f.length) for f in frags], [(2, 4), (6, 4)])
        self.assertEqual(frags[0].sequence, "CCGG")
        self.assertEqual(frags[1].sequence, "TTAA")
        self.assertTrue(frags[1].wraps_origin)
        self.assertFalse(frags[0].wraps_origin)

    def test_single_cut_full_length_linear_fragment(self):
        dna = circular("ATCGAT")
        sites = cut_sites(dna, [ECORV])  # single cut at 0
        frags = build_fragments(dna, sites)
        self.assertEqual(len(frags), 1)
        self.assertEqual(frags[0].length, 6)
        self.assertEqual(frags[0].sequence, "ATCGAT")

    def test_no_cut_keeps_circular_molecule(self):
        dna = circular("ACGT")
        frags = build_fragments(dna, [])
        self.assertEqual(len(frags), 1)
        self.assertEqual(frags[0].length, 4)
        self.assertFalse(frags[0].wraps_origin)

    def test_no_zero_length_from_multiple_enzymes_same_site(self):
        dna = circular("AAGATATCGG")
        enz2 = Enzyme(id="Other", recognition="GATATC")
        sites, frags = digest_one(dna, [ECORV, enz2])
        self.assertEqual(len(sites), 1)
        self.assertEqual(len(frags), 1)
        self.assertEqual(frags[0].length, 10)

    def test_fragment_lengths_sum_to_dna_length(self):
        dna = circular("GATATCGATATCAA")
        _, frags = digest_one(dna, [ECORV])
        self.assertEqual(sum(f.length for f in frags), len(dna.sequence))
