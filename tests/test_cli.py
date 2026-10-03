import csv
import json
import os
import tempfile
import unittest

from digest_map.__main__ import main


TASK = {
    "dna": [
        {"id": "lin", "sequence": "AAAGATATCTTT", "topology": "linear"},
        {"id": "circ", "sequence": "AAGGCCAAGGCC", "topology": "circular"},
    ],
    "enzymes": [
        {"id": "ecorv", "recognition": "GATATC"},
        {"id": "haeiii", "recognition": "GGCC"},
    ],
    "digests": [
        {"id": "d1", "dna": "lin", "enzymes": ["ecorv"]},
        {"id": "d2", "dna": "circ", "enzymes": ["haeiii"]},
        {"id": "d3", "dna": "circ", "enzymes": []},
    ],
}


class CliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.input_path = os.path.join(self.tmp, "task.json")
        with open(self.input_path, "w", encoding="utf-8") as fh:
            json.dump(TASK, fh)
        self.out = os.path.join(self.tmp, "out")

    def run_cli(self):
        return main([self.input_path, "--out", self.out])

    def test_successful_run_outputs(self):
        self.assertEqual(self.run_cli(), 0)
        with open(os.path.join(self.out, "summary.json"), encoding="utf-8") as fh:
            summary = json.load(fh)
        self.assertEqual(len(summary["digests"]), 3)
        d1 = summary["digests"][0]
        self.assertEqual(d1["cut_sites"][0]["position"], 6)
        self.assertEqual(d1["length_histogram"], {"6": 2})
        d2 = summary["digests"][1]
        self.assertEqual(d2["fragment_count"], 2)
        self.assertEqual(d2["total_fragment_length"], 12)
        self.assertEqual(d2["cut_sites"][0]["position"], 4)
        self.assertEqual(d2["cut_sites"][1]["position"], 10)
        d3 = summary["digests"][2]
        self.assertEqual(d3["fragment_count"], 1)
        self.assertEqual(d3["length_histogram"], {"12": 1})
        with open(os.path.join(self.out, "fragments_d1.csv"),
                  encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0],
                         ["fragment_index", "start", "length", "sequence",
                          "wraps_origin"])
        self.assertEqual(rows[1], ["0", "0", "6", "AAAGAT", "no"])
        self.assertEqual(rows[2], ["1", "6", "6", "ATCTTT", "no"])

    def test_invalid_input_fails_without_output(self):
        with open(self.input_path, "w", encoding="utf-8") as fh:
            json.dump({"dna": [], "enzymes": [], "digests": [
                {"id": "x", "dna": "ghost", "enzymes": []}]}, fh)
        self.assertEqual(self.run_cli(), 1)
        self.assertFalse(os.path.exists(self.out))

    def test_missing_input_file_fails(self):
        self.assertEqual(main([os.path.join(self.tmp, "nope.json"),
                               "--out", self.out]), 1)
        self.assertFalse(os.path.exists(self.out))


if __name__ == "__main__":
    unittest.main()
