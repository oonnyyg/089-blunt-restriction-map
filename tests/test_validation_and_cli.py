import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from digest_map.core import DigestError, load_task


def write_task(tmp: Path, data: dict) -> Path:
    path = tmp / "task.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


VALID_TASK = {
    "dnas": [
        {"id": "lin", "topology": "linear", "sequence": "aaGATATCgg"},
        {"id": "circ", "topology": "circular", "sequence": "GATCCATC"},
    ],
    "enzymes": [{"id": "EcoRV", "recognition": "gatatc"}],
    "digests": [
        {"id": "d1", "dna": "lin", "enzymes": ["EcoRV"]},
        {"id": "d2", "dna": "circ", "enzymes": ["EcoRV"]},
        {"id": "d3", "dna": "lin", "enzymes": []},
    ],
}


class TestValidation(unittest.TestCase):
    def check_invalid(self, mutate, msg_part):
        with tempfile.TemporaryDirectory() as td:
            data = json.loads(json.dumps(VALID_TASK))
            mutate(data)
            path = write_task(Path(td), data)
            with self.assertRaises(DigestError) as ctx:
                load_task(path)
            self.assertIn(msg_part, str(ctx.exception))

    def test_valid_task_loads(self):
        with tempfile.TemporaryDirectory() as td:
            task = load_task(write_task(Path(td), VALID_TASK))
            self.assertEqual(len(task.digests), 3)
            self.assertEqual(task.dnas["lin"].sequence, "AAGATATCGG")
            self.assertEqual(task.enzymes["EcoRV"].recognition, "GATATC")

    def test_empty_sequence_rejected(self):
        self.check_invalid(lambda d: d["dnas"][0].update(sequence=""), "non-empty")

    def test_invalid_base_rejected(self):
        self.check_invalid(lambda d: d["dnas"][0].update(sequence="ACGN"), "only A, C, G, T")

    def test_bad_topology_rejected(self):
        self.check_invalid(lambda d: d["dnas"][0].update(topology="weird"), "topology")

    def test_odd_recognition_rejected(self):
        self.check_invalid(lambda d: d["enzymes"][0].update(recognition="ATA"), "even")

    def test_non_self_reverse_complement_rejected(self):
        self.check_invalid(lambda d: d["enzymes"][0].update(recognition="AAAA"), "reverse complement")

    def test_duplicate_ids_rejected(self):
        self.check_invalid(lambda d: d["dnas"].append(dict(d["dnas"][0])), "duplicate")

    def test_unknown_dna_reference_rejected(self):
        self.check_invalid(lambda d: d["digests"][0].update(dna="nope"), "unknown dna")

    def test_unknown_enzyme_reference_rejected(self):
        self.check_invalid(lambda d: d["digests"][0].update(enzymes=["nope"]), "unknown enzyme")

    def test_duplicate_enzyme_in_digest_rejected(self):
        self.check_invalid(lambda d: d["digests"][0].update(enzymes=["EcoRV", "EcoRV"]),
                           "duplicate enzyme")


class TestCli(unittest.TestCase):
    def run_cli(self, task_path, out_dir):
        return subprocess.run(
            [sys.executable, "-m", "digest_map", str(task_path), "--out", str(out_dir)],
            capture_output=True, text=True,
            cwd=Path(__file__).resolve().parent.parent)

    def test_cli_writes_csv_and_summary(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            task_path = write_task(tmp, VALID_TASK)
            out = tmp / "out"
            proc = self.run_cli(task_path, out)
            self.assertEqual(proc.returncode, 0, proc.stderr)

            summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(len(summary["digests"]), 3)
            d1 = summary["digests"][0]
            self.assertEqual(d1["cut_sites"], [
                {"position": 5, "hits": [{"enzyme": "EcoRV", "match_start": 2}]}])
            self.assertEqual(d1["length_histogram"], {"5": 2})
            d2 = summary["digests"][1]
            self.assertEqual(d2["fragment_count"], 1)  # single cut -> full length
            d3 = summary["digests"][2]
            self.assertEqual(d3["cut_sites"], [])
            self.assertEqual(d3["length_histogram"], {"10": 1})

            with (out / "d1.fragments.csv").open(encoding="utf-8", newline="") as fh:
                rows = list(csv.reader(fh))
            self.assertEqual(rows[0], ["start", "length", "sequence", "wraps_origin"])
            self.assertEqual(rows[1], ["0", "5", "AAGAT", "false"])
            self.assertEqual(rows[2], ["5", "5", "ATCGG", "false"])

    def test_cli_invalid_input_fails_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            bad = json.loads(json.dumps(VALID_TASK))
            bad["dnas"][0]["sequence"] = "ACGX"
            task_path = write_task(tmp, bad)
            out = tmp / "out"
            proc = self.run_cli(task_path, out)
            self.assertEqual(proc.returncode, 1)
            self.assertIn("error:", proc.stderr)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
