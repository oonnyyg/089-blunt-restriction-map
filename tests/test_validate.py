import json
import os
import tempfile
import unittest

from digest_map.errors import InputError
from digest_map.validate import load_task


def write_task(tmp, data):
    path = os.path.join(tmp, "task.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    return path


VALID = {
    "dna": [{"id": "d1", "sequence": "aaagatatcttt", "topology": "linear"}],
    "enzymes": [{"id": "eco", "recognition": "gatatc"}],
    "digests": [{"id": "x", "dna": "d1", "enzymes": ["eco"]}],
}


class ValidateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def load(self, data):
        return load_task(write_task(self.tmp, data))

    def test_valid_task_uppercases(self):
        dnas, enzymes, digests = self.load(VALID)
        self.assertEqual(dnas[0].sequence, "AAAGATATCTTT")
        self.assertEqual(enzymes[0].recognition, "GATATC")
        self.assertEqual(digests[0].enzyme_ids, ("eco",))

    def test_empty_enzyme_set_allowed(self):
        data = dict(VALID, digests=[{"id": "x", "dna": "d1", "enzymes": []}])
        _, _, digests = self.load(data)
        self.assertEqual(digests[0].enzyme_ids, ())

    def test_empty_sequence_rejected(self):
        data = dict(VALID, dna=[{"id": "d1", "sequence": "", "topology": "linear"}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_invalid_base_rejected(self):
        data = dict(VALID, dna=[{"id": "d1", "sequence": "ACGTN", "topology": "linear"}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_bad_topology_rejected(self):
        data = dict(VALID, dna=[{"id": "d1", "sequence": "ACGT", "topology": "x"}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_odd_recognition_rejected(self):
        data = dict(VALID, enzymes=[{"id": "e", "recognition": "ATATA"}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_non_palindromic_recognition_rejected(self):
        data = dict(VALID, enzymes=[{"id": "e", "recognition": "GATACC"}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_duplicate_ids_rejected(self):
        data = dict(VALID, dna=VALID["dna"] + [dict(VALID["dna"][0])])
        with self.assertRaises(InputError):
            self.load(data)

    def test_duplicate_enzymes_in_digest_rejected(self):
        data = dict(VALID, digests=[{"id": "x", "dna": "d1",
                                     "enzymes": ["eco", "eco"]}])
        with self.assertRaises(InputError):
            self.load(data)

    def test_unknown_references_rejected(self):
        data = dict(VALID, digests=[{"id": "x", "dna": "nope", "enzymes": []}])
        with self.assertRaises(InputError):
            self.load(data)
        data = dict(VALID, digests=[{"id": "x", "dna": "d1", "enzymes": ["nope"]}])
        with self.assertRaises(InputError):
            self.load(data)


if __name__ == "__main__":
    unittest.main()
