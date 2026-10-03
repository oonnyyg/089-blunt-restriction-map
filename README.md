# blunt-restriction-map

Predict blunt-end restriction fragments for linear and circular DNA.
Pure Python standard library; no sequence databases are accessed.

## Usage

```
python -m digest_map input.json --out output_dir
```

The output directory receives one `fragments_<digest_id>.csv` per digest
plus a `summary.json` covering all digests. Any invalid input fails the
whole batch and nothing is written.

## Input format

```json
{
  "dna": [
    {"id": "plasmid", "sequence": "AAGGCCAAGGCC", "topology": "circular"},
    {"id": "chr1", "sequence": "AAAGATATCTTT", "topology": "linear"}
  ],
  "enzymes": [
    {"id": "EcoRV", "recognition": "GATATC"},
    {"id": "HaeIII", "recognition": "GGCC"}
  ],
  "digests": [
    {"id": "d1", "dna": "chr1", "enzymes": ["EcoRV"]},
    {"id": "d2", "dna": "plasmid", "enzymes": ["HaeIII", "EcoRV"]},
    {"id": "d3", "dna": "plasmid", "enzymes": []}
  ]
}
```
    {"id": "d3", "dna": "plasmid", "enzymes": []}
  ]
}
```

Validation rules:

- `sequence` is non-empty, A/C/G/T only, case-insensitive; `topology` is
  `linear` or `circular`.
- `recognition` has positive even length, A/C/G/T only, and must equal its
  own reverse complement (blunt cut at the midpoint; no methylation or
  partial-digest modeling).
- Ids are unique within each of `dna`, `enzymes`, `digests`.
- Each digest references one existing DNA and a duplicate-free enzyme set
  (the empty set is allowed).

## Coordinate model

- Coordinates are zero-based positions *between* bases: position `p` lies
  immediately before base `p`, and base `i` occupies interval `[i, i+1)`.
- A recognition match starting at `s` with length `k` cuts at
  `s + k/2`. On circular DNA the cut position is taken modulo the DNA
  length `n`.
- Linear DNA: matches must lie fully inside the molecule (no end
  spanning). Circular DNA: matches may wrap across the origin; if the
  recognition sequence is longer than the DNA, that enzyme has no matches.
  Overlapping matches are all reported.
- Multiple enzymes cutting at the same position produce a single cut; the
  summary keeps every hit enzyme and its match start. No zero-length
  fragments are created.
- Linear fragments span consecutive cut positions including the molecule
  ends. Circular fragments run from each cut to the next cut, the last one
  wrapping across the origin. A circular molecule with exactly one cut
  yields one full-length linear fragment; with no cuts it is reported as
  the intact circle (start 0, length n, `wraps_origin = no`).
- Fragment `start` is the coordinate of its first base. `wraps_origin` is
  `yes` only when `start + length > n`, i.e. the fragment sequence extends
  past the original end of the sequence. Fragment sequences are always
  reported 5'->3' on the forward strand.
- Fragment lengths always sum to the DNA length.

## Output

`fragments_<digest_id>.csv` columns:
`fragment_index, start, length, sequence, wraps_origin`.

`summary.json` contains, per digest: cut sites (position, enzymes,
match_starts), fragment count, total fragment length, a length histogram
(fragment count per length), and the CSV file name.

## Tests

```
python -m unittest discover -s tests
```
