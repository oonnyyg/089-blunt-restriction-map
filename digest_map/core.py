"""Core validation and digest logic (standard library only)."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path


class DigestError(ValueError):
    """Raised for any invalid task input."""


_BASES = set("ACGT")
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass
class Dna:
    id: str
    topology: str  # "linear" or "circular"
    sequence: str  # upper-case ACGT


@dataclass
class Enzyme:
    id: str
    recognition: str  # upper-case, even length, self reverse-complement

    @property
    def half(self) -> int:
        return len(self.recognition) // 2


@dataclass
class Digest:
    id: str
    dna_id: str
    enzyme_ids: list[str]


@dataclass
class Task:
    dnas: dict[str, Dna] = field(default_factory=dict)
    enzymes: dict[str, Enzyme] = field(default_factory=dict)
    digests: list[Digest] = field(default_factory=list)


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise DigestError(msg)


def _check_id(value: object, kind: str) -> str:
    _require(isinstance(value, str) and value.strip() != "",
             f"{kind} id must be a non-empty string")
    return value


def _check_sequence(value: object, kind: str) -> str:
    _require(isinstance(value, str) and value != "",
             f"{kind} sequence must be a non-empty string")
    seq = value.upper()
    _require(all(b in _BASES for b in seq),
             f"{kind} sequence must contain only A, C, G, T (case-insensitive)")
    return seq


def _parse_dna(item: object) -> Dna:
    _require(isinstance(item, dict), "each dna entry must be an object")
    dna_id = _check_id(item.get("id"), "dna")
    topology = item.get("topology")
    _require(topology in ("linear", "circular"),
             f"dna {dna_id!r}: topology must be 'linear' or 'circular'")
    seq = _check_sequence(item.get("sequence"), f"dna {dna_id!r}")
    return Dna(id=dna_id, topology=topology, sequence=seq)


def _reverse_complement(seq: str) -> str:
    return "".join(_COMPLEMENT[b] for b in reversed(seq))


def _parse_enzyme(item: object) -> Enzyme:
    _require(isinstance(item, dict), "each enzyme entry must be an object")
    enz_id = _check_id(item.get("id"), "enzyme")
    rec = _check_sequence(item.get("recognition"), f"enzyme {enz_id!r}")
    _require(len(rec) % 2 == 0,
             f"enzyme {enz_id!r}: recognition length must be a positive even number")
    _require(rec == _reverse_complement(rec),
             f"enzyme {enz_id!r}: recognition must equal its own reverse complement")
    return Enzyme(id=enz_id, recognition=rec)


def _parse_digest(item: object, task: Task) -> Digest:
    _require(isinstance(item, dict), "each digest entry must be an object")
    dig_id = _check_id(item.get("id"), "digest")
    dna_id = _check_id(item.get("dna"), f"digest {dig_id!r} dna reference")
    _require(dna_id in task.dnas,
             f"digest {dig_id!r}: unknown dna {dna_id!r}")
    enzyme_ids = item.get("enzymes", [])
    _require(isinstance(enzyme_ids, list),
             f"digest {dig_id!r}: enzymes must be a list")
    seen: set[str] = set()
    for eid in enzyme_ids:
        eid = _check_id(eid, f"digest {dig_id!r} enzyme reference")
        _require(eid in task.enzymes,
                 f"digest {dig_id!r}: unknown enzyme {eid!r}")
        _require(eid not in seen,
                 f"digest {dig_id!r}: duplicate enzyme {eid!r}")
        seen.add(eid)
    return Digest(id=dig_id, dna_id=dna_id, enzyme_ids=list(enzyme_ids))


def load_task(path: str | Path) -> Task:
    """Load and fully validate a task JSON file. Raises DigestError on any problem."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DigestError(f"cannot read task file: {exc}") from exc
    _require(isinstance(data, dict), "task root must be an object")

    task = Task()
    for key, parser, store in (("dnas", _parse_dna, task.dnas),
                               ("enzymes", _parse_enzyme, task.enzymes)):
        entries = data.get(key, [])
        _require(isinstance(entries, list), f"{key} must be a list")
        for entry in entries:
            obj = parser(entry)
            _require(obj.id not in store, f"duplicate {key} id {obj.id!r}")
            store[obj.id] = obj

    digests = data.get("digests", [])
    _require(isinstance(digests, list), "digests must be a list")
    seen_ids: set[str] = set()
    for entry in digests:
        dig = _parse_digest(entry, task)
        _require(dig.id not in seen_ids, f"duplicate digest id {dig.id!r}")
        seen_ids.add(dig.id)
        task.digests.append(dig)
    return task


def find_matches(dna: Dna, enzyme: Enzyme) -> list[int]:
    """Return all match start positions (0-based) of the enzyme on the DNA.

    Linear DNA: no matching across the ends. Circular DNA: matches may wrap
    the origin, but if the recognition site is longer than the DNA there are
    no matches at all.
    """
    rec = enzyme.recognition
    n = len(dna.sequence)
    m = len(rec)
    if m > n:
        return []
    haystack = dna.sequence if dna.topology == "linear" else dna.sequence * 2
    limit = n - m + 1 if dna.topology == "linear" else n
    starts: list[int] = []
    pos = 0
    while True:
        idx = haystack.find(rec, pos)
        if idx < 0 or idx >= limit:
            break
        starts.append(idx)
        pos = idx + 1  # allow overlapping matches
    return starts


@dataclass
class CutSite:
    position: int  # between-base coordinate, 0-based; 0..n
    hits: list[tuple[str, int]]  # (enzyme_id, match start) contributing here


@dataclass
class Fragment:
    start: int
    length: int
    sequence: str
    wraps_origin: bool


def cut_sites(dna: Dna, enzymes: list[Enzyme]) -> list[CutSite]:
    """Merge cuts at identical positions; keep every contributing hit."""
    n = len(dna.sequence)
    by_pos: dict[int, list[tuple[str, int]]] = {}
    for enz in enzymes:
        for start in find_matches(dna, enz):
            pos = (start + enz.half) % n if dna.topology == "circular" else start + enz.half
            by_pos.setdefault(pos, []).append((enz.id, start))
    sites = [CutSite(position=p, hits=sorted(by_pos[p])) for p in sorted(by_pos)]
    return sites


def build_fragments(dna: Dna, sites: list[CutSite]) -> list[Fragment]:
    """Build fragments from sorted unique cut sites."""
    n = len(dna.sequence)
    seq = dna.sequence
    if not sites:
        return [Fragment(start=0, length=n, sequence=seq, wraps_origin=False)]
    positions = [s.position for s in sites]
    fragments: list[Fragment] = []
    if dna.topology == "linear":
        bounds = [0] + positions + [n]
        for a, b in zip(bounds, bounds[1:]):
            fragments.append(Fragment(start=a, length=b - a,
                                      sequence=seq[a:b], wraps_origin=False))
    else:
        count = len(positions)
        for i, a in enumerate(positions):
            b = positions[(i + 1) % count]
            length = (b - a) % n
            if length == 0:
                length = n  # single cut site: full-length linear fragment
            doubled = seq + seq
            fragments.append(Fragment(start=a, length=length,
                                      sequence=doubled[a:a + length],
                                      wraps_origin=a + length > n))
    return fragments


def digest_one(dna: Dna, enzymes: list[Enzyme]) -> tuple[list[CutSite], list[Fragment]]:
    sites = cut_sites(dna, enzymes)
    fragments = build_fragments(dna, sites)
    total = sum(f.length for f in fragments)
    _require(total == len(dna.sequence),
             f"internal error: fragment lengths sum to {total}, expected {len(dna.sequence)}")
    return sites, fragments


def run_task(task: Task, out_dir: str | Path) -> dict:
    """Run all digests, write per-digest fragment CSVs and a summary JSON.

    Returns the summary dict. Writes nothing if any digest fails.
    """
    results: list[dict] = []
    csv_rows: list[tuple[str, list[list[object]]]] = []
    for dig in task.digests:
        dna = task.dnas[dig.dna_id]
        enzymes = [task.enzymes[eid] for eid in dig.enzyme_ids]
        sites, fragments = digest_one(dna, enzymes)
        rows = [[f.start, f.length, f.sequence, str(f.wraps_origin).lower()]
                for f in fragments]
        csv_rows.append((dig.id, rows))
        length_histogram: dict[str, int] = {}
        for f in fragments:
            length_histogram[str(f.length)] = length_histogram.get(str(f.length), 0) + 1
        results.append({
            "id": dig.id,
            "dna": dna.id,
            "topology": dna.topology,
            "dna_length": len(dna.sequence),
            "enzymes": list(dig.enzyme_ids),
            "cut_sites": [
                {"position": s.position,
                 "hits": [{"enzyme": e, "match_start": ms} for e, ms in s.hits]}
                for s in sites
            ],
            "fragment_count": len(fragments),
            "fragments_csv": f"{dig.id}.fragments.csv",
            "length_histogram": dict(sorted(length_histogram.items(), key=lambda kv: int(kv[0]))),
        })

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for dig_id, rows in csv_rows:
        with (out / f"{dig_id}.fragments.csv").open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["start", "length", "sequence", "wraps_origin"])
            writer.writerows(rows)
    summary = {"digests": results}
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary
