"""Input loading and validation."""

from __future__ import annotations

import json
from dataclasses import dataclass

from .errors import InputError

_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Dna:
    id: str
    sequence: str
    topology: str  # "linear" or "circular"


@dataclass(frozen=True)
class Enzyme:
    id: str
    recognition: str


@dataclass(frozen=True)
class Digest:
    id: str
    dna_id: str
    enzyme_ids: tuple[str, ...]


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _require(cond: bool, message: str) -> None:
    if not cond:
        raise InputError(message)


def _require_str_id(value: object, what: str) -> str:
    _require(isinstance(value, str) and value != "", f"{what} must be a non-empty string")
    return value


def _parse_dna(item: object, index: int) -> Dna:
    _require(isinstance(item, dict), f"dna[{index}] must be an object")
    dna_id = _require_str_id(item.get("id"), f"dna[{index}].id")
    seq = item.get("sequence")
    _require(isinstance(seq, str), f"dna[{index}].sequence must be a string")
    seq = seq.upper()
    _require(len(seq) > 0, f"dna[{index}].sequence must not be empty")
    _require(set(seq) <= _BASES, f"dna[{index}].sequence must contain only A/C/G/T")
    topology = item.get("topology")
    _require(topology in ("linear", "circular"),
             f"dna[{index}].topology must be 'linear' or 'circular'")
    return Dna(dna_id, seq, topology)


def _parse_enzyme(item: object, index: int) -> Enzyme:
    _require(isinstance(item, dict), f"enzymes[{index}] must be an object")
    enz_id = _require_str_id(item.get("id"), f"enzymes[{index}].id")
    rec = item.get("recognition")
    _require(isinstance(rec, str), f"enzymes[{index}].recognition must be a string")
    rec = rec.upper()
    _require(len(rec) > 0 and len(rec) % 2 == 0,
             f"enzymes[{index}].recognition length must be a positive even number")
    _require(set(rec) <= _BASES,
             f"enzymes[{index}].recognition must contain only A/C/G/T")
    _require(rec == _reverse_complement(rec),
             f"enzymes[{index}].recognition must equal its own reverse complement")
    return Enzyme(enz_id, rec)


def _parse_digest(item: object, index: int, dna_ids: set[str],
                  enzyme_ids: set[str]) -> Digest:
    _require(isinstance(item, dict), f"digests[{index}] must be an object")
    dig_id = _require_str_id(item.get("id"), f"digests[{index}].id")
    dna_id = item.get("dna")
    _require(isinstance(dna_id, str) and dna_id in dna_ids,
             f"digests[{index}].dna must reference a known dna id")
    enzymes = item.get("enzymes")
    _require(isinstance(enzymes, list), f"digests[{index}].enzymes must be a list")
    seen: set[str] = set()
    for e in enzymes:
        _require(isinstance(e, str) and e in enzyme_ids,
                 f"digests[{index}].enzymes entries must reference known enzyme ids")
        _require(e not in seen, f"digests[{index}].enzymes must not contain duplicates")
        seen.add(e)
    return Digest(dig_id, dna_id, tuple(enzymes))


def _parse_list(data: dict, key: str) -> list:
    value = data.get(key)
    _require(isinstance(value, list), f"'{key}' must be a list")
    return value


def _check_unique(ids: list[str], what: str) -> None:
    _require(len(ids) == len(set(ids)), f"{what} ids must be unique")


def load_task(path: str) -> tuple[list[Dna], list[Enzyme], list[Digest]]:
    """Load and fully validate a task JSON file. Raises InputError on any problem."""
    with open(path, "r", encoding="utf-8") as fh:
        try:
            data = json.load(fh)
        except json.JSONDecodeError as exc:
            raise InputError(f"invalid JSON: {exc}") from exc
    _require(isinstance(data, dict), "top-level JSON must be an object")

    dnas = [_parse_dna(x, i) for i, x in enumerate(_parse_list(data, "dna"))]
    enzymes = [_parse_enzyme(x, i) for i, x in enumerate(_parse_list(data, "enzymes"))]
    _check_unique([d.id for d in dnas], "dna")
    _check_unique([e.id for e in enzymes], "enzyme")
    digests = [_parse_digest(x, i, {d.id for d in dnas}, {e.id for e in enzymes})
               for i, x in enumerate(_parse_list(data, "digests"))]
    _check_unique([d.id for d in digests], "digest")
    return dnas, enzymes, digests
