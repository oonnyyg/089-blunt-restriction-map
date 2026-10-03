"""Core digest logic: site matching, cutting, fragment generation."""

from __future__ import annotations

from dataclasses import dataclass, field

from .validate import Dna, Enzyme


@dataclass(frozen=True)
class CutSite:
    position: int
    enzymes: tuple[str, ...]  # enzyme ids, sorted
    match_starts: tuple[int, ...]  # recognition start per enzyme, same order


@dataclass(frozen=True)
class Fragment:
    start: int
    length: int
    sequence: str
    wraps_origin: bool


@dataclass
class DigestResult:
    digest_id: str
    dna_id: str
    topology: str
    dna_length: int
    cut_sites: list[CutSite] = field(default_factory=list)
    fragments: list[Fragment] = field(default_factory=list)


def _find_matches(seq: str, recognition: str, circular: bool) -> list[int]:
    """Return all recognition start positions, including overlapping matches."""
    n = len(seq)
    k = len(recognition)
    if k > n:
        return []
    haystack = seq + seq[: k - 1] if circular else seq
    last_start = n - 1 if circular else n - k
    starts: list[int] = []
    pos = haystack.find(recognition, 0, last_start + k)
    while pos != -1 and pos <= last_start:
        starts.append(pos)
        pos = haystack.find(recognition, pos + 1, last_start + k)
    return starts


def compute_cut_sites(dna: Dna, enzymes: list[Enzyme]) -> list[CutSite]:
    """Cut position = match start + half recognition length (mod length if circular)."""
    n = len(dna.sequence)
    circular = dna.topology == "circular"
    by_position: dict[int, list[tuple[str, int]]] = {}
    for enzyme in enzymes:
        half = len(enzyme.recognition) // 2
        for start in _find_matches(dna.sequence, enzyme.recognition, circular):
            cut = (start + half) % n if circular else start + half
            by_position.setdefault(cut, []).append((enzyme.id, start))
    sites = []
    for position in sorted(by_position):
        hits = sorted(by_position[position])
        sites.append(CutSite(position,
                             tuple(h[0] for h in hits),
                             tuple(h[1] for h in hits)))
    return sites


def _circular_subseq(seq: str, start: int, length: int) -> str:
    n = len(seq)
    start %= n
    if start + length <= n:
        return seq[start:start + length]
    return (seq[start:] + seq)[:length]


def compute_fragments(dna: Dna, cut_positions: list[int]) -> list[Fragment]:
    seq = dna.sequence
    n = len(seq)
    if dna.topology == "linear":
        boundaries = [0, *cut_positions, n]
        return [Fragment(boundaries[i], boundaries[i + 1] - boundaries[i],
                         seq[boundaries[i]:boundaries[i + 1]], False)
                for i in range(len(boundaries) - 1)]
    if not cut_positions:
        return [Fragment(0, n, seq, False)]
    fragments = []
    for i, start in enumerate(cut_positions):
        end = cut_positions[i + 1] if i + 1 < len(cut_positions) else cut_positions[0] + n
        length = end - start
        fragments.append(Fragment(start, length,
                                  _circular_subseq(seq, start, length),
                                  start + length > n))
    return fragments


def run_digest(digest_id: str, dna: Dna, enzymes: list[Enzyme]) -> DigestResult:
    sites = compute_cut_sites(dna, enzymes)
    fragments = compute_fragments(dna, [s.position for s in sites])
    assert sum(f.length for f in fragments) == len(dna.sequence)
    return DigestResult(digest_id, dna.id, dna.topology, len(dna.sequence),
                        sites, fragments)
