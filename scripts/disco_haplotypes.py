#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
disco_haplotypes.py

Locus-level post-processing of DiscoSnp++ / DiscoSnpRad bubbles.  The de Bruijn
graph discovery (kissnp2) is untouched: every bubble stays a pair of paths.
This script groups the bubbles of a locus on a common coordinate system and
works on SITES instead of bubbles, which makes it possible to

    - merge the pairwise bubbles of a tri/tetra-allelic site into one record,
    - cope with windows holding three or more close SNPs,
    - phase the sites with the reads and report more than two haplotypes.

Three sub-commands:

  augment   BEFORE kissreads2.  For every SNP of a locus with several sites,
            writes single-SNP bubbles in the other sequence contexts seen at
            the neighbouring sites ("synthetic" bubbles).  kissreads2 only
            anchors a read with an exact seed of k-5 nt, so without this a read
            carrying other alleles at the close SNPs is simply never counted.
  call      AFTER kissreads2 (run with -phasing).  Builds the loci, merges the
            per-bubble read counts and the phased facts into per-site allele
            depths, calls multi-allelic genotypes, phases them and writes
            <out>.vcf, <out>.tsv (haplotypes, loci with >= 2 sites only),
            <out>_loci.tsv and <out>_loci.fa
  strip     Removes the synthetic bubbles from a fasta file, so that the usual
            DiscoSnp VCF is unchanged.

Scaling: bubbles are kept in flat numpy arrays (about 130 bytes + 2 x S bytes
per bubble for S read sets).  Python objects are created only for the bubbles
that belong to a locus with several bubbles.  10 million bubbles need a few GB.

INDEL bubbles are ignored (kissreads2 does not phase them); they remain in the
usual DiscoSnp VCF.
"""

import argparse
import array
import itertools
import math
import os
import re
import sys
from collections import defaultdict, Counter

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: disco_haplotypes.py needs numpy (python3 -m pip install numpy)")

PAD = 255                                   # padding of the 2D path array
CODE = np.full(256, 4, dtype=np.uint8)      # A C G T -> 0 1 2 3, anything else -> 4
for _nt, _code in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    CODE[_nt] = _code
COMP = np.full(256, PAD, dtype=np.uint8)
COMP[:5] = (3, 2, 1, 0, 4)
DECODE = bytes.maketrans(bytes(range(5)), b"ACGTN")
LOWER = b"abcdefghijklmnopqrstuvwxyz\r\n"
HEADER_RE = re.compile(rb">(SNP|INDEL)_(higher|lower)_path_(\d+)")
RANK_RE = re.compile(rb"rank_([0-9.eE+-]+)")
COUNT_RE = re.compile(rb"\|C\d+_(\d+)")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
ERROR_RATE = 0.01
MAX_PATH_LENGTH = 1000
CHUNK_ROWS = 50000                           # rows (paths) handled at once by the numpy steps


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def complement(nucleotide):
    return nucleotide.translate(COMPLEMENT)


def log(message):
    sys.stderr.write(message + "\n")
    sys.stderr.flush()


###############################################################################
# Bubble store: every SNP bubble of a fasta file in flat arrays
###############################################################################

class BubbleStore:
    """Bubble i: higher path = row 2i, lower path = row 2i+1 of P (codes, PAD padded)."""

    def __init__(self, ids, P, lens, counts, ranks):
        self.ids = ids                  # int64 [n]   bubble ids of the fasta headers
        self.P = P                      # uint8 [2n, Lmax]
        self.lens = lens                # int32 [2n]
        self.counts = counts            # uint16 [2n, S] read counts of kissreads2, or None
        self.ranks = ranks              # float32 [n]
        self.n = len(ids)
        self.nsamples = 0 if counts is None else counts.shape[1]
        self.max_len = P.shape[1] if self.n else 0
        max_id = int(ids.max()) if self.n else 0
        self.by_id = np.full(max_id + 1, -1, dtype=np.int64)
        self.by_id[ids] = np.arange(self.n)

    def index_of(self, bubble_id):
        """Index of a bubble id, -1 when absent."""
        if bubble_id < 0 or bubble_id >= len(self.by_id):
            return -1
        return int(self.by_id[bubble_id])

    def path_str(self, row):
        return bytes(self.P[row, :self.lens[row]]).translate(DECODE).decode()

    def snp_positions(self, indices):
        """For a chunk of bubble indices: list of arrays of SNP positions (higher != lower)."""
        H = self.P[2 * indices]
        L = self.P[2 * indices + 1]
        diff = (H != L) & (H != PAD) & (L != PAD)
        rows, cols = np.nonzero(diff)
        splits = np.searchsorted(rows, np.arange(1, len(indices)))
        return np.split(cols, splits)


def parse_store(fasta_file, with_counts=False, only_ids=None):
    """Read a DiscoSnp fasta file into a BubbleStore (SNP bubbles only).

    only_ids: sorted int64 array, keep these bubble ids only.
    """
    ids, ranks, lens = array.array("q"), array.array("f"), array.array("i")
    counts, buf = array.array("H"), bytearray()
    nsamples = None
    pending = None                  # (id, path, header) of a higher path waiting for its lower path
    header = None
    n_indel = 0
    with open(fasta_file, "rb") as handle:
        for line in handle:
            if line[:1] == b">":
                header = line
                continue
            if header is None:
                continue
            match = HEADER_RE.match(header)
            if match is None:
                sys.exit(f"ERROR: unexpected header in {fasta_file}: {header.decode(errors='replace').strip()}")
            this_header, header = header, None
            kind, level, bubble_id = match.group(1), match.group(2), int(match.group(3))
            if kind != b"SNP":
                n_indel += level == b"higher"
                continue
            if only_ids is not None:
                position = np.searchsorted(only_ids, bubble_id)
                if position >= len(only_ids) or only_ids[position] != bubble_id:
                    continue
            path = line.translate(None, LOWER)       # upper-case part only, fast C-level deletion
            if level == b"higher":
                pending = (bubble_id, path, this_header)
                continue
            if pending is None or pending[0] != bubble_id:
                sys.exit(f"ERROR: lower path {bubble_id} is not preceded by its higher path in {fasta_file}")
            if max(len(path), len(pending[1])) > MAX_PATH_LENGTH:
                sys.exit(f"ERROR: bubble {bubble_id} has an upper-case path longer than {MAX_PATH_LENGTH} nt")
            ids.append(bubble_id)
            for sequence in (pending[1], path):
                buf += sequence
                lens.append(len(sequence))
            rank = RANK_RE.search(pending[2])
            ranks.append(float(rank.group(1)) if rank else float("nan"))
            if with_counts:
                for text in (pending[2], this_header):
                    values = list(map(int, COUNT_RE.findall(text)))
                    if values and max(values) > 65535:
                        values = [min(v, 65535) for v in values]
                    if nsamples is None:
                        nsamples = len(values)
                    elif len(values) != nsamples:
                        sys.exit(f"ERROR: bubble {bubble_id}: {len(values)} read counts, {nsamples} expected")
                    counts.extend(values)
            pending = None
    if pending is not None:
        sys.exit(f"ERROR: higher path {pending[0]} has no lower path in {fasta_file}")

    n = len(ids)
    ids = np.array(ids, dtype=np.int64)
    lens = np.array(lens, dtype=np.int32)
    ranks = np.array(ranks, dtype=np.float32)
    max_len = int(lens.max()) if n else 0
    P = np.full((2 * n, max_len), PAD, dtype=np.uint8)
    codes = CODE[np.frombuffer(buf, dtype=np.uint8)]
    starts = np.zeros(2 * n + 1, dtype=np.int64)
    np.cumsum(lens, out=starts[1:])
    if n and (lens == max_len).all():
        P[:] = codes.reshape(2 * n, max_len)              # the usual case: closed bubbles, one length
    else:
        for r0 in range(0, 2 * n, 20000):
            r1 = min(2 * n, r0 + 20000)
            chunk_lens = lens[r0:r1]
            total = int(starts[r1] - starts[r0])
            rows = np.repeat(np.arange(r0, r1, dtype=np.int32), chunk_lens)
            cols = (np.arange(total, dtype=np.int32)
                    - np.repeat((starts[r0:r1] - starts[r0]).astype(np.int32), chunk_lens))
            P[rows, cols] = codes[starts[r0]:starts[r1]]
    del buf, codes
    if with_counts:
        if n and not nsamples:
            sys.exit(f"ERROR: no read counts (C1_, C2_...) in the headers of {fasta_file}: is this a kissreads2 output?")
        counts = np.frombuffer(counts, dtype=np.uint16).reshape(2 * n, nsamples or 0)
    else:
        counts = None
    store = BubbleStore(ids, P, lens, counts, ranks)
    log(f"[{os.path.basename(fasta_file)}] {n} SNP bubbles read ({n_indel} INDEL bubbles ignored),"
        f" longest path {max_len} nt" + (f", {nsamples} read sets" if with_counts else ""))
    return store


def read_synthetic_map(map_file):
    """{synthetic_id: (parent_id, start of the synthetic path in the parent path)}"""
    parents = {}
    if map_file and os.path.exists(map_file):
        with open(map_file) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                synthetic_id, parent_id, start = line.split()[:3]
                parents[int(synthetic_id)] = (int(parent_id), int(start))
    return parents


###############################################################################
# Placement edges from the sequences (numpy)
#
# An edge (b1, b2, shift, relative) reads: "position 0 of the path of b2 lies at
# position <shift> of the path of b1, and b2 runs in the same (+1) or in the
# opposite (-1) direction".  Position 0 is shared by the two paths of a bubble.
###############################################################################

def window_codes(rows, k):
    """2-bit codes of every k-window of every row: (codes uint32 [m, L-k+1], valid mask)."""
    m, L = rows.shape
    nw = L - k + 1
    if nw <= 0:
        return np.zeros((m, 0), dtype=np.uint32), np.zeros((m, 0), dtype=bool)
    codes = np.zeros((m, nw), dtype=np.uint32)
    valid = np.ones((m, nw), dtype=bool)
    for j in range(k):
        column = rows[:, j:j + nw]
        valid &= column < 4
        codes = (codes << np.uint32(2)) | (column & 3).astype(np.uint32)
    return codes, valid


def revcomp_rows(rows, lens):
    """Reverse complement of every row (its first `lens` codes), same padding."""
    m, L = rows.shape
    index = lens[:, None].astype(np.int64) - 1 - np.arange(L)[None, :]
    ok = index >= 0
    out = np.full_like(rows, PAD)
    r, c = np.nonzero(ok)
    out[r, c] = COMP[rows[r, index[r, c]]]
    return out


def pack_keys(qb, tb, o, rel):
    """(query bubble, target bubble, oriented offset, orientation) -> one uint64."""
    flipped = np.broadcast_to(np.asarray(rel) < 0, qb.shape).astype(np.uint64)
    return ((qb.astype(np.uint64) << np.uint64(38)) | (tb.astype(np.uint64) << np.uint64(12))
            | ((o + 1024).astype(np.uint64) << np.uint64(1)) | flipped)


def unpack_keys(keys):
    qb = (keys >> np.uint64(38)).astype(np.int64)
    tb = ((keys >> np.uint64(12)) & np.uint64((1 << 26) - 1)).astype(np.int64)
    o = ((keys >> np.uint64(1)) & np.uint64(2047)).astype(np.int64) - 1024
    rel = np.where(keys & np.uint64(1), -1, 1).astype(np.int64)
    return qb, tb, o, rel


def sequence_edges(store, seed_size, max_mismatches, min_overlap, max_candidates=64):
    """Edges between bubbles whose paths overlap, from the sequences alone.

    Links the pairwise bubbles of a multi-allelic site (they start at the same
    place, and kissreads2 never reports two bubbles starting at the same read
    position in a fact) and is the only source of links before kissreads2 runs.
    Seeds are taken every seed_size positions (and at the end) of every path
    and indexed in sorted arrays; every path is then scanned in both directions
    and a candidate placement is accepted when the best of the 2 x 2 path pairs
    agrees over the whole overlap, up to max_mismatches (other alleles of the
    locus) and 10 % of the overlap.
    Returns (b1, b2, shift, rel) int64 arrays with b1 < b2 (bubble indices).
    """
    empty = tuple(np.zeros(0, dtype=np.int64) for _ in range(4))
    n, P, lens, k = store.n, store.P, store.lens, seed_size
    if n < 2 or store.max_len < k:
        return empty
    if store.n >= (1 << 25):
        sys.exit("ERROR: more than 33 million bubbles are not supported by the edge search")

    # ---- index: seeds at offsets 0, k, 2k... and at the end of every path
    seeds, rows, offs = [], [], []
    for r0 in range(0, 2 * n, CHUNK_ROWS):
        r1 = min(2 * n, r0 + CHUNK_ROWS)
        codes, valid = window_codes(P[r0:r1], k)
        last = lens[r0:r1].astype(np.int64) - k
        for offset in range(0, codes.shape[1], k):
            selected = (last >= offset) & valid[:, offset]
            seeds.append(codes[selected, offset])
            rows.append(np.nonzero(selected)[0] + r0)
            offs.append(np.full(int(selected.sum()), offset, dtype=np.int16))
        selected = (last >= 0) & (last % k != 0)
        r = np.nonzero(selected)[0]
        l = last[r]
        ok = valid[r, l]
        seeds.append(codes[r[ok], l[ok]])
        rows.append(r[ok] + r0)
        offs.append(l[ok].astype(np.int16))
    seeds, rows, offs = (np.concatenate(a) for a in (seeds, rows, offs))
    order = np.argsort(seeds, kind="stable")
    seeds, rows, offs = seeds[order], rows[order].astype(np.int32), offs[order]
    del order
    _, counts = np.unique(seeds, return_counts=True)
    keep = np.repeat(counts <= max_candidates, counts)
    n_repetitive = int((counts > max_candidates).sum())
    seeds, rows, offs = seeds[keep], rows[keep], offs[keep]
    del keep, counts

    # ---- query: every window of every path, both orientations
    all_keys = []
    n_hits = 0
    for b0 in range(0, n, CHUNK_ROWS // 2):
        r0, r1 = 2 * b0, 2 * min(n, b0 + CHUNK_ROWS // 2)
        chunk, chunk_lens = P[r0:r1], lens[r0:r1]
        for rel, oriented in ((1, chunk), (-1, revcomp_rows(chunk, chunk_lens))):
            codes, valid = window_codes(oriented, k)
            qrow, qpos = np.nonzero(valid)
            qcodes = codes[qrow, qpos]
            lo = np.searchsorted(seeds, qcodes, "left")
            hi = np.searchsorted(seeds, qcodes, "right")
            hits = hi - lo
            has = hits > 0
            qrow, qpos, lo, hits = qrow[has], qpos[has], lo[has], hits[has]
            if not len(hits):
                continue
            cumulative = np.cumsum(hits)
            n_hits += int(cumulative[-1])
            bounds = [0] + [int(i) + 1 for i in np.nonzero(np.diff(cumulative // 20_000_000))[0]] + [len(hits)]
            for i0, i1 in zip(bounds, bounds[1:]):
                if i1 <= i0:
                    continue
                c = hits[i0:i1]
                total = int(c.sum())
                index = np.repeat(lo[i0:i1], c) + (np.arange(total) - np.repeat(np.cumsum(c) - c, c))
                qb = (np.repeat(qrow[i0:i1], c) + r0) // 2
                tb = rows[index].astype(np.int64) // 2
                keep = qb != tb
                o = (np.repeat(qpos[i0:i1], c) - offs[index]).astype(np.int64)
                all_keys.append(np.unique(pack_keys(qb[keep], tb[keep], o[keep], rel)))
    del seeds, rows, offs
    if not all_keys:
        return empty
    keys = np.unique(np.concatenate(all_keys))
    del all_keys
    qb, tb, o, rel = unpack_keys(keys)
    del keys

    # ---- verification: best of the 2 x 2 path pairs, grouped by (offset, orientation)
    best_mm = np.full(len(qb), 10 ** 6, dtype=np.int32)
    best_ov = np.zeros(len(qb), dtype=np.int32)
    groups = np.unique(np.stack([rel, o], axis=1), axis=0)
    for rel_v, o_v in groups.tolist():
        selected = np.nonzero((rel == rel_v) & (o == o_v))[0]
        a, b = max(0, o_v), max(0, -o_v)
        width = store.max_len - max(a, b)
        if width <= 0:
            continue
        group_mm = np.full(len(selected), 10 ** 6, dtype=np.int32)
        group_ov = np.zeros(len(selected), dtype=np.int32)
        for qi in (0, 1):
            Q = P[2 * qb[selected] + qi]
            if rel_v == -1:
                Q = revcomp_rows(Q, lens[2 * qb[selected] + qi])
            QS = Q[:, a:a + width]
            for ti in (0, 1):
                TS = P[2 * tb[selected] + ti][:, b:b + width]
                valid = (QS != PAD) & (TS != PAD)
                overlap = valid.sum(axis=1).astype(np.int32)
                mismatches = ((QS != TS) & valid).sum(axis=1).astype(np.int32)
                better = (mismatches < group_mm) | ((mismatches == group_mm) & (overlap > group_ov))
                group_mm = np.where(better, mismatches, group_mm)
                group_ov = np.where(better, overlap, group_ov)
        best_mm[selected], best_ov[selected] = group_mm, group_ov
    accepted = (best_ov >= min_overlap) & (best_mm <= max_mismatches) & (best_mm * 10 <= best_ov)
    qb, tb, o, rel, best_mm, best_ov = (a[accepted] for a in (qb, tb, o, rel, best_mm, best_ov))
    if not len(qb):
        return empty
    shift = np.where(rel == 1, o, lens[2 * qb].astype(np.int64) - 1 - o)

    # ---- one placement per bubble pair, normalised to b1 < b2
    order = np.lexsort((-best_ov, best_mm, tb, qb))
    qb, tb, shift, rel = qb[order], tb[order], shift[order], rel[order]
    first = np.ones(len(qb), dtype=bool)
    first[1:] = (qb[1:] != qb[:-1]) | (tb[1:] != tb[:-1])
    qb, tb, shift, rel = qb[first], tb[first], shift[first], rel[first]
    swap = qb > tb
    b1, b2 = np.where(swap, tb, qb), np.where(swap, qb, tb)
    shift = np.where(swap, -rel * shift, shift)
    keys = np.unique(pack_keys(b1, b2, shift, rel))
    b1, b2, shift, rel = unpack_keys(keys)
    log(f"[edges] {len(b1)} sequence overlaps between bubbles ({n_hits} seed hits,"
        f" {n_repetitive} repetitive seeds ignored)")
    return b1, b2, shift, rel


###############################################################################
# Facts (kissreads2 -phasing)
###############################################################################

FACT_TOKEN = re.compile(r"(-?)(\d+)([hl])_(-?\d+)$")


def parse_fact_part(text):
    """'-1h_0;-3h_-51;' -> [(-1, 1, 0, 0), (-1, 3, 0, -51)] = (sign, bubble id, path, gap)"""
    tokens = []
    for token in text.strip().split(";"):
        if not token:
            continue
        match = FACT_TOKEN.match(token)
        if match is None:
            return None
        sign, bubble_id, level, gap = match.groups()
        tokens.append((-1 if sign else 1, int(bubble_id), 0 if level == "h" else 1, int(gap)))
    return tokens


def sample_of_fact_file(fact_file):
    match = re.search(r"read_set_id_(\d+)", os.path.basename(fact_file))
    if match is None:
        sys.exit(f"ERROR: cannot find the read set index in the name of {fact_file}")
    return int(match.group(1)) - 1            # 0-based sample index (C1_ -> 0)


def read_set_name(fact_file):
    """The read set name written by kissreads2 on the first line of a fact file."""
    with open(fact_file) as handle:
        first = handle.readline()
    return first[1:].strip() if first.startswith("#") else "."


def iter_facts(fact_file):
    """Yields (parts, support); parts = tuple of parts (one per mate), a part = tuple of tokens."""
    with open(fact_file) as handle:
        for line in handle:
            if line.startswith("#") or "=>" not in line:
                continue
            text, support = line.rsplit("=>", 1)
            parts = tuple(tuple(p) for p in map(parse_fact_part, text.split()) if p)
            if parts:
                yield parts, int(support)


SITE_TOKEN = re.compile(r"(\d+)([hl]):([01]+)$")


def site_fact_observations(text, by_id):
    """'12h:101;15l:1;' (kissreads2 -phasing_sites) -> {locus name: {coordinate: nucleotide}}.

    Every SNP path mapped by the read(s) is listed with the SNPs of the path the
    read covers (1) or not (0), so every observation is real: no multi-SNP path
    is trusted beyond what the read saw, and no bubble is lost because it starts
    at the same place as, or lies inside, another one.
    Returns (observations, number of contradictions).
    """
    per_locus = {}
    n_conflicts = 0
    for token in text.strip().split(";"):
        if not token:
            continue
        match = SITE_TOKEN.match(token)
        if match is None:
            continue
        bubble = by_id.get(int(match.group(1)))
        if bubble is None or len(match.group(3)) != len(bubble.snps):
            continue
        path_index = 0 if match.group(2) == "h" else 1
        observed = per_locus.setdefault(bubble.locus.name, {})
        for (position, _, _), covered in zip(bubble.snps, match.group(3)):
            if covered != "1":
                continue
            coordinate = bubble.coordinate(position)
            nucleotide = bubble.nucleotide(path_index, coordinate)
            if observed.get(coordinate, nucleotide) != nucleotide:
                nucleotide = None          # two paths disagree: drop this site
                n_conflicts += 1
            observed[coordinate] = nucleotide
    return per_locus, n_conflicts


def fact_edges(fact_files, store):
    """Placement edges given by the reads.

    In the frame of the read, bubble i covers [start_i, start_i + len_i) with
    start_(i+1) = start_i + len_i + gap_(i+1)   (kissreads2 'shift'; negative
    when the paths overlap).  With sign -1 the path runs leftwards, so its
    position 0 is at start_i + len_i - 1.
    Returns (b1, b2, shift, rel) arrays, b1 < b2.
    """
    keys = set()
    seen = set()                  # fact texts already processed: the same facts recur in every read set
    lens = store.lens
    n_unknown = 0
    for fact_file in fact_files:
        with open(fact_file) as handle:
            lines = (line.rsplit("=>", 1)[0] for line in handle if not line.startswith("#") and "=>" in line)
            texts = [text for text in lines if text not in seen and not seen.add(text)]
        for text in texts:
            parts = tuple(tuple(p) for p in map(parse_fact_part, text.split()) if p)
            for part in parts:
                previous = None
                start = 0
                for sign, bubble_id, path_index, gap in part:
                    index = store.index_of(bubble_id)
                    if index < 0:
                        n_unknown += 1
                        previous = None
                        continue
                    length = int(lens[2 * index + path_index])
                    if previous is not None:
                        start = previous[0] + previous[1] + gap
                    zero = start if sign == 1 else start + length - 1
                    if previous is not None:
                        p_sign, p_index, p_zero = previous[2], previous[3], previous[4]
                        shift, rel = p_sign * (zero - p_zero), p_sign * sign
                        # the edge is stored with its lower index first; the chain itself must keep
                        # following the read (do not overwrite index)
                        b1, b2 = (p_index, index) if p_index < index else (index, p_index)
                        if p_index > index:
                            shift = -rel * shift
                        keys.add((b1 << 38) | (b2 << 12) | ((shift + 1024) << 1) | (rel < 0))
                    previous = (start, length, sign, index, zero)
    if n_unknown:
        log(f"[facts] {n_unknown} fact tokens name bubbles absent from the fasta files (ignored)")
    if not keys:
        return tuple(np.zeros(0, dtype=np.int64) for _ in range(4))
    return unpack_keys(np.array(sorted(keys), dtype=np.uint64))


###############################################################################
# Loci
###############################################################################

class Bubble:
    """A bubble materialised from the store (only for bubbles of multi-bubble loci)."""

    __slots__ = ("id", "index", "store", "paths", "snps", "headers", "sequences",
                 "parent", "anchor", "orient", "locus")

    def __init__(self, store, index):
        self.id = int(store.ids[index])
        self.index = index
        self.store = store
        self.paths = [store.path_str(2 * index), store.path_str(2 * index + 1)]
        high, low = self.paths
        self.snps = [(p, high[p], low[p]) for p in range(min(len(high), len(low))) if high[p] != low[p]]
        self.headers = None             # filled by augment (second pass over the fasta)
        self.sequences = None
        self.parent = None              # id of the real bubble if this one is synthetic
        self.anchor = None              # locus coordinate of path position 0
        self.orient = 1                 # +1 / -1 : direction of the path in the locus
        self.locus = None

    @property
    def rank(self):
        value = float(self.store.ranks[self.index])
        return None if math.isnan(value) else value

    def count(self, path_index, sample):
        return int(self.store.counts[2 * self.index + path_index, sample])

    def count_vector(self, path_index):
        return self.store.counts[2 * self.index + path_index]

    def coordinate(self, position):
        return self.anchor + self.orient * position

    def position(self, coordinate):
        return self.orient * (coordinate - self.anchor)

    def nucleotide(self, path_index, coordinate):
        """Nucleotide of a path at a locus coordinate, in locus orientation (None if outside)."""
        position = self.position(coordinate)
        path = self.paths[path_index]
        if position < 0 or position >= len(path):
            return None
        return path[position] if self.orient == 1 else complement(path[position])


class Locus:
    __slots__ = ("name", "bubbles", "conflicts", "sites", "length")

    def __init__(self):
        self.name = None
        self.bubbles = []
        self.conflicts = 0
        self.sites = {}      # {coordinate: {nucleotide: [(bubble, path_index), ...]}}
        self.length = 0


def connected_components(b1, b2):
    """{root: [bubble indices]} for the bubbles appearing in at least one edge."""
    parent = {}

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for a, b in zip(b1.tolist(), b2.tolist()):
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    components = defaultdict(list)
    for node in parent:
        components[find(node)].append(node)
    return components


def build_loci(bubbles, edges):
    """Place every materialised bubble (anchor, orient) in its locus; returns the loci.

    bubbles: {index: Bubble};  edges: (b1, b2, shift, rel) arrays, restricted to
    those bubbles.  Edges are walked from the smallest index of each component.
    """
    graph = defaultdict(list)
    for a, b, shift, rel in zip(*(e.tolist() for e in edges)):
        if a in bubbles and b in bubbles:
            graph[a].append((b, shift, rel))
            graph[b].append((a, -rel * shift, rel))

    loci = []
    for seed_index in sorted(bubbles):
        seed = bubbles[seed_index]
        if seed.locus is not None:
            continue
        locus = Locus()
        seed.locus, seed.anchor, seed.orient = locus, 0, 1
        stack = [seed]
        while stack:
            current = stack.pop()
            locus.bubbles.append(current)
            for neighbour_index, shift, relative in graph[current.index]:
                neighbour = bubbles[neighbour_index]
                anchor = current.anchor + current.orient * shift
                orient = current.orient * relative
                if neighbour.locus is None:
                    neighbour.locus, neighbour.anchor, neighbour.orient = locus, anchor, orient
                    stack.append(neighbour)
                elif (neighbour.anchor, neighbour.orient) != (anchor, orient):
                    locus.conflicts += 1
        loci.append(locus)

    for locus in loci:
        locus.conflicts //= 2   # every edge is walked from both sides
        locus.bubbles.sort(key=lambda b: b.index)
        lowest = min(min(b.anchor, b.coordinate(max(map(len, b.paths)) - 1)) for b in locus.bubbles)
        highest = 0
        for bubble in locus.bubbles:
            bubble.anchor -= lowest
            highest = max(highest, bubble.anchor, bubble.coordinate(max(map(len, bubble.paths)) - 1))
        locus.length = highest + 1
        sites = defaultdict(lambda: defaultdict(list))
        for bubble in locus.bubbles:
            for position, high, low in bubble.snps:
                for path_index, nucleotide in ((0, high), (1, low)):
                    if bubble.orient == -1:
                        nucleotide = complement(nucleotide)
                    sites[bubble.coordinate(position)][nucleotide].append((bubble, path_index))
        locus.sites = {c: dict(sites[c]) for c in sorted(sites)}
    return loci


def locus_sequence(locus, reference_alleles):
    """Consensus of the locus: the paths laid on the coordinates, REF allele at the sites."""
    sequence = ["N"] * locus.length
    for bubble in locus.bubbles:
        for position in range(len(bubble.paths[0])):
            coordinate = bubble.coordinate(position)
            if sequence[coordinate] == "N":
                sequence[coordinate] = bubble.nucleotide(0, coordinate)
    for coordinate, nucleotide in reference_alleles.items():
        sequence[coordinate] = nucleotide
    return "".join(sequence)


def snp_counts(store):
    """Number of SNPs (positions where the two paths differ) of every bubble."""
    counts = np.zeros(store.n, dtype=np.int32)
    for i0 in range(0, store.n, CHUNK_ROWS // 2):
        chunk = np.arange(i0, min(store.n, i0 + CHUNK_ROWS // 2))
        H, L = store.P[2 * chunk], store.P[2 * chunk + 1]
        counts[chunk] = ((H != L) & (H != PAD) & (L != PAD)).sum(axis=1)
    return counts


def multi_bubble_loci(store, edges, max_locus_bubbles, lone_multi_snp=False):
    """Materialise the bubbles of the components (2..max_locus_bubbles bubbles), build their loci.

    lone_multi_snp: also make a locus of every bubble linked to no other one but
    holding several SNPs (needed by 'augment': the recombinant haplotypes of a
    lone multi-SNP bubble are only counted through its single-SNP versions).
    Returns (loci, {index: Bubble}).
    """
    components = connected_components(edges[0], edges[1])
    n_oversized = 0
    materialised = {}
    if lone_multi_snp:
        linked = np.zeros(store.n, dtype=bool)
        linked[edges[0]] = True
        linked[edges[1]] = True
        for index in np.nonzero(~linked & (snp_counts(store) > 1))[0].tolist():
            materialised[index] = Bubble(store, index)
    for members in components.values():
        if len(members) > max_locus_bubbles:
            n_oversized += len(members)
            continue
        for index in members:
            materialised[index] = Bubble(store, index)
    loci = build_loci(materialised, edges)
    if n_oversized:
        log(f"[loci] {n_oversized} bubbles belong to {sum(1 for m in components.values() if len(m) > max_locus_bubbles)}"
            f" components larger than --max_locus_bubbles ({max_locus_bubbles}): probably repeats, treated as isolated bubbles")
    return loci, materialised


###############################################################################
# augment
###############################################################################

def canonical_pair(path_1, path_2):
    return frozenset((min(path_1, revcomp(path_1)), min(path_2, revcomp(path_2))))


def sub_bubble(bubble, snp_index, locus, max_contexts):
    """Single-SNP versions of one SNP of a bubble, in every context seen in the locus.

    The paths are cropped to the flanks kissnp2 gives to an isolated SNP (k-1 on
    each side when the bubble is closed), so that every read mapped by kissreads2
    (overlap of at least k) covers the SNP: the read counts become site specific.
    This is not the case for a multi-SNP bubble, where a read covering the last
    SNP only is still counted for the whole path.
    Yields (start_in_parent_path, context_label, [sequence_higher, sequence_lower], capped).
    """
    position = bubble.snps[snp_index][0]
    start = position - bubble.snps[0][0]
    site = bubble.coordinate(position)
    templates = []
    for path_index in (0, 1):
        path = bubble.paths[path_index]
        stop = position + (len(path) - 1 - bubble.snps[-1][0])
        full = bubble.sequences[path_index]
        upper = [i for i, c in enumerate(full) if c.isupper()]
        left = full[:upper[0]] if start == 0 else ""
        right = full[upper[-1] + 1:] if stop == len(path) - 1 else ""
        templates.append((left, list(path[start:stop + 1]), right))
    shortest = min(len(t[1]) for t in templates)
    inside = [c for c in locus.sites
              if c != site and 0 <= bubble.position(c) - start < shortest]
    inside.sort(key=lambda c: abs(c - site))   # nearest first: they are the ones breaking the seeds
    used, n_contexts, capped = [], 1, False
    for coordinate in inside:
        if n_contexts * len(locus.sites[coordinate]) > max_contexts:
            capped = True
            break
        n_contexts *= len(locus.sites[coordinate])
        used.append(coordinate)
    for context in itertools.product(*(sorted(locus.sites[c]) for c in used)):
        sequences = []
        for path_index, (left, template, right) in enumerate(templates):
            path = list(template)
            for coordinate, nucleotide in zip(used, context):
                if bubble.orient == -1:
                    nucleotide = complement(nucleotide)
                path[bubble.position(coordinate) - start] = nucleotide
            path[position - start] = bubble.snps[snp_index][1 + path_index]
            sequences.append(left + "".join(path) + right)
        label = ",".join(f"{c + 1}:{n}" for c, n in zip(used, context))
        yield start, label, sequences, capped


def augment(args):
    store = parse_store(args.input)
    edges = sequence_edges(store, args.seed_size, args.max_mismatches, args.min_overlap)
    loci, materialised = multi_bubble_loci(store, edges, args.max_locus_bubbles, lone_multi_snp=True)
    loci = [locus for locus in loci if len(locus.sites) > 1]
    for number, locus in enumerate(loci, 1):
        locus.name = f"locus_{number}"
    wanted = {bubble.id: bubble for locus in loci for bubble in locus.bubbles}
    del store, edges
    next_id = 0

    # second pass: copy the input, keep the full records of the bubbles of multi-site loci
    with open(args.output, "w") as fasta, open(args.input, "rb") as handle:
        header = None
        for line in handle:
            fasta.write(line.decode())
            if line[:1] == b">":
                header = line
                continue
            if header is None:
                continue
            match = HEADER_RE.match(header)
            if match:
                # kissnp2 numbers SNPs and INDELs with one counter: synthetic ids must follow both
                next_id = max(next_id, int(match.group(3)) + 1)
            bubble = wanted.get(int(match.group(3))) if match else None
            if bubble is not None:
                path_index = 0 if match.group(2) == b"higher" else 1
                if bubble.headers is None:
                    bubble.headers, bubble.sequences = [None, None], [None, None]
                bubble.headers[path_index] = header.decode().rstrip("\n")
                bubble.sequences[path_index] = line.decode().rstrip("\n")
            header = None

        n_synthetic = n_capped = 0
        with open(args.map, "w") as mapping:
            mapping.write("#synthetic_id\tparent_id\tstart_in_parent\tlocus\tsite\tcontext\n")
            for locus in loci:
                existing = {canonical_pair(*b.paths) for b in locus.bubbles}
                for bubble in locus.bubbles:
                    for snp_index, (position, _, _) in enumerate(bubble.snps):
                        capped = False
                        for start, label, sequences, capped in sub_bubble(bubble, snp_index, locus,
                                                                          args.max_contexts):
                            pair = canonical_pair(*("".join(c for c in s if c.isupper()) for s in sequences))
                            if pair in existing:
                                continue
                            existing.add(pair)
                            for path_index, level in ((0, "higher"), (1, "lower")):
                                header = re.sub(r"^>SNP_(higher|lower)_path_\d+",
                                                f">SNP_{level}_path_{next_id}", bubble.headers[path_index])
                                fasta.write(f"{header}\n{sequences[path_index]}\n")
                            mapping.write(f"{next_id}\t{bubble.id}\t{start}\t{locus.name}"
                                          f"\t{bubble.coordinate(position) + 1}\t{label or '.'}\n")
                            next_id += 1
                            n_synthetic += 1
                        n_capped += capped
    log(f"[augment] {len(loci)} loci with several sites ({len(wanted)} bubbles)\n"
        f"[augment] {n_synthetic} synthetic single-SNP context bubbles written"
        f" ({n_capped} SNPs limited by --max_contexts)")


###############################################################################
# call: genotyping
###############################################################################

def genotype_model(n_alleles):
    """VCF-ordered genotypes and the log10 weight matrix W [G, A] with LL = D @ W.T"""
    genotypes = [(a, b) for b in range(n_alleles) for a in range(b + 1)]
    W = np.zeros((len(genotypes), n_alleles))
    for g, (a, b) in enumerate(genotypes):
        for i in range(n_alleles):
            p_a = 1 - ERROR_RATE if i == a else ERROR_RATE / 3
            p_b = 1 - ERROR_RATE if i == b else ERROR_RATE / 3
            W[g, i] = math.log10((p_a + p_b) / 2)
    return genotypes, W


MODELS = {A: genotype_model(A) for A in (2, 3, 4)}


def call_genotypes(depths, n_alleles, min_depth):
    """Vectorised diploid calling.  depths: [n, A] -> (gt [n, 2] int8, PL [n, G] uint16, GQ [n] uint8).

    gt is -1 where the depth is below min_depth.
    """
    genotypes, W = MODELS[n_alleles]
    LL = depths.astype(np.float64) @ W.T
    best = LL.max(axis=1, keepdims=True)
    PL = np.minimum(np.rint(-10 * (LL - best)), 9999).astype(np.uint16)
    which = PL.argmin(axis=1)
    gt = np.array(genotypes, dtype=np.int8)[which]
    second = np.partition(PL, 1, axis=1)[:, 1] if PL.shape[1] > 1 else np.full(len(PL), 99)
    GQ = np.minimum(second, 99).astype(np.uint8)
    missing = depths.sum(axis=1) < min_depth
    gt[missing] = -1
    GQ[missing] = 0
    return gt, PL, GQ


def format_field(gt, phased, ps, depths, gq, pl):
    """One VCF sample field GT:PS:DP:AD:GQ:PL from python ints."""
    ad = ",".join(map(str, depths))
    if gt[0] < 0:
        return f"./.:.:{sum(depths)}:{ad}:.:."
    if phased:
        return f"{gt[0]}|{gt[1]}:{ps}:{sum(depths)}:{ad}:{gq}:{','.join(map(str, pl))}"
    return f"{min(gt)}/{max(gt)}:.:{sum(depths)}:{ad}:{gq}:{','.join(map(str, pl))}"


###############################################################################
# call: read-backed phasing
###############################################################################

class ParityUnionFind:
    """Union-find keeping, for every site, whether it is flipped with respect to its root."""

    def __init__(self, items):
        self.parent = {i: i for i in items}
        self.parity = {i: 0 for i in items}

    def find(self, item):
        path = []
        while self.parent[item] != item:
            path.append(item)
            item = self.parent[item]
        parity = 0
        for node in reversed(path):
            parity ^= self.parity[node]
            self.parent[node], self.parity[node] = item, parity
        return item, parity

    def union(self, item_1, item_2, flipped):
        root_1, parity_1 = self.find(item_1)
        root_2, parity_2 = self.find(item_2)
        if root_1 == root_2:
            return (parity_1 ^ parity_2) == flipped
        self.parent[root_2] = root_1
        self.parity[root_2] = parity_1 ^ parity_2 ^ flipped
        return True


def add_pair_votes(observed, support, votes):
    """observed {coordinate: nucleotide} of one fragment -> votes[(c1, c2, n1, n2)] += support"""
    coordinates = sorted(observed)
    for i, c1 in enumerate(coordinates):
        for c2 in coordinates[i + 1:]:
            votes[(c1, c2, observed[c1], observed[c2])] += support


def phase_from_votes(genotypes, pair_votes, min_support):
    """Phase the heterozygous sites of one sample at one locus.

    genotypes  : {coordinate: (nt_1, nt_2)}          (called sites only)
    pair_votes : {(c1, c2, nt1, nt2): support}       alleles seen together on fragments
    Every fragment covering two heterozygous sites votes for 'cis' or 'trans';
    the pairs are joined from the best supported to the least supported.
    Returns ({coordinate: (nt_hap1, nt_hap2)}, set of phased coordinates).
    """
    het = sorted(c for c, (a, b) in genotypes.items() if a != b)
    votes = Counter()
    if len(het) > 1:
        het_set = set(het)
        for (c1, c2, n1, n2), support in pair_votes.items():
            if c1 in het_set and c2 in het_set and n1 in genotypes[c1] and n2 in genotypes[c2]:
                same = (n1 == genotypes[c1][0]) == (n2 == genotypes[c2][0])
                votes[(c1, c2)] += support if same else -support
    forest = ParityUnionFind(het)
    for (c1, c2), vote in sorted(votes.items(), key=lambda kv: -abs(kv[1])):
        if abs(vote) >= min_support:
            forest.union(c1, c2, 0 if vote > 0 else 1)
    blocks = defaultdict(list)
    for coordinate in het:
        blocks[forest.find(coordinate)[0]].append(coordinate)
    main = max(blocks.values(), key=len) if blocks else []      # one phase set: the largest block
    haplotypes = {}
    for coordinate, (a, b) in genotypes.items():
        if a != b and coordinate in main and forest.find(coordinate)[1] == 1:
            a, b = b, a
        haplotypes[coordinate] = (a, b)
    phased = {c for c in genotypes if genotypes[c][0] == genotypes[c][1]} | set(main)
    return haplotypes, phased


###############################################################################
# call: loci with several bubbles
###############################################################################

class LocusData:
    """Sample-independent data of a multi-bubble locus, and its slots in the global arrays."""

    __slots__ = ("locus", "coordinates", "alleles", "site_base", "sa_base", "pl_base",
                 "site_specific", "multi_snp", "n_real", "n_synthetic")

    def __init__(self, locus):
        self.locus = locus
        self.coordinates = list(locus.sites)
        self.alleles = {}
        for coordinate, carriers in locus.sites.items():
            first = min((b for cs in carriers.values() for b, _ in cs), key=lambda b: b.index)
            reference = first.nucleotide(0, coordinate)
            self.alleles[coordinate] = [reference] + sorted(n for n in carriers if n != reference)
        self.site_specific = {c for c, carriers in locus.sites.items()
                              if any(len(b.snps) == 1 for cs in carriers.values() for b, _ in cs)}
        self.multi_snp = [b for b in locus.bubbles if len(b.snps) > 1]
        self.n_real = sum(1 for b in locus.bubbles if b.parent is None)
        self.n_synthetic = len(locus.bubbles) - self.n_real


def bubble_depth_vector(coordinate, carriers):
    """Read depth of one allele of one site in every sample, from the kissreads2 counts.

    Several paths may carry the allele: the pairwise bubbles of a multi-allelic
    site, the same SNP in several contexts...  A read may be counted by several
    of them, so they cannot simply be added.  Two paths are counted apart only
    when they differ at a SNP position of one of the two bubbles, because
    kissreads2 refuses any mismatch there: a read cannot map on both.  Inside a
    group the largest count is kept.
    """
    groups = []
    for bubble, path_index in carriers:
        own = {bubble.coordinate(p) for p, _, _ in bubble.snps}
        placed = False
        for group in groups:
            for other, other_index, other_own in group:
                differ = False
                for c in (own | other_own) - {coordinate}:
                    n1, n2 = bubble.nucleotide(path_index, c), other.nucleotide(other_index, c)
                    if n1 is not None and n2 is not None and n1 != n2:
                        differ = True
                        break
                if not differ:
                    group.append((bubble, path_index, own))
                    placed = True
                    break
            if placed:
                break
        if not placed:
            groups.append([(bubble, path_index, own)])
    total = None
    for group in groups:
        stacked = np.stack([b.count_vector(i) for b, i, _ in group]).max(axis=0).astype(np.int64)
        total = stacked if total is None else total + stacked
    return np.minimum(total, 65535).astype(np.uint16)


def fact_fragments(parts, materialised, site_specific_of):
    """One read (pair) -> {LocusData: {coordinate: nucleotide}}: at most ONE observation per site."""
    per_locus = {}
    n_conflicts = 0
    for part in parts:
        for _, bubble_id, path_index, _ in part:
            index = materialised.get(bubble_id)
            if index is None:
                continue
            bubble = index
            data = bubble.locus.name
            observed = per_locus.setdefault(data, {})
            site_specific = site_specific_of[data]
            for position, _, _ in bubble.snps:
                coordinate = bubble.coordinate(position)
                if len(bubble.snps) > 1 and coordinate in site_specific:
                    # kissreads2 reports a multi-SNP path even when the read covers one of
                    # its SNPs only: not an observation of this site, a single-SNP bubble is
                    continue
                nucleotide = bubble.nucleotide(path_index, coordinate)
                if observed.get(coordinate, nucleotide) != nucleotide:
                    nucleotide = None          # two bubbles disagree: drop this site
                    n_conflicts += 1
                observed[coordinate] = nucleotide
    return per_locus, n_conflicts


FRAGMENT_CACHE_SIZE = 2_000_000


def cached_fragments(text, cache, by_id, site_specific_of, stats, site_facts=False):
    """Observations of one fact text, cached: the same texts recur in every read set.

    The cached value is a tuple of (locus name, ((coordinate, nucleotide), ...)).
    site_facts: the text comes from a phased_sites file (kissreads2 -phasing_sites).
    """
    hit = cache.get(text)
    if hit is not None:
        return hit
    if site_facts:
        per_locus, conflicts = site_fact_observations(text, by_id)
    else:
        parts = tuple(tuple(p) for p in map(parse_fact_part, text.split()) if p)
        per_locus, conflicts = fact_fragments(parts, by_id, site_specific_of)
    stats[0] += conflicts
    value = tuple((name, tuple((c, n) for c, n in observed.items() if n is not None))
                  for name, observed in per_locus.items())
    value = tuple(v for v in value if v[1])
    if len(cache) < FRAGMENT_CACHE_SIZE:
        cache[text] = value
    return value


def call_multi_loci(loci, materialised, store, fact_files, args, vcf, hap_file, locus_file, fasta, site_files=()):
    """Everything for the loci with several bubbles.  Returns the number of loci written."""
    data = [LocusData(locus) for locus in loci]
    site_specific_of = {d.locus.name: d.site_specific for d in data}
    by_id = {b.id: b for d in data for b in d.locus.bubbles}   # the bubbles of the loci kept (--min_sites)
    S = store.nsamples

    # ---- global slots
    site_base = sa_base = pl_base = 0
    for d in data:
        d.site_base, d.sa_base, d.pl_base = site_base, sa_base, pl_base
        for coordinate in d.coordinates:
            A = len(d.alleles[coordinate])
            site_base, sa_base, pl_base = site_base + 1, sa_base + A, pl_base + A * (A + 1) // 2
    n_sites, n_sa, n_pl = site_base, sa_base, pl_base
    AD = np.zeros((n_sa, S), dtype=np.uint16)
    GT = np.full((n_sites, S, 2), -1, dtype=np.int8)
    PH = np.zeros((n_sites, S), dtype=bool)
    PS = np.zeros((n_sites, S), dtype=np.int32)
    GQ = np.zeros((n_sites, S), dtype=np.uint8)
    PL = np.zeros((n_pl, S), dtype=np.uint16)
    site_n_alleles = np.zeros(n_sites, dtype=np.int8)
    site_sa = np.zeros(n_sites, dtype=np.int64)
    site_pl = np.zeros(n_sites, dtype=np.int64)
    sa_index = {}                                       # (locus name, coordinate, nt) -> row of AD
    for d in data:
        site = d.site_base
        sa = d.sa_base
        pl = d.pl_base
        for coordinate in d.coordinates:
            alleles = d.alleles[coordinate]
            site_n_alleles[site], site_sa[site], site_pl[site] = len(alleles), sa, pl
            # the counts of a multi-SNP bubble are not site specific (see sub_bubble):
            # they are only used for an allele that no single-SNP bubble carries
            carriers = d.locus.sites[coordinate]
            for nucleotide in alleles:
                counted = [c for c in carriers[nucleotide] if len(c[0].snps) == 1] or carriers[nucleotide]
                AD[sa] = bubble_depth_vector(coordinate, counted)
                sa_index[(d.locus.name, coordinate, nucleotide)] = sa
                sa += 1
            site += 1
            pl += len(alleles) * (len(alleles) + 1) // 2
    sites_by_A = {A: np.nonzero(site_n_alleles == A)[0] for A in (2, 3, 4)}
    fact_of_sample = {sample_of_fact_file(f): (f, False) for f in fact_files}
    # the phased_sites files, when present, replace the classic facts as the source of observations
    fact_of_sample.update({sample_of_fact_file(f): (f, True) for f in site_files})
    names = {sample: read_set_name(f) for sample, (f, _) in fact_of_sample.items()}
    stats = [0]                                          # contradictory observations in facts
    cache, classic_cache = {}, {}
    locus_site_base = np.array([d.site_base for d in data] + [n_sites], dtype=np.int64)
    locus_of_site = np.repeat(np.arange(len(data)), np.diff(locus_site_base))
    site_coordinate = np.array([c for d in data for c in d.coordinates], dtype=np.int64)

    # ---- sample by sample: facts -> depths, genotypes, phasing
    for sample in range(S):
        votes = {}                        # locus name -> Counter((c1, c2, n1, n2))
        fact_file, site_facts = fact_of_sample.get(sample, (None, False))
        if fact_file is not None:
            fact_depth = Counter()
            with open(fact_file) as handle:
                for line in handle:
                    if line.startswith("#") or "=>" not in line:
                        continue
                    text, support = line.rsplit("=>", 1)
                    support = int(support)
                    for name, observed in cached_fragments(text, cache if site_facts else classic_cache, by_id,
                                                           site_specific_of, stats, site_facts):
                        for coordinate, nucleotide in observed:
                            row = sa_index.get((name, coordinate, nucleotide))
                            if row is not None:
                                fact_depth[row] += support
                        add_pair_votes(dict(observed), support, votes.setdefault(name, Counter()))
            if fact_depth:
                rows = np.fromiter(fact_depth.keys(), dtype=np.int64, count=len(fact_depth))
                values = np.fromiter(fact_depth.values(), dtype=np.int64, count=len(fact_depth))
                AD[rows, sample] = np.maximum(AD[rows, sample], np.minimum(values, 65535))
            del fact_depth
        # genotypes of every site of this sample, vectorised by number of alleles
        for A, sites in sites_by_A.items():
            if not len(sites):
                continue
            depths = AD[site_sa[sites][:, None] + np.arange(A)[None, :], sample]
            gt, pl, gq = call_genotypes(depths, A, args.min_depth)
            GT[sites, sample] = gt
            GQ[sites, sample] = gq
            PL[site_pl[sites][:, None] + np.arange(pl.shape[1])[None, :], sample] = pl
        # phasing: only the loci where this sample has two or more heterozygous sites need it;
        # elsewhere every called site is trivially phased, phase set = first called site
        gt_sample = GT[:, sample]
        called_sites = gt_sample[:, 0] >= 0
        het_sites = called_sites & (gt_sample[:, 0] != gt_sample[:, 1])
        n_het = np.add.reduceat(het_sites.astype(np.int64), locus_site_base[:-1])
        n_het[np.diff(locus_site_base) == 0] = 0
        first_called = np.minimum.reduceat(np.where(called_sites, site_coordinate, 1 << 40), locus_site_base[:-1])
        trivial = (n_het < 2)[locus_of_site] & called_sites
        PH[trivial, sample] = True
        PS[trivial, sample] = first_called[locus_of_site[trivial]] + 1
        for locus_index in np.nonzero(n_het >= 2)[0].tolist():
            d = data[locus_index]
            called = {}
            for i, coordinate in enumerate(d.coordinates):
                a, b = gt_sample[d.site_base + i]
                if a >= 0:
                    alleles = d.alleles[coordinate]
                    called[coordinate] = (alleles[a], alleles[b])
            if not called:
                continue
            pair_votes = votes.get(d.locus.name)
            if pair_votes is None:
                # without facts (lone bubble, no -phasing) the two paths of a multi-SNP bubble
                # are the only phase information
                pair_votes = Counter()
                for bubble in d.multi_snp:
                    for path_index in (0, 1):
                        observed = {bubble.coordinate(p): bubble.nucleotide(path_index, bubble.coordinate(p))
                                    for p, _, _ in bubble.snps}
                        add_pair_votes(observed, bubble.count(path_index, sample), pair_votes)
            haplotypes, phased = phase_from_votes(called, pair_votes, args.min_support)
            phase_set = min(phased) + 1 if phased else 0
            for i, coordinate in enumerate(d.coordinates):
                if coordinate in haplotypes:
                    alleles = d.alleles[coordinate]
                    site = d.site_base + i
                    GT[site, sample] = (alleles.index(haplotypes[coordinate][0]),
                                        alleles.index(haplotypes[coordinate][1]))
                    if coordinate in phased:
                        PH[site, sample] = True
                        PS[site, sample] = phase_set
        del votes
    n_conflicts = stats[0]

    # ---- output
    n_multiallelic = 0
    for number, d in enumerate(data, 1):
        locus = d.locus
        locus.name = f"locus_{number}"
        for i, coordinate in enumerate(d.coordinates):
            site = d.site_base + i
            alleles = d.alleles[coordinate]
            A = len(alleles)
            n_multiallelic += A > 2
            involved = {b for cs in locus.sites[coordinate].values() for b, _ in cs}
            real = sorted(b.id for b in involved if b.parent is None)
            ranks = [b.rank for b in involved if b.rank is not None]
            info = (f"NB={len(real)};NX={len(involved) - len(real)};BUB={','.join(map(str, real))};"
                    f"RK={max(ranks) if ranks else '.'};NSITES={len(d.coordinates)};PC={locus.conflicts}")
            columns = [locus.name, str(coordinate + 1), f"{locus.name}_{coordinate + 1}", alleles[0],
                       ",".join(alleles[1:]), ".", ".", info, "GT:PS:DP:AD:GQ:PL"]
            ad = AD[site_sa[site]:site_sa[site] + A].T.tolist()
            pl = PL[site_pl[site]:site_pl[site] + A * (A + 1) // 2].T.tolist()
            gt, ph, ps, gq = GT[site].tolist(), PH[site].tolist(), PS[site].tolist(), GQ[site].tolist()
            columns.extend(format_field(gt[s], ph[s], ps[s], ad[s], gq[s], pl[s]) for s in range(S))
            vcf.write("\t".join(columns) + "\n")
        seen = Counter()
        positions = ",".join(str(c + 1) for c in d.coordinates)
        for sample in range(S):
            hap_1, hap_2, n_missing, n_unphased = "", "", 0, 0
            for i, coordinate in enumerate(d.coordinates):
                site = d.site_base + i
                a, b = GT[site, sample]
                if a < 0:
                    hap_1, hap_2, n_missing = hap_1 + "N", hap_2 + "N", n_missing + 1
                elif not PH[site, sample]:
                    hap_1, hap_2, n_unphased = hap_1 + "?", hap_2 + "?", n_unphased + 1
                else:
                    hap_1 += d.alleles[coordinate][a]
                    hap_2 += d.alleles[coordinate][b]
            if n_missing == len(d.coordinates):
                status = "missing"
            elif n_missing or n_unphased:
                status = "partial"
            else:
                status = "resolved"
                hap_1, hap_2 = sorted((hap_1, hap_2))
                seen.update((hap_1, hap_2))
            hap_file.write("\t".join(map(str, (locus.name, len(d.coordinates), positions, f"G{sample + 1}",
                                               names.get(sample, "."), status, hap_1, hap_2,
                                               n_missing, n_unphased))) + "\n")
        locus_file.write("\t".join(map(str, (
            locus.name, locus.length, d.n_real, d.n_synthetic, len(d.coordinates),
            max(len(a) for a in d.alleles.values()), len(seen),
            ",".join(f"{h}:{n}" for h, n in sorted(seen.items())) or ".", locus.conflicts))) + "\n")
        fasta.write(f">{locus.name}\n{locus_sequence(locus, {c: a[0] for c, a in d.alleles.items()})}\n")
    log(f"[call] {len(data)} loci with several bubbles: {n_sites} sites, {n_multiallelic} with more than two alleles,"
        f" {sum(l.conflicts for l in loci)} placement conflicts, {n_conflicts} contradictory observations in facts")
    return len(data)


###############################################################################
# call: isolated bubbles (vectorised)
###############################################################################

def call_single_bubbles(store, indices, first_number, args, vcf, hap_file, names, locus_file, fasta):
    """Loci made of one bubble: biallelic sites sharing the bubble's read counts.

    A lone multi-SNP bubble is phased by construction (its two paths): it is
    written to the haplotype table like the loci with several bubbles.
    """
    S = store.nsamples
    number = first_number
    n_sites = 0
    field_cache = {}
    for i0 in range(0, len(indices), CHUNK_ROWS // 2):
        chunk = indices[i0:i0 + CHUNK_ROWS // 2]
        depths = np.stack([store.counts[2 * chunk], store.counts[2 * chunk + 1]], axis=2)   # [m, S, 2]
        gt, pl, gq = call_genotypes(depths.reshape(-1, 2), 2, args.min_depth)
        gt, pl, gq = gt.reshape(len(chunk), S, 2), pl.reshape(len(chunk), S, 3), gq.reshape(len(chunk), S)
        # a biallelic field is a function of (phase set, depth higher, depth lower): format each once
        keys = (depths[:, :, 0].astype(np.int64) << 16) | depths[:, :, 1].astype(np.int64)
        snps = store.snp_positions(chunk)
        ranks = store.ranks[chunk].tolist()
        ids = store.ids[chunk].tolist()
        lens = store.lens[2 * chunk].tolist()
        for j, index in enumerate(chunk.tolist()):
            positions = snps[j].tolist()
            if not positions:
                continue
            higher = store.path_str(2 * index)
            lower = store.path_str(2 * index + 1)
            name = f"locus_{number}"
            number += 1
            rank = "." if math.isnan(ranks[j]) else ranks[j]
            gt_j = gt[j].tolist()
            table = field_cache.setdefault(positions[0] + 1, {})
            fields = []
            for s, key in enumerate(keys[j].tolist()):
                field = table.get(key)
                if field is None:
                    field = format_field(gt_j[s], True, positions[0] + 1, depths[j, s].tolist(),
                                         int(gq[j, s]), pl[j, s].tolist())
                    if len(table) < 1_000_000:
                        table[key] = field
                fields.append(field)
            sample_columns = "\t".join(fields)
            for position in positions:
                info = f"NB=1;NX=0;BUB={ids[j]};RK={rank};NSITES={len(positions)};PC=0"
                vcf.write(f"{name}\t{position + 1}\t{name}_{position + 1}\t{higher[position]}\t{lower[position]}"
                          f"\t.\t.\t{info}\tGT:PS:DP:AD:GQ:PL\t{sample_columns}\n")
            n_sites += len(positions)
            hap_h = "".join(higher[p] for p in positions)
            hap_l = "".join(lower[p] for p in positions)
            called = gt_j and [g for g in gt_j if g[0] >= 0]
            n_h = sum(g.count(0) for g in called)      # haplotype copies, as for the other loci
            n_l = sum(g.count(1) for g in called)
            haplotypes = [f"{h}:{n}" for h, n in ((hap_h, n_h), (hap_l, n_l)) if n]
            if len(positions) > 1:
                pos_text = ",".join(str(p + 1) for p in positions)
                for s, (a, b) in enumerate(gt_j):
                    if a < 0:
                        row = ("missing", "N" * len(positions), "N" * len(positions), len(positions))
                    else:
                        row = ("resolved", (hap_h, hap_l)[a], (hap_h, hap_l)[b], 0)
                    hap_file.write("\t".join(map(str, (name, len(positions), pos_text, f"G{s + 1}",
                                                       names.get(s, "."), *row, 0))) + "\n")
            locus_file.write("\t".join(map(str, (name, lens[j], 1, 0, len(positions), 2, len(haplotypes),
                                                 ",".join(sorted(haplotypes)) or ".", 0))) + "\n")
            fasta.write(f">{name}\n{higher}\n")
    log(f"[call] {number - first_number} isolated bubbles: {n_sites} sites")
    return number - first_number


###############################################################################
# call: main
###############################################################################

VCF_HEADER = """##fileformat=VCFv4.2
##source=disco_haplotypes.py
##INFO=<ID=NB,Number=1,Type=Integer,Description="Number of kissnp2 bubbles describing this site">
##INFO=<ID=NX,Number=1,Type=Integer,Description="Number of synthetic context bubbles describing this site">
##INFO=<ID=BUB,Number=.,Type=String,Description="Ids of the kissnp2 bubbles describing this site">
##INFO=<ID=RK,Number=1,Type=Float,Description="Best DiscoSnp rank among these bubbles">
##INFO=<ID=NSITES,Number=1,Type=Integer,Description="Number of sites of the locus">
##INFO=<ID=PC,Number=1,Type=Integer,Description="Bubble placement conflicts in this locus (repeats, paralogs)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; | = phased by the reads inside the locus">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set (position of its first site)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sum of the allele depths">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths (lower bounds, see the documentation)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
"""


def call(args):
    parents = read_synthetic_map(args.map)
    synthetic_ids = np.array(sorted(parents), dtype=np.int64)
    store = parse_store(args.coherent, with_counts=True)
    if args.uncoherent and len(synthetic_ids) and os.path.exists(args.uncoherent):
        # a context existing for one allele only is 'uncoherent' for kissreads2,
        # its read counts are nevertheless valid
        extra = parse_store(args.uncoherent, with_counts=True, only_ids=synthetic_ids)
        if extra.n:
            keep = np.array([store.index_of(int(i)) < 0 for i in extra.ids], dtype=bool)
            if keep.any():
                if extra.nsamples != store.nsamples:
                    sys.exit("ERROR: the coherent and uncoherent files have different numbers of read sets")
                width = max(store.max_len, extra.max_len)
                P = np.full((2 * (store.n + int(keep.sum())), width), PAD, dtype=np.uint8)
                P[:2 * store.n, :store.max_len] = store.P
                rows = np.repeat(2 * np.nonzero(keep)[0], 2) + np.tile([0, 1], int(keep.sum()))
                P[2 * store.n:, :extra.max_len] = extra.P[rows]
                store = BubbleStore(np.concatenate([store.ids, extra.ids[keep]]), P,
                                    np.concatenate([store.lens, extra.lens[rows]]),
                                    np.concatenate([store.counts, extra.counts[rows]]),
                                    np.concatenate([store.ranks, extra.ranks[keep]]))
        del extra
    if not parents:
        log("[call] WARNING: no synthetic bubbles (no 'augment' step before kissreads2). Reads carrying\n"
            "[call]          other alleles at close SNPs are under-counted by kissreads2: genotypes and\n"
            "[call]          haplotypes of loci with several SNPs in one k-window are not reliable.")
    fact_files = [f for f in args.phased if os.path.exists(f)]
    site_files = [f for f in args.sites if os.path.exists(f)]
    if not site_files:
        log("[call] WARNING: no phased_sites files (kissreads2 -phasing_sites): reads mapping a multi-SNP bubble report\n"
            "[call]          only some of its sites, so some heterozygous sites may stay unphased.")
    if not fact_files:
        log("[call] WARNING: no phased facts (kissreads2 -phasing): sites of different bubbles stay unphased.")

    # ---- edges: synthetic -> parent, sequence overlaps, reads
    parent_edges = [[], [], [], []]
    for synthetic_id, (parent_id, start) in parents.items():
        s, p = store.index_of(synthetic_id), store.index_of(parent_id)
        if s >= 0 and p >= 0:
            b1, b2, shift = (p, s, start) if p < s else (s, p, -start)
            parent_edges[0].append(b1)
            parent_edges[1].append(b2)
            parent_edges[2].append(shift)
            parent_edges[3].append(1)
    edges = [np.array(e, dtype=np.int64) for e in parent_edges]
    for source in (sequence_edges(store, args.seed_size, args.max_mismatches, args.min_overlap),
                   fact_edges(fact_files, store)):
        edges = [np.concatenate([a, b]) for a, b in zip(edges, source)]
    edges = tuple(edges)
    loci, materialised = multi_bubble_loci(store, edges, args.max_locus_bubbles)
    for bubble in materialised.values():
        if bubble.id in parents:
            bubble.parent = parents[bubble.id][0]
    for number, locus in enumerate(loci, 1):
        locus.name = f"locus_{number}"
    in_locus = np.zeros(store.n, dtype=bool)
    in_locus[list(materialised)] = True
    is_synthetic = np.isin(store.ids, synthetic_ids)
    singles = np.nonzero(~in_locus & ~is_synthetic)[0]
    n_dropped = int((~in_locus & is_synthetic).sum())
    if n_dropped:
        log(f"[call] {n_dropped} synthetic bubbles without their parent bubble were dropped")
    if args.min_sites > 1:
        loci = [locus for locus in loci if len(locus.sites) >= args.min_sites]
        singles = singles[:0]
    del edges

    with open(args.out + ".vcf", "w") as vcf, \
            open(args.out + ".tsv", "w") as hap_file, \
            open(args.out + "_loci.tsv", "w") as locus_file, \
            open(args.out + "_loci.fa", "w") as fasta:
        vcf.write(VCF_HEADER)
        # contig lines: multi-bubble loci first, then isolated bubbles, in the order the records follow
        for number, locus in enumerate(loci, 1):
            vcf.write(f"##contig=<ID=locus_{number},length={locus.length}>\n")
        single_lens = store.lens[2 * singles].tolist()
        snp_counts = np.zeros(len(singles), dtype=np.int64)
        for i0 in range(0, len(singles), CHUNK_ROWS):
            snp_counts[i0:i0 + CHUNK_ROWS] = [len(p) for p in store.snp_positions(singles[i0:i0 + CHUNK_ROWS])]
        number = len(loci) + 1
        for length, n_snp in zip(single_lens, snp_counts.tolist()):
            if n_snp:
                vcf.write(f"##contig=<ID=locus_{number},length={length}>\n")
                number += 1
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(f"G{s + 1}" for s in range(store.nsamples)) + "\n")
        hap_file.write("locus\tn_sites\tpositions\tsample\tread_set\tstatus\thaplotype_1\thaplotype_2"
                       "\tn_missing\tn_unphased\n")
        locus_file.write("locus\tlength\tn_bubbles\tn_synthetic\tn_sites\tmax_alleles\tn_haplotypes"
                         "\thaplotypes\tplacement_conflicts\n")
        n_multi = call_multi_loci(loci, materialised, store, fact_files, args, vcf, hap_file, locus_file, fasta,
                                  site_files)
        if len(singles):
            names = {sample_of_fact_file(f): read_set_name(f) for f in fact_files + site_files}
            call_single_bubbles(store, singles, n_multi + 1, args, vcf, hap_file, names, locus_file, fasta)


###############################################################################
# strip
###############################################################################

def strip(args):
    synthetic = set(read_synthetic_map(args.map))
    keep = True
    with open(args.input, "rb") as handle, open(args.output, "wb") as out:
        for line in handle:
            if line[:1] == b">":
                match = re.match(rb">SNP_\w+_path_(\d+)", line)
                keep = match is None or int(match.group(1)) not in synthetic
            if keep:
                out.write(line)


###############################################################################

def add_placement_options(parser):
    parser.add_argument("--seed_size", type=int, default=16,
                        help="exact seed used to find overlapping bubble paths [16]")
    parser.add_argument("--max_mismatches", type=int, default=4,
                        help="mismatches tolerated in the overlap of two paths (other contexts) [4]")
    parser.add_argument("--min_overlap", type=int, default=25,
                        help="minimal overlap between two paths [25]")
    parser.add_argument("--max_locus_bubbles", type=int, default=500,
                        help="components with more bubbles are repeats: their bubbles are kept isolated [500]")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    sub = commands.add_parser("augment", help="add the missing sequence contexts (before kissreads2)")
    sub.add_argument("-i", "--input", required=True, help="kissnp2 fasta file")
    sub.add_argument("-o", "--output", required=True, help="augmented fasta file")
    sub.add_argument("-m", "--map", required=True, help="output table: synthetic bubble -> parent bubble")
    sub.add_argument("--max_contexts", type=int, default=32,
                     help="maximal number of contexts written for one SNP [32]")
    add_placement_options(sub)
    sub.set_defaults(function=augment)

    sub = commands.add_parser("call", help="sites, genotypes and haplotypes (after kissreads2)")
    sub.add_argument("-c", "--coherent", required=True, help="kissreads2 *_coherent.fa")
    sub.add_argument("-u", "--uncoherent", help="kissreads2 *_uncoherent.fa (synthetic bubbles only are read)")
    sub.add_argument("-m", "--map", help="table written by 'augment'")
    sub.add_argument("-p", "--phased", nargs="*", default=[],
                     help="phased_alleles_read_set_id_*.txt files (kissreads2 -phasing)")
    sub.add_argument("-s", "--sites", nargs="*", default=[],
                     help="phased_sites_read_set_id_*.txt files (kissreads2 -phasing_sites): exact per-read site observations")
    sub.add_argument("-o", "--out", required=True, help="prefix of the output files")
    sub.add_argument("--min_depth", type=int, default=3, help="minimal depth to call a genotype [3]")
    sub.add_argument("--min_support", type=int, default=2,
                     help="minimal |cis - trans| read support to phase two sites [2]")
    sub.add_argument("--min_sites", type=int, default=1,
                     help="write only the loci with at least this many sites (2: skip the isolated bubbles,"
                          " which the usual DiscoSnp VCF already describes) [1]")
    add_placement_options(sub)
    sub.set_defaults(function=call)

    sub = commands.add_parser("strip", help="remove the synthetic bubbles from a fasta file")
    sub.add_argument("-i", "--input", required=True)
    sub.add_argument("-o", "--output", required=True)
    sub.add_argument("-m", "--map", required=True)
    sub.set_defaults(function=strip)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
