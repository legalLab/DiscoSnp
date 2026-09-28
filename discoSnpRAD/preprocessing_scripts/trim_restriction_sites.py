#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
trim_restriction_sites.py

Removes the restriction site remnant at the 5' end of RAD / ddRAD reads before
discoSnpRad (option of discoSnpRAD/run_discoSnpRad.sh, on by default).

In RAD / ddRAD data every read starting at a restriction site begins with the
same few nucleotides (the remnant of the site, e.g. TGCAG for PstI, CGG for
MspI, AATTC for EcoRI), whatever the locus. These nucleotides are not
polymorphic and are shared by all the loci.

Input: the file of files (fof) given to discoSnpRad. Each line is one sample:
either a read file (fasta/fastq, gzipped or not) or a fof listing the read files
of the sample (typically R1 and R2). The output fof has the same structure and
points to the trimmed files. A file needing no trimming is not copied.

Trimmed length, per file
------------------------
  --trim_r1 / --trim_r2 auto (default): detected on the file itself, see below.
  --trim_r1 / --trim_r2 <int>         : fixed length for the reads 1 / reads 2.

Reads 1 and reads 2 are told apart by the order of the files in a sample fof
(first = R1, second = R2) or by their names (_R1 _R2, _1. _2., .1. .2. ...).
Files that cannot be told apart are treated as reads 1.

Automatic detection
-------------------
Every read starting at a restriction site begins with its remnant, but a
library also holds reads that do not (adapter dimers, contaminants, repeated
artefacts, organelle reads...): commonly 5 to 20 %, sometimes more. The
detection therefore does not ask every read to carry the site. The first
--n_reads reads of the file (reads from thousands of loci) are profiled
position by position from the 5' end:
  - N are ignored: the frequencies are computed on A, C, G and T only. A
    position where more than half of the reads have an N (a failed sequencing
    cycle, often the second one) is 'unknown': it neither extends nor stops the
    site, but is trimmed when site positions follow it.
  - the first position starts a site when one nucleotide makes at least
    --min_site_fraction of the reads (default 0.5): the 'level' of the site,
    about the fraction of the reads carrying it.
  - a next position continues the site while its major nucleotide stays at
    --plateau x the level of the site (default 0.9). Two nucleotides together
    are accepted (degenerate sites such as ApeKI G^CWGC) when the second one
    makes at least 30 % of the pair. Behind the remnant the reads enter the
    loci and the frequency drops (e.g. 0.82 -> 0.41, or 1.00 -> 0.84 in a low
    diversity library): the site ends there.
  - the trimmed length is the end of the site, if it holds at least 2 positions.
Then, since all the files of a RAD / ddRAD run share their enzymes, the site
found in most files of the same read (reads 1, or reads 2) is looked for at the
5' end of the reads of the files where nothing was detected: when it starts at
least --rescue_fraction of their reads (default 0.1; by chance a 3 to 5 nt site
starts 1.6 to 0.1 % of the reads), it is trimmed from them too. A file where a
third of the reads or more carry the site is found by the plateau rule itself.
The detection is done per file, so that
  - single digest RAD: reads 2 start at a random (sheared) position, nothing
    is conserved and nothing is trimmed,
  - an inline barcode left after demultiplexing (constant in a sample) is
    removed with the remnant.
Profiling more reads does not change the result: the fractions are stable
along a file. At most --n_reads reads are read per file for the detection.
Limits: variable length spacers ("heterogeneity spacers", staggered adapters)
shift the remnant from read to read: nothing is conserved, the file is not
trimmed and a warning is printed. Remove the spacers first (e.g. with the
demultiplexing tool) or give the lengths with --trim_r1 / --trim_r2.

A report (trimming_report.tsv in the output directory) gives, per file, the
trimmed length, the consensus of the trimmed nucleotides and the profile of the
first positions.

Reads without the site (default: removed)
-----------------------------------------
Reads that do not start with the site of their file (adapter dimers,
contaminants, repeated artefacts, organelle reads...) are removed. The site is
the detected one (with --trim_r1 / --trim_r2 <int>, the detected site within
the given length; a position without site matches anything, and a file where
no site is detected is not filtered); N in the reads and degenerate positions of the site
(e.g. W of ApeKI) match. The files of a sample fof (typically R1 and R2) are
read together, as pairs: a pair is removed when one of its reads lacks the site
of its file, and both files are rewritten so that they stay paired (the reads 2
of a single digest RAD have no site and are not checked, but lose the pairs of
the removed reads 1). Files without a site are not filtered. With
--keep_reads_without_site nothing is removed: the reads 1 and 2 of a pair stay
in the same order, and a read not longer than the trimmed length is replaced by
a single N (this is also the case for the reads kept by the filter).
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S", frozenset("AT"): "W",
         frozenset("GT"): "K", frozenset("AC"): "M"}
SEQUENCE_EXTENSIONS = re.compile(r"(\.(fastq|fq|fasta|fa|fna|txt))?(\.gz)?$", re.IGNORECASE)
MATE_RE = re.compile(r"(?:^|[._-])R?([12])(?=[._-]|$)", re.IGNORECASE)
MIN_COUNT = 100              # nucleotides (A, C, G, T) needed to judge a position


def log(message):
    sys.stderr.write(message + "\n")
    sys.stderr.flush()


###############################################################################
# files
###############################################################################

def is_gzip(path):
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


class Reader:
    """Binary line reader of a (gzipped) file, decompressed by pigz / gzip when available."""

    def __init__(self, path):
        self.process = None
        if is_gzip(path):
            tool = shutil.which("pigz") or shutil.which("gzip")
            if tool:
                self.process = subprocess.Popen([tool, "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 20)
                self.handle = self.process.stdout
            else:
                self.handle = gzip.open(path, "rb")
        else:
            self.handle = open(path, "rb")

    def close(self):
        self.handle.close()
        if self.process is not None:
            self.process.kill()
            self.process.wait()


class Writer:
    """Binary writer of a gzipped file, compressed by pigz / gzip when available."""

    def __init__(self, path, threads=2):
        self.process = None
        self.output = open(path, "wb")
        pigz, gzip_tool = shutil.which("pigz"), shutil.which("gzip")
        command = [pigz, "-1", "-p", str(threads)] if pigz else ([gzip_tool, "-1"] if gzip_tool else None)
        if command:
            self.process = subprocess.Popen(command + ["-c"], stdin=subprocess.PIPE, stdout=self.output,
                                            bufsize=1 << 20)
            self.handle = self.process.stdin
        else:
            self.handle = gzip.open(self.output, "wb", compresslevel=1)

    def close(self):
        self.handle.close()
        if self.process is not None:
            if self.process.wait() != 0:
                raise IOError("compression failed")
        self.output.close()


def records(handle):
    """(header, sequence, quality or None) of a fasta or fastq file (bytes, no end of line)."""
    line = handle.readline()
    while line and not line.strip():
        line = handle.readline()
    if not line:
        return
    if line[:1] == b"@":                                 # fastq: 4 lines per record
        while line:
            if line.strip():
                sequence = handle.readline().rstrip(b"\r\n")
                handle.readline()
                quality = handle.readline().rstrip(b"\r\n")
                yield line.rstrip(b"\r\n"), sequence, quality
            line = handle.readline()
    elif line[:1] == b">":                               # fasta, possibly multi-line
        header, parts = line.rstrip(b"\r\n"), []
        for line in handle:
            if line[:1] == b">":
                yield header, b"".join(parts), None
                header, parts = line.rstrip(b"\r\n"), []
            else:
                parts.append(line.strip())
        yield header, b"".join(parts), None
    else:
        raise ValueError("not a fasta or fastq file")


def is_sequence_file(path):
    if is_gzip(path):
        return True
    with open(path, "rb") as handle:
        for line in handle:
            if line.strip():
                return line[:1] in (b">", b"@")
    return False


def resolve(path, fof):
    """A path of a fof: as given (relative to the working directory, as discoSnpRad reads it),
    else relative to the directory of the fof."""
    path = path.strip()
    if os.path.exists(path):
        return os.path.abspath(path)
    other = os.path.join(os.path.dirname(os.path.abspath(fof)), path)
    if os.path.exists(other):
        return other
    sys.exit(f"ERROR: file {path} (listed in {fof}) does not exist")


def mate_from_name(path):
    stem = SEQUENCE_EXTENSIONS.sub("", os.path.basename(path))
    matches = MATE_RE.findall(stem)
    return int(matches[-1]) if matches else None


def read_fof(fof):
    """[(sample fof or None, [(file, mate), ...]), ...], one entry per line of the fof."""
    samples = []
    with open(fof) as handle:
        lines = [line.strip() for line in handle if line.strip()]
    for line in lines:
        path = resolve(line, fof)
        if is_sequence_file(path):
            samples.append((None, [(path, mate_from_name(path) or 1)]))
            continue
        with open(path) as handle:
            files = [resolve(l, path) for l in handle if l.strip()]
        mates = []
        for position, file_path in enumerate(files):
            by_name = mate_from_name(file_path)
            by_order = position + 1 if len(files) == 2 else None
            if by_name and by_order and by_name != by_order:
                log(f"WARNING: {file_path} is file {by_order} of {path} but its name says R{by_name}: R{by_name} used")
            mates.append((file_path, by_name or by_order or 1))
        samples.append((path, mates))
    return samples


###############################################################################
# detection
###############################################################################

MOTIF = 4                    # staggered remnants: 4-mers looked for in the first MOTIF_WINDOW nucleotides
MOTIF_WINDOW = 12


def profile(path, n_reads, max_trim):
    """For the first n_reads reads: nucleotide counts of the first max_trim + 1 positions, the 5'
    ends (first max_trim nucleotides) and the number of reads holding each 4-mer in their first
    12 nucleotides (staggered remnants)."""
    counts = [Counter() for _ in range(max_trim + 1)]
    starts = Counter()
    motifs = Counter()
    reader = Reader(path)
    n = 0
    try:
        for _, sequence, _ in records(reader.handle):
            start = sequence[:max(max_trim + 1, MOTIF_WINDOW)].upper()
            for position, nucleotide in enumerate(start[:max_trim + 1]):
                counts[position][nucleotide] += 1
            starts[start[:max_trim]] += 1
            window = start[:MOTIF_WINDOW]
            motifs.update({window[i:i + MOTIF] for i in range(len(window) - MOTIF + 1)} - {b""})
            n += 1
            if n >= n_reads:
                break
    finally:
        reader.close()
    return counts, starts, motifs, n


def position_frequencies(counter):
    """(major nucleotide, its frequency, second nucleotide, its frequency, frequency of the two others)
    on A, C, G, T, or None when the position is unknown (mostly N, or too few nucleotides)."""
    total = sum(counter.values())
    acgt = sorted(((counter.get(ord(b), 0), b) for b in "ACGT"), reverse=True)
    informative = sum(c for c, _ in acgt)
    if not total or counter.get(ord("N"), 0) > total / 2 or informative < MIN_COUNT:
        return None
    return (acgt[0][1], acgt[0][0] / informative, acgt[1][1], acgt[1][0] / informative,
            (acgt[2][0] + acgt[3][0]) / informative)


def site_position(frequencies, threshold):
    """(IUPAC code, frequency) if the position reaches the threshold with one nucleotide, or with two
    (degenerate site: the second one makes at least 30 % of the pair and the two others together at
    most 15 % of it; a random position has its four nucleotides at ~25 % each); else None."""
    top, f_top, second, f_second, f_others = frequencies
    if f_top >= threshold:
        return top, f_top
    pair = f_top + f_second
    if pair >= threshold and f_second >= 0.3 * pair and f_others <= 0.15 * pair:
        return IUPAC[frozenset((top, second))], f_top + f_second
    return None


def detect(counts, args):
    """(trimmed length, consensus of the trimmed nucleotides, profile text)"""
    level = None                    # frequency of the site (about the fraction of the reads carrying it)
    codes, trim, n_site = [], 0, 0
    for counter in counts:
        frequencies = position_frequencies(counter)
        if frequencies is None:                          # failed cycle: transparent
            codes.append("N")
            continue
        position = site_position(frequencies, args.min_site_fraction if level is None else args.plateau * level)
        if position is None:
            break
        code, frequency = position
        codes.append(code)
        level = frequency if level is None else max(level, frequency)
        n_site += 1
        trim = len(codes)
    if n_site < 2:
        trim = 0
    consensus = "".join(codes[:trim])
    profile_text = " ".join("N:-" if f is None else f"{f[0]}:{f[1]:.2f}"
                            for f in map(position_frequencies, counts))
    return trim, consensus, profile_text


def matching_fraction(starts, consensus, n_reads):
    """Fraction of the reads whose 5' end matches the consensus (IUPAC codes; N match anything)."""
    allowed = []
    for code in consensus:
        if code == "N":
            allowed.append(None)
        else:
            bases = next((set(k) for k, v in IUPAC.items() if v == code), {code})
            allowed.append({ord(b) for b in bases} | {ord("N")})
    n = 0
    for start, count in starts.items():
        if len(start) >= len(consensus) and all(a is None or c in a for c, a in zip(start, allowed)):
            n += count
    return n / n_reads if n_reads else 0.0


def staggered_warning(motifs, starts, n_reads):
    """A warning when one 4-mer is near the 5' end of most reads but at changing positions."""
    if not motifs or not n_reads:
        return None
    motif, n_with = max(((m, c) for m, c in motifs.items() if b"N" not in m), key=lambda mc: mc[1],
                        default=(b"", 0))
    if n_with <= 0.5 * n_reads:
        return None
    at_start = sum(c for start, c in starts.items() if start.startswith(motif))
    if at_start >= 0.5 * n_with:
        return (f"{motif.decode()} starts {100 * at_start / n_reads:.0f} % of the reads only:"
                f" too few for a site? Nothing trimmed (give --trim_r1 / --trim_r2)")
    return (f"{motif.decode()} is in the first {MOTIF_WINDOW} nt of {100 * n_with / n_reads:.0f} % of the"
            f" reads, at variable positions: variable length spacers before the restriction site?"
            f" Nothing trimmed: remove the spacers first, or give --trim_r1 / --trim_r2")


def detect_file(job):
    """Detection on one read file. Returns a report dict.
    fixed: trimmed length given by the user (manual mode): the site is then the major nucleotides
    of the trimmed positions (used by the filter of the reads without the site)."""
    path, mate, fixed, args = job
    counts, starts, motifs, n = profile(path, args.n_reads, max(args.max_trim, fixed or 0))
    trim, consensus, profile_text = detect(counts, args)
    report = {"file": path, "mate": mate, "mode": "auto", "reads_profiled": n, "trim": trim,
              "consensus": consensus or ".", "profile": profile_text, "warning": None,
              "starts": starts, "motifs": motifs}
    if fixed is not None:
        # the site checked by the filter: the detected site positions within the given length; the other
        # positions (no site detected there, e.g. the sheared reads 2 of a single digest RAD) match anything
        codes = (list(consensus) + ["N"] * fixed)[:fixed]
        site = "".join(codes)
        report.update(mode="manual", trim=fixed, consensus=site if site.strip("N") else ".")
    elif n < MIN_COUNT:
        report.update(trim=0, consensus=".", warning=f"only {n} reads: nothing trimmed")
    return report


def site_pattern(consensus):
    """Compiled pattern of the 5' end of a read carrying the site (N of the read and degenerate codes match)."""
    parts = []
    for code in consensus:
        if code == "N":
            parts.append(b".")
        else:
            bases = next((sorted(k) for k, v in IUPAC.items() if v == code), [code])
            parts.append(b"[" + "".join(bases).encode() + b"Nn]")
    return re.compile(b"".join(parts))


def trim_unit(job):
    """One sample line of the fof: its files (1, or R1 and R2...) read together, as pairs.
    Every file is trimmed; a read (pair) is removed if one of its reads lacks the site of its file.
    Returns (output paths, reads (pairs) read, reads (pairs) removed, reads too short per file)."""
    files, out_dir, threads, filtering = job
    patterns = [site_pattern(consensus) if filtering and trim > 0 and consensus != "." else None
                for _, _, trim, consensus in files]
    rewrite = [trim > 0 or any(patterns) for _, _, trim, _ in files]
    if not any(rewrite):
        return [path for _, path, _, _ in files], 0, 0, [0] * len(files)
    readers = [Reader(path) for _, path, _, _ in files]
    targets, writers = [], []
    try:
        for (index, path, _, _), needed in zip(files, rewrite):
            if needed:
                stem = SEQUENCE_EXTENSIONS.sub("", os.path.basename(path))
                kind = "fastq" if first_character(path) == b"@" else "fasta"
                targets.append(os.path.join(out_dir, f"{index}_{stem}.trimmed.{kind}.gz"))
                writers.append(Writer(targets[-1], threads))
            else:
                targets.append(path)
                writers.append(None)
        iterators = [records(reader.handle) for reader in readers]
        trims = [trim for _, _, trim, _ in files]
        n = n_removed = 0
        n_short = [0] * len(files)
        for reads in zip(*iterators):
            n += 1
            if any(pattern is not None and not pattern.match(sequence)
                   for pattern, (_, sequence, _) in zip(patterns, reads)):
                n_removed += 1
                continue
            for i, ((header, sequence, quality), writer) in enumerate(zip(reads, writers)):
                if writer is None:
                    continue
                sequence = sequence[trims[i]:]
                if not sequence:
                    sequence, n_short[i] = b"N", n_short[i] + 1
                    quality = b"!" if quality is not None else None
                elif quality is not None:
                    quality = quality[trims[i]:]
                if quality is None:
                    writer.handle.write(header + b"\n" + sequence + b"\n")
                else:
                    writer.handle.write(header + b"\n" + sequence + b"\n+\n" + quality + b"\n")
        leftovers = [path for iterator, (_, path, _, _) in zip(iterators, files) if next(iterator, None) is not None]
        if leftovers and len(files) > 1:
            raise ValueError(f"the files of a sample do not have the same number of reads: {', '.join(files[i][1] for i in range(len(files)))}")
    finally:
        for reader in readers:
            reader.close()
        for writer in writers:
            if writer is not None:
                writer.close()
    return targets, n, n_removed, n_short


def rescue_with_consensus(reports, rescue_fraction):
    """Files of a read (R1 or R2) where nothing was detected: trim the site found in most files of
    the same read if it starts at least rescue_fraction of their reads."""
    for mate in sorted({r["mate"] for r in reports}):
        detected = Counter((r["trim"], r["consensus"]) for r in reports
                           if r["mate"] == mate and r["mode"] == "auto" and r["trim"] > 0)
        if not detected:
            continue
        # the most frequent length, then its most frequent consensus without failed cycles
        length = Counter({t: 0 for t, _ in detected})
        for (t, _), n in detected.items():
            length[t] += n
        trim = length.most_common(1)[0][0]
        candidates = sorted(((n, c.count("N") == 0, c) for (t, c), n in detected.items() if t == trim), reverse=True)
        consensus = max(candidates, key=lambda x: (x[1], x[0]))[2]
        for r in reports:
            if r["mate"] != mate or r["mode"] != "auto" or r["trim"] > 0 or not r["reads_profiled"]:
                continue
            fraction = matching_fraction(r["starts"], consensus, r["reads_profiled"])
            if fraction >= rescue_fraction:
                r.update(trim=trim, consensus=consensus, mode=f"auto (R{mate} consensus)",
                         note=f"{consensus} starts {100 * fraction:.0f} % of the reads")
            else:
                r["note"] = f"the R{mate} site {consensus} starts {100 * fraction:.1f} % of the reads only"


def first_character(path):
    reader = Reader(path)
    try:
        for line in reader.handle:
            if line.strip():
                return line[:1]
    finally:
        reader.close()
    return b""


def parse_length(text):
    if text == "auto":
        return None
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("an integer or 'auto' is expected")
    if value < 0:
        raise argparse.ArgumentTypeError("a length cannot be negative")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--fof", required=True, help="file of files given to discoSnpRad")
    parser.add_argument("-o", "--out_dir", required=True, help="directory of the trimmed files")
    parser.add_argument("--out_fof", required=True,
                        help="fof of the trimmed files (written only if at least one file is trimmed)")
    parser.add_argument("--trim_r1", type=parse_length, default=None, metavar="INT|auto",
                        help="nucleotides removed at the 5' end of the reads 1 [auto]")
    parser.add_argument("--trim_r2", type=parse_length, default=None, metavar="INT|auto",
                        help="nucleotides removed at the 5' end of the reads 2 [auto]")
    parser.add_argument("--n_reads", type=int, default=50000, help="reads profiled per file [50000]")
    parser.add_argument("--max_trim", type=int, default=15, help="longest remnant looked for [15]")
    parser.add_argument("--min_site_fraction", type=float, default=0.5,
                        help="fraction of the reads starting with the major nucleotide for a site to start [0.5]")
    parser.add_argument("--plateau", type=float, default=0.9,
                        help="a site position keeps at least this ratio of the frequency of the site [0.9]")
    parser.add_argument("--rescue_fraction", type=float, default=0.1,
                        help="files without a detected site: trim the site of the other files of the same read"
                             " if it starts at least this fraction of the reads [0.1]")
    parser.add_argument("--keep_reads_without_site", action="store_true",
                        help="keep the reads (pairs) not starting with the site of their file (default: removed)")
    parser.add_argument("--threads", type=int, default=0, help="files processed at once [0: all cores]")
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)          # the fofs written list absolute paths
    threads = args.threads or os.cpu_count() or 1
    args.compression_threads = 2

    samples = read_fof(args.fof)
    os.makedirs(args.out_dir, exist_ok=True)
    files = [(f"{sample_index}_{file_index}", path, mate)
             for sample_index, (_, sample_files) in enumerate(samples, 1)
             for file_index, (path, mate) in enumerate(sample_files, 1)]
    workers = min(threads, len(files))

    # 1. detection (auto mode), 2. site shared by the other files of the same read, 3. trimming
    detection = [(path, mate, args.trim_r1 if mate == 1 else args.trim_r2, args) for _, path, mate in files]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        reports = list(pool.map(detect_file, detection))
    rescue_with_consensus(reports, args.rescue_fraction)
    for r in reports:
        if r["mode"] == "auto" and r["trim"] == 0 and r["warning"] is None:
            r["warning"] = staggered_warning(r.get("motifs"), r.get("starts"), r["reads_profiled"])
        r.pop("starts", None)
        r.pop("motifs", None)
        r["output"] = r["file"]
    # one job per sample line: its files are filtered together (pairs)
    units, first = [], 0
    for _, sample_files in samples:
        members = list(range(first, first + len(sample_files)))
        first += len(sample_files)
        units.append((members, ([(files[i][0], reports[i]["file"], reports[i]["trim"], reports[i]["consensus"])
                                 for i in members], args.out_dir, args.compression_threads,
                                not args.keep_reads_without_site)))
    with ProcessPoolExecutor(max_workers=max(1, min(threads, len(units)))) as pool:
        results = list(pool.map(trim_unit, [job for _, job in units]))
    unit_summaries = []
    for (members, _), (targets, n, n_removed, n_short) in zip(units, results):
        for i, target, short in zip(members, targets, n_short):
            reports[i].update(output=target, reads=n, reads_removed=n_removed, reads_too_short=short)
        unit_summaries.append((members, n, n_removed))

    report_file = os.path.join(args.out_dir, "trimming_report.tsv")
    with open(report_file, "w") as handle:
        handle.write("file\tmate\tmode\treads_profiled\ttrimmed\tconsensus\tprofile_of_the_first_positions"
                     "\treads\treads_removed_without_site\ttrimmed_file\tnote\twarning\n")
        for r in reports:
            reads = r.get("reads") or "."
            removed = r.get("reads_removed", ".") if r.get("reads") else "."
            handle.write(f"{r['file']}\tR{r['mate']}\t{r['mode']}\t{r['reads_profiled']}\t{r['trim']}\t"
                         f"{r['consensus']}\t{r['profile']}\t{reads}\t{removed}\t{r['output']}\t"
                         f"{r.get('note') or '.'}\t{r['warning'] or '.'}\n")
    for r in reports:
        text = f"[trim] {os.path.basename(r['file'])} (R{r['mate']}, {r['mode']}): {r['trim']} nt"
        if r["mode"].startswith("auto") or r["trim"]:
            text += f" [{r['consensus']}]" if r["trim"] else " (no conserved 5' end)"
        if r.get("note"):
            text += f" ({r['note']})"
        if r.get("reads_too_short"):
            text += f", {r['reads_too_short']} reads not longer than that replaced by N"
        log(text)
        if r["warning"]:
            log(f"[trim] WARNING {os.path.basename(r['file'])}: {r['warning']}")
    for members, n, n_removed in unit_summaries:
        if not n or args.keep_reads_without_site:
            continue
        what = "pairs" if len(members) > 1 else "reads"
        name = " + ".join(os.path.basename(reports[i]["file"]) for i in members)
        text = f"[filter] {name}: {n_removed} of {n} {what} without the site removed ({100 * n_removed / n:.1f} %)"
        log(text)
        if n_removed > 0.5 * n:
            log(f"[filter] WARNING {name}: more than half of the {what} removed: check the sites in {report_file}"
                f" (or use --keep_reads_without_site)")
    by_mate = {}
    for r in reports:
        by_mate.setdefault(r["mate"], Counter())[r["trim"]] += 1
    for mate, lengths in sorted(by_mate.items()):
        if len(lengths) > 1:
            log(f"[trim] WARNING: the R{mate} files are not trimmed to the same length: "
                + ", ".join(f"{n} files {t} nt" for t, n in sorted(lengths.items()))
                + f" (see {report_file})")

    if all(r["output"] == r["file"] for r in reports):
        if os.path.exists(args.out_fof):
            os.remove(args.out_fof)
        log("[trim] nothing to trim: the input files are used as they are")
        return

    # the new fof, same structure as the input one
    lines = []
    job_index = 0
    for sample_index, (sample_fof, files) in enumerate(samples, 1):
        outputs = []
        for path, _ in files:
            outputs.append(reports[job_index]["output"])
            job_index += 1
        if sample_fof is None:
            lines.append(outputs[0])
        else:
            sample_out = os.path.abspath(os.path.join(args.out_dir, f"{sample_index}_{os.path.basename(sample_fof)}"))
            with open(sample_out, "w") as handle:
                handle.write("\n".join(outputs) + "\n")
            lines.append(sample_out)
    with open(args.out_fof, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    log(f"[trim] trimmed read files listed in {args.out_fof} (report: {report_file})")


if __name__ == "__main__":
    main()
