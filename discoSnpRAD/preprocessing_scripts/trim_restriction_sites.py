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
The first --n_reads reads of the file (reads from thousands of loci) are
profiled position by position from the 5' end:
  - N are ignored: the frequencies are computed on A, C, G and T only. A
    position where more than half of the reads have an N (a failed sequencing
    cycle, often the second one) is 'unknown': it neither extends nor stops the
    conserved prefix, but is trimmed when conserved positions follow it.
  - a position is 'conserved' when one nucleotide makes at least
    --min_fraction of the reads (default 0.9), or two nucleotides make at least
    --min_fraction + (1 - --min_fraction) / 2 of them with at least 20 % each
    (degenerate sites such as ApeKI G^CWGC).
  - the trimmed length is the end of the leading run of conserved (or unknown)
    positions, if it holds at least 2 conserved positions; 0 otherwise.
Behind the remnant the reads enter the loci and no position is conserved any
more: the run stops there. The detection is done per file, so that
  - single digest RAD: reads 2 start at a random (sheared) position, nothing
    is conserved and nothing is trimmed,
  - an inline barcode left after demultiplexing (constant in a sample) is
    removed with the remnant.
Limits: variable length spacers ("heterogeneity spacers", staggered adapters)
shift the remnant from read to read: nothing is conserved, the file is not
trimmed and a warning is printed. Remove the spacers first (e.g. with the
demultiplexing tool) or give the lengths with --trim_r1 / --trim_r2.

A report (trimming_report.tsv in the output directory) gives, per file, the
trimmed length, the consensus of the trimmed nucleotides and the profile of the
first positions.

Reads are never removed (the reads 1 and 2 of a pair stay in the same order);
a read not longer than the trimmed length is replaced by a single N.
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
    """Nucleotide counts of the first max_trim + 1 positions of the first n_reads reads, and the
    number of reads holding each 4-mer in their first 12 nucleotides (staggered remnants)."""
    counts = [Counter() for _ in range(max_trim + 1)]
    motifs = Counter()
    reader = Reader(path)
    n = 0
    try:
        for _, sequence, _ in records(reader.handle):
            start = sequence[:max(max_trim + 1, MOTIF_WINDOW)].upper()
            for position, nucleotide in enumerate(start[:max_trim + 1]):
                counts[position][nucleotide] += 1
            window = start[:MOTIF_WINDOW]
            motifs.update({window[i:i + MOTIF] for i in range(len(window) - MOTIF + 1)} - {b""})
            n += 1
            if n >= n_reads:
                break
    finally:
        reader.close()
    return counts, motifs, n


def position_status(counter, min_fraction):
    """('conserved', code) / ('unknown', 'N') / ('variable', '.') / ('empty', '-'), and the top frequency."""
    total = sum(counter.values())
    if not total:
        return "empty", "-", 0.0
    acgt = sorted(((counter.get(ord(b), 0), b) for b in "ACGT"), reverse=True)
    informative = sum(c for c, _ in acgt)
    if counter.get(ord("N"), 0) > total / 2 or informative < MIN_COUNT:
        return "unknown", "N", 0.0
    top, second = acgt[0][0] / informative, acgt[1][0] / informative
    if top >= min_fraction:
        return "conserved", acgt[0][1], top
    if top + second >= min_fraction + (1 - min_fraction) / 2 and second >= 0.2:
        return "conserved", IUPAC[frozenset((acgt[0][1], acgt[1][1]))], top + second
    return "variable", ".", top


def detect(counts, min_fraction, motifs=None, n_reads=0):
    """(trimmed length, consensus of the trimmed nucleotides, profile text, warning or None)"""
    statuses = [position_status(c, min_fraction) for c in counts]
    last_conserved, n_conserved = -1, 0
    for position, (status, _, _) in enumerate(statuses):
        if status == "conserved":
            last_conserved, n_conserved = position, n_conserved + 1
        elif status != "unknown":
            break
    trim = last_conserved + 1 if n_conserved >= 2 else 0
    consensus = "".join(code for _, code, _ in statuses[:trim])
    profile_text = " ".join(f"{code}:{fraction:.2f}" if status != "unknown" else "N:-"
                            for status, code, fraction in statuses)
    warning = None
    if trim == len(counts):
        warning = "every profiled position is conserved: the remnant may be longer (see --max_trim)"
    elif trim == 0 and motifs and n_reads:
        # a remnant behind variable length spacers: one 4-mer near the start of most reads, at
        # changing positions (a given 4-mer is in the first 12 nt of ~3 % of random reads)
        motif, n_with = max(((m, c) for m, c in motifs.items() if b"N" not in m), key=lambda mc: mc[1],
                            default=(b"", 0))
        if n_with > 0.5 * n_reads:
            warning = (f"{motif.decode()} is in the first {MOTIF_WINDOW} nt of {100 * n_with / n_reads:.0f} % of the"
                       f" reads, at variable positions: variable length spacers before the restriction site?"
                       f" Nothing trimmed: remove the spacers first, or give --trim_r1 / --trim_r2")
    return trim, consensus, profile_text, warning


###############################################################################
# trimming
###############################################################################

def trim_file(source, target, trim, threads):
    reader = Reader(source)
    writer = Writer(target, threads)
    out = writer.handle
    n = n_short = 0
    try:
        for header, sequence, quality in records(reader.handle):
            n += 1
            sequence = sequence[trim:]
            if not sequence:
                sequence, n_short = b"N", n_short + 1
                quality = b"!" if quality is not None else None
            elif quality is not None:
                quality = quality[trim:]
            if quality is None:
                out.write(header + b"\n" + sequence + b"\n")
            else:
                out.write(header + b"\n" + sequence + b"\n+\n" + quality + b"\n")
    finally:
        reader.close()
        writer.close()
    return n, n_short


def process_file(job):
    """One read file: detection (if needed) and trimming. Returns a report dict."""
    index, path, mate, fixed, args, out_dir = job
    report = {"file": path, "mate": mate, "mode": "manual" if fixed is not None else "auto",
              "reads_profiled": 0, "consensus": ".", "profile": ".", "warning": None}
    if fixed is None:
        counts, motifs, n = profile(path, args.n_reads, args.max_trim)
        report["reads_profiled"] = n
        trim, consensus, profile_text, warning = detect(counts, args.min_fraction, motifs, n)
        report.update(consensus=consensus or ".", profile=profile_text, warning=warning)
        if n < MIN_COUNT:
            report["warning"] = f"only {n} reads: nothing trimmed"
            trim = 0
    else:
        trim = fixed
    report["trim"] = trim
    report["output"] = path
    if trim > 0:
        stem = SEQUENCE_EXTENSIONS.sub("", os.path.basename(path))
        kind = "fastq" if first_character(path) == b"@" else "fasta"
        target = os.path.join(out_dir, f"{index}_{stem}.trimmed.{kind}.gz")
        n, n_short = trim_file(path, target, trim, args.compression_threads)
        report.update(output=target, reads=n, reads_too_short=n_short)
    return report


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
    parser.add_argument("--min_fraction", type=float, default=0.9,
                        help="frequency of the major nucleotide of a conserved position [0.9]")
    parser.add_argument("--threads", type=int, default=0, help="files processed at once [0: all cores]")
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)          # the fofs written list absolute paths
    threads = args.threads or os.cpu_count() or 1
    args.compression_threads = 2

    samples = read_fof(args.fof)
    os.makedirs(args.out_dir, exist_ok=True)
    jobs = []
    for sample_index, (_, files) in enumerate(samples, 1):
        for file_index, (path, mate) in enumerate(files, 1):
            fixed = args.trim_r1 if mate == 1 else args.trim_r2
            jobs.append((f"{sample_index}_{file_index}", path, mate, fixed, args, args.out_dir))
    with ProcessPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
        reports = list(pool.map(process_file, jobs))

    report_file = os.path.join(args.out_dir, "trimming_report.tsv")
    with open(report_file, "w") as handle:
        handle.write("file\tmate\tmode\treads_profiled\ttrimmed\tconsensus\tprofile_of_the_first_positions"
                     "\ttrimmed_file\twarning\n")
        for r in reports:
            handle.write(f"{r['file']}\tR{r['mate']}\t{r['mode']}\t{r['reads_profiled']}\t{r['trim']}\t"
                         f"{r['consensus']}\t{r['profile']}\t{r['output']}\t{r['warning'] or '.'}\n")
    for r in reports:
        text = f"[trim] {os.path.basename(r['file'])} (R{r['mate']}, {r['mode']}): {r['trim']} nt"
        if r["mode"] == "auto":
            text += f" [{r['consensus']}]" if r["trim"] else " (no conserved 5' end)"
        if r.get("reads_too_short"):
            text += f", {r['reads_too_short']} reads not longer than that replaced by N"
        log(text)
        if r["warning"]:
            log(f"[trim] WARNING {os.path.basename(r['file'])}: {r['warning']}")
    by_mate = {}
    for r in reports:
        by_mate.setdefault(r["mate"], Counter())[r["trim"]] += 1
    for mate, lengths in sorted(by_mate.items()):
        if len(lengths) > 1:
            log(f"[trim] WARNING: the R{mate} files are not trimmed to the same length: "
                + ", ".join(f"{n} files {t} nt" for t, n in sorted(lengths.items()))
                + f" (see {report_file})")

    if not any(r["trim"] for r in reports):
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
