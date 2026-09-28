#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
remove_site_variants.py

Removes the false variants lying inside the restriction sites from the bubbles
of kissnp2 (discoSnpRAD/run_discoSnpRad.sh runs it after the redundancy
removal, before kissreads2, so that the read counts, the genotypes, the
clustering and every VCF are computed on the corrected bubbles).

The reads keep their restriction site (trim_restriction_sites.py without
--trim): a RAD locus starts with the site remnant (e.g. CTAC, TGCAG), and
its reverse complement ends it. The site is the same in every read of the
locus, so a "variant" inside it is a sequencing error or an artefact (a
genuine mutation of the site would prevent the cut).

A SNP of a bubble is inside a site when one of the two paths of the bubble
(lower-case extensions included) starts with a site found by
trim_restriction_sites.py (sites.txt) and the SNP is within its first
nucleotides, or ends with the reverse complement of a site and the SNP is
within its last nucleotides (N and degenerate codes of the site match).
  - a bubble whose SNPs are all inside a site is removed,
  - otherwise the SNPs inside a site are removed from the bubble: the lower
    path takes the nucleotide of the higher path there, and the SNP is removed
    from the P_ list of the headers (nb_pol_ updated).
INDEL bubbles are copied unchanged.
"""

import argparse
import re
import sys

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
P_FIELD = re.compile(r"^P_\d+:")


def read_sites(path):
    sites = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            site = line.split("\t")[1].strip().upper()
            if site and site != ".":
                sites.add(site)
    return sorted(sites)


def site_regex(site):
    return re.compile("".join(f"[{IUPAC.get(code, 'ACGT')}N]" for code in site))


def records(handle):
    """(header, sequence) of a fasta file of one-line sequences."""
    header = None
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            header = line
        elif header is not None:
            yield header, line
            header = None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="kissnp2 fasta file (bubbles)")
    parser.add_argument("-s", "--sites", required=True, help="sites.txt written by trim_restriction_sites.py")
    parser.add_argument("-o", "--output", required=True, help="corrected fasta file")
    args = parser.parse_args()

    sites = read_sites(args.sites)
    starts = [(len(s), site_regex(s)) for s in sites]
    ends = [(len(s), site_regex(s.translate(COMPLEMENT)[::-1])) for s in sites]

    def inside(paths, position):
        for path in paths:
            for length, pattern in starts:
                if position < length and pattern.match(path):
                    return True
            for length, pattern in ends:
                if position >= len(path) - length and pattern.fullmatch(path[-length:]):
                    return True
        return False

    n_bubbles = n_removed_bubbles = n_edited = n_snps = n_removed_snps = 0
    with open(args.input) as handle, open(args.output, "w") as out:
        pending = None
        for header, sequence in records(handle):
            if not header.startswith(">SNP_") or not sites:
                out.write(f"{header}\n{sequence}\n")
                continue
            if "_higher_path_" in header:
                pending = (header, sequence)
                continue
            higher_header, higher = pending
            lower_header, lower = header, sequence
            pending = None
            n_bubbles += 1
            fields = higher_header.split("|")
            p_index = next(i for i, f in enumerate(fields) if P_FIELD.match(f))
            polymorphisms = fields[p_index].split(",")
            left = len(higher) - len(higher.lstrip("acgtn"))
            paths = (higher.upper(), lower.upper())
            keep, removed = [], []
            for polymorphism in polymorphisms:
                position = left + int(polymorphism.split(":")[1].split("_")[0])
                (removed if inside(paths, position) else keep).append((polymorphism, position))
            n_snps += len(polymorphisms)
            n_removed_snps += len(removed)
            if not removed:
                out.write(f"{higher_header}\n{higher}\n{lower_header}\n{lower}\n")
                continue
            if not keep:
                n_removed_bubbles += 1
                continue
            n_edited += 1
            lower = list(lower)
            for _, position in removed:
                lower[position] = higher[position]
            lower = "".join(lower)
            new_p = ",".join(f"P_{i}:{p.split(':', 1)[1]}" for i, (p, _) in enumerate(keep, 1))

            def edit(h):
                f = h.split("|")
                f[p_index] = new_p
                f = [re.sub(r"^nb_pol_\d+$", f"nb_pol_{len(keep)}", x) for x in f]
                return "|".join(f)

            out.write(f"{edit(higher_header)}\n{higher}\n{edit(lower_header)}\n{lower}\n")
    sys.stderr.write(f"[site variants] sites {', '.join(sites) or 'none'}: {n_removed_snps} of {n_snps} SNPs inside a"
                     f" site removed ({n_removed_bubbles} bubbles removed, {n_edited} bubbles edited, of {n_bubbles})\n")


if __name__ == "__main__":
    main()
