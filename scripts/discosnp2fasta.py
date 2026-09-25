#! /usr/bin/env python3

#######################################
# A script to generate fasta, nexus and phylip files from
# DiscoRad/DiscoSnpRad generated pseudo-fasta files
# Also generates multiline fasta like file (.pseudo.fasta)
# Also generates partition list
#
# Two inputs are possible:
#   -i   the DiscoSnpRad/DiscoSnp++ fasta file: one locus per bubble, the two
#        alleles of a sample are the two paths of the bubble (only correct when
#        the locus has two haplotypes in the whole data set)
#   -H   the prefix of the disco_haplotypes.py outputs (DiscoSnpRad run with -H):
#        one locus per DiscoSnp locus, sites merged and phased with the reads,
#        any number of haplotypes and alleles per site
#
# No dependencies outside the python standard library
#
# Author:  Tomas Hrbek
# Email: hrbek@evoamazon.net
# Date: 24.09.2026
# Version: 1.0
#######################################

__author__ = 'legal'

import os
import re
import sys
import argparse
from collections import defaultdict

# IUPAC codes of two different nucleotides; indel heterozygotes become the transversion
# heterozygote of the nucleotide, as in the previous versions
IUPAC = {}
for pair, code in (('AG', 'R'), ('CT', 'Y'), ('CG', 'S'), ('AT', 'W'), ('GT', 'K'), ('AC', 'M'),
                   ('-A', 'W'), ('-G', 'S'), ('-C', 'M'), ('-T', 'K')):
    IUPAC[(pair[0], pair[1])] = code
    IUPAC[(pair[1], pair[0])] = code

GENOTYPE_PATTERN = re.compile(r'^(G[0-9]+)_([0-9.])[/|]([0-9.])')


def str2bool(value):
    # solving issues of parsing booleans: https://from-locals.com/python-argparse-bool/
    return value in ('True', 'TRUE', 'T', 'true', 't', '1', 'yes')


def ambiguity(a, b):
    """IUPAC code of two characters (identical characters are returned unchanged)."""
    if a == b:
        return a
    return IUPAC.get((a.upper(), b.upper()), 'N')


def read_lookup(lookup_file):
    """Gx -> taxon, in the order of the lookup table. Duplicated names are an error."""
    include = {}
    with open(lookup_file) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                sys.exit('ERROR: lookup table line without two tab separated columns: {0}'.format(line.strip()))
            sample, taxon = fields[0].strip(), fields[1].strip()
            if sample in include:
                sys.exit('ERROR: {0} appears twice in the lookup table'.format(sample))
            include[sample] = taxon
    taxa = list(include.values())
    duplicated = {x for x in taxa if taxa.count(x) > 1}
    if duplicated:
        sys.exit('ERROR: taxon names used for several samples in the lookup table: {0}'.format(', '.join(sorted(duplicated))))
    return include


class Alignment:
    """Per taxon lists of sequence chunks (one chunk per locus), joined only at the end.

    The chunks of samples with the same genotype are the same string object, so a
    locus costs one reference per taxon, whatever its length.
    """

    def __init__(self, taxa, cons):
        self.taxa = taxa
        self.cons = cons
        self.chunks = {x: ([], []) for x in taxa}
        self.partitions = []
        self.length = 0

    def add_locus(self, sequences):
        """sequences: {taxon: (allele_0, allele_1)} or {taxon: (consensus,)} for every taxon."""
        n = len(next(iter(sequences.values()))[0])
        for taxon in self.taxa:
            for k, chunk in enumerate(sequences[taxon]):
                self.chunks[taxon][k].append(chunk)
        self.partitions.append('CHARSET p{0}={1}-{2};'.format(len(self.partitions) + 1, self.length + 1, self.length + n))
        self.length += n

    def rows(self):
        """[(name, sequence)] in the order of the lookup table."""
        rows = []
        for taxon in self.taxa:
            if self.cons:
                rows.append((taxon, ''.join(self.chunks[taxon][0])))
            else:
                rows.append((taxon + '_0', ''.join(self.chunks[taxon][0])))
                rows.append((taxon + '_1', ''.join(self.chunks[taxon][1])))
        return rows


########
# input 1: the DiscoSnpRad/DiscoSnp++ fasta file

def read_bubbles(fasta_in):
    """Yields (header_higher, sequence_higher, header_lower, sequence_lower), streaming."""
    with open(fasta_in) as handle:
        pending = []
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            pending.append(line)
            if len(pending) == 4:
                if not (pending[0].startswith('>') and pending[2].startswith('>')):
                    sys.exit('ERROR: {0} is not a DiscoSnp fasta file (one line per sequence, paths in pairs)'.format(fasta_in))
                if pending[0].split('|')[0].replace('higher', 'lower') != pending[2].split('|')[0]:
                    sys.exit('ERROR: {0} is not followed by its lower path'.format(pending[0].split('|')[0]))
                yield pending[0][1:], pending[1], pending[2][1:], pending[3]
                pending = []
        if pending:
            sys.exit('ERROR: {0} ends with an incomplete bubble'.format(fasta_in))


def parse_bubble_header(header):
    """rank, number of polymorphisms and {Gx: (allele, allele)} from a DiscoSnp header."""
    rank, nb_pol, genotypes = None, None, {}
    for field in header.split('|'):
        if field.startswith('G') and field[1:2].isdigit():
            match = GENOTYPE_PATTERN.match(field)
            if match:
                genotypes[match.group(1)] = (match.group(2), match.group(3))
        elif field.startswith('nb_pol_'):
            nb_pol = int(field[7:])
        elif field.startswith('rank_'):
            rank = float(field[5:])
    return rank, nb_pol, genotypes


def bubbles_to_alignment(args, include, output_handle, stats):
    taxa = list(include.values())
    alignment = Alignment(taxa, args.cons)
    for header_a, seq_a, header_b, seq_b in read_bubbles(os.path.join(args.path, args.infile)):
        stats['available'] += 1
        # INDEL paths have different lengths: their positions cannot be aligned by padding
        if header_a.startswith('INDEL'):
            stats['indel'] += 1
            continue
        rank, nb_pol, genotypes = parse_bubble_header(header_a)
        # filter by rank
        if rank is not None and rank < args.minimum_rank:
            stats['rank'] += 1
            continue
        # filter on number of polymorphic sites per locus
        if nb_pol is not None and nb_pol < args.minimum_polymorphisms:
            stats['poly'] += 1
            continue
        # filter by missingness (over all the samples of the header, as before)
        nb_missing = sum(1 for g in genotypes.values() if g[0] == '.')
        if not genotypes or nb_missing / len(genotypes) >= args.maximum_missingness:
            stats['miss'] += 1
            continue
        # check if sequences of different length and fix
        if len(seq_a) < len(seq_b):
            seq_a = seq_a + seq_b[len(seq_a):]
        elif len(seq_a) > len(seq_b):
            seq_b = seq_b + seq_a[len(seq_b):]

        # the four possible sequences of a sample, computed once per locus
        missing = 'N' * len(seq_a)
        paths = {'0': seq_a, '1': seq_b}
        if args.cons:
            het = ''.join(a if a == b else ambiguity(a, b) for a, b in zip(seq_a, seq_b))
        locus = header_a.split('|')[0]
        # in allele mode a sequence is named after the path it comes from (as before)
        path_names = {'0': locus, '1': header_b.split('|')[0]}
        sequences = {}
        for sample, (g1, g2) in genotypes.items():
            taxon = include.get(sample)
            if taxon is None:
                continue
            if args.cons:
                if g1 == '.' or g2 == '.':
                    seq = missing
                elif g1 == g2:
                    seq = paths[g1]
                else:
                    seq = het
                sequences[taxon] = (seq,)
                # write fasta to a common file
                if args.loc_info:
                    output_handle.write('>{0}_c_{1}\n{2}\n'.format(taxon, locus, seq))
                else:
                    output_handle.write('>{0}_{1}\n{2}\n'.format(taxon, locus, seq))
            else:
                alleles = (paths.get(g1, missing) if g2 != '.' else missing,
                           paths.get(g2, missing) if g1 != '.' else missing)
                sequences[taxon] = alleles
                # write fasta to a common file
                for k, g, default in ((0, g1, '0'), (1, g2, '1')):
                    if args.loc_info:
                        output_handle.write('>{0}_{1}_{2}\n{3}\n'.format(taxon, k, path_names.get(g, path_names[default]), alleles[k]))
                    else:
                        output_handle.write('>{0}_{1}\n{2}\n'.format(taxon, k, alleles[k]))
        if not sequences:
            stats['no_sample'] += 1
            continue
        # taxa of the lookup table absent from this locus are missing data
        for taxon in taxa:
            if taxon not in sequences:
                sequences[taxon] = (missing,) if args.cons else (missing, missing)
        output_handle.write('\n')
        alignment.add_locus(sequences)
    return alignment


########
# input 2: the outputs of disco_haplotypes.py (<prefix>.vcf and <prefix>_loci.fa)

def read_loci_fasta(loci_fasta):
    """Yields (locus name, consensus sequence), streaming."""
    with open(loci_fasta) as handle:
        name = None
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:]
            elif name is not None:
                yield name, line
                name = None


def read_vcf_loci(vcf_file):
    """Yields (locus name, samples, [(position, [alleles], info, [sample fields])]) grouped by locus."""
    samples = None
    locus, sites = None, []
    with open(vcf_file) as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = fields[9:]
                continue
            if fields[0] != locus:
                if sites:
                    yield locus, samples, sites
                locus, sites = fields[0], []
            info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
            sites.append((int(fields[1]), [fields[3]] + fields[4].split(','), info, fields[9:]))
    if sites:
        yield locus, samples, sites


def haplotypes_to_alignment(args, include, output_handle, stats):
    taxa = list(include.values())
    alignment = Alignment(taxa, args.cons)
    prefix = os.path.join(args.path, args.haplotypes)
    consensus = read_loci_fasta(prefix + '_loci.fa')
    for locus, samples, sites in read_vcf_loci(prefix + '.vcf'):
        stats['available'] += 1
        # the loci are in the same order in both files
        for name, reference in consensus:
            if name == locus:
                break
        else:
            sys.exit('ERROR: {0} is missing from {1}_loci.fa'.format(locus, prefix))
        # filter by rank: the lowest rank of the bubbles of the locus
        ranks = [float(info['RK']) for _, _, info, _ in sites if info.get('RK', '.') != '.']
        if ranks and min(ranks) < args.minimum_rank:
            stats['rank'] += 1
            continue
        # filter on number of polymorphic sites per locus
        if len(sites) < args.minimum_polymorphisms:
            stats['poly'] += 1
            continue
        # genotype of every sample at every site: (allele index, allele index, phased)
        genotypes = []
        for column in range(len(samples)):
            genotype = []
            for _, _, _, fields in sites:
                gt = fields[column].split(':', 1)[0]
                phased = '|' in gt
                a, b = gt.replace('|', '/').split('/')
                genotype.append((a, b, phased))
            genotypes.append(tuple(genotype))
        # filter by missingness: a sample is missing when no site of the locus is called
        nb_missing = sum(1 for g in genotypes if all(a == '.' for a, _, _ in g))
        if nb_missing / len(samples) >= args.maximum_missingness:
            stats['miss'] += 1
            continue

        # sequences of every distinct multi-site genotype, computed once per locus
        cache = {}
        missing = 'N' * len(reference)
        sequences = {}
        for sample, genotype in zip(samples, genotypes):
            taxon = include.get(sample)
            if taxon is None:
                continue
            if genotype not in cache:
                hap_0, hap_1 = list(reference), list(reference)
                for (position, alleles, _, _), (a, b, phased) in zip(sites, genotype):
                    if a == '.' or b == '.':
                        hap_0[position - 1] = hap_1[position - 1] = 'N'
                        continue
                    a, b = alleles[int(a)], alleles[int(b)]
                    if phased or a == b:
                        hap_0[position - 1], hap_1[position - 1] = a, b
                    else:
                        # unphased heterozygote: ambiguous in both haplotypes
                        hap_0[position - 1] = hap_1[position - 1] = ambiguity(a, b)
                if all(a == '.' for a, _, _ in genotype):
                    hap_0, hap_1 = missing, missing
                else:
                    hap_0, hap_1 = ''.join(hap_0), ''.join(hap_1)
                if args.cons:
                    cache[genotype] = (''.join(ambiguity(x, y) for x, y in zip(hap_0, hap_1)),)
                else:
                    cache[genotype] = (hap_0, hap_1)
            sequences[taxon] = cache[genotype]
            # write fasta to a common file
            if args.cons:
                if args.loc_info:
                    output_handle.write('>{0}_c_{1}\n{2}\n'.format(taxon, locus, cache[genotype][0]))
                else:
                    output_handle.write('>{0}_{1}\n{2}\n'.format(taxon, locus, cache[genotype][0]))
            else:
                for k in (0, 1):
                    if args.loc_info:
                        output_handle.write('>{0}_{1}_{2}\n{3}\n'.format(taxon, k, locus, cache[genotype][k]))
                    else:
                        output_handle.write('>{0}_{1}\n{2}\n'.format(taxon, k, cache[genotype][k]))
        if not sequences:
            stats['no_sample'] += 1
            continue
        for taxon in taxa:
            if taxon not in sequences:
                sequences[taxon] = (missing,) if args.cons else (missing, missing)
        output_handle.write('\n')
        alignment.add_locus(sequences)
    return alignment


########
# output

def nexus_name(name):
    return name if re.match(r'^[A-Za-z0-9_.-]+$', name) else "'" + name.replace("'", "''") + "'"


def make_fasta(rows, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.fas'), 'w') as output_handle:
        for name, seq in rows:
            output_handle.write('>{0}\n{1}\n'.format(name, seq))


def make_nexus(rows, partitions, path, fasta_out):
    names = [nexus_name(name) for name, _ in rows]
    width = max(len(x) for x in names) + 5
    with open(os.path.join(path, fasta_out + '.nex'), 'w') as output_handle:
        output_handle.write('#NEXUS\n')
        output_handle.write('BEGIN DATA;\n')
        output_handle.write('  DIMENSIONS NTAX={0} NCHAR={1};\n'.format(len(rows), len(rows[0][1])))
        output_handle.write('  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n')
        output_handle.write('  MATRIX\n')
        for name, (_, seq) in zip(names, rows):
            output_handle.write('{0}{1}\n'.format(name.ljust(width), seq))
        output_handle.write(';\n')
        output_handle.write('END;\n\n')
        # add partition information
        output_handle.write('BEGIN SETS;\n')
        output_handle.write('\n'.join(partitions) + '\n')
        output_handle.write('END;\n')


def make_phylip(rows, path, fasta_out):
    width = max(len(name) for name, _ in rows) + 5
    with open(os.path.join(path, fasta_out + '.phy'), 'w') as output_handle:
        output_handle.write('{0} {1}\n'.format(len(rows), len(rows[0][1])))
        for name, seq in rows:
            output_handle.write('{0}{1}\n'.format(name.ljust(width), seq))


########
# main

def main():
    parser = argparse.ArgumentParser(description='Script to generate a fasta file from DiscoSnpRad/DiscoSnp++ pseudo-fasta file')
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('-i', '--infile', help='DiscoSnpRad/DiscoSnp++ pseudo-fasta file', type=str)
    source.add_argument('-H', '--haplotypes', help='prefix of the disco_haplotypes.py outputs (<prefix>.vcf and <prefix>_loci.fa); '
                        'DiscoSnpRad run with -H', type=str)
    parser.add_argument('-o', '--outfile', help='Output file base name (fasta and nexus)', type=str, required=True)
    parser.add_argument('-l', '--lookup', help='Gx to sample lookup table (tab separated); in project directory', type=str, required=True)
    parser.add_argument('-path', help='Path to project directory; default ./', type=str, required=False, default='./')
    parser.add_argument('-cons', help='Generate consensus sequence; default True', type=str, required=False, default='True')
    parser.add_argument('-loc_info', help='Include locus name in header; default False', type=str, required=False, default='False')
    parser.add_argument('-min_rank', '--minimum_rank', help='minimum rank; parallog metric, default .4 (decimal)', type=float, required=False, default=.4)
    parser.add_argument('-max_miss', '--maximum_missingness', help='maximum missing data per locus, default .5 (decimal)', type=float, required=False, default=.5)
    parser.add_argument('-min_poly', '--minimum_polymorphisms', help='minimum number of SNPs per locus, default 3 (integer)', type=int, required=False, default=3)
    args = parser.parse_args()
    args.cons = str2bool(args.cons)
    args.loc_info = str2bool(args.loc_info)

    include = read_lookup(os.path.join(args.path, args.lookup))
    stats = defaultdict(int)
    with open(os.path.join(args.path, args.outfile + '.pseudo.fasta'), 'w') as output_handle:
        if args.infile:
            alignment = bubbles_to_alignment(args, include, output_handle, stats)
        else:
            alignment = haplotypes_to_alignment(args, include, output_handle, stats)
    with open(os.path.join(args.path, args.outfile + '.partitions'), 'w') as partition_handle:
        partition_handle.write(''.join(x + '\n' for x in alignment.partitions))

    print('********')
    print('Number of loci available: {0}'.format(stats['available']))
    for key, text in (('indel', 'INDEL bubbles skipped'), ('rank', 'loci below the minimum rank'),
                      ('poly', 'loci with too few polymorphisms'), ('miss', 'loci with too much missing data'),
                      ('no_sample', 'loci without any sample of the lookup table')):
        if stats[key]:
            print('  {0}: {1}'.format(text, stats[key]))
    print('Number of loci extracted: {0}'.format(len(alignment.partitions)))
    if not alignment.partitions:
        print('No locus passed the filters: no fasta, nexus or phylip file written')
        return

    # write out fasta, nexus and phylip formats
    rows = alignment.rows()
    make_fasta(rows, args.path, args.outfile)
    make_nexus(rows, alignment.partitions, args.path, args.outfile)
    make_phylip(rows, args.path, args.outfile)


if __name__ == '__main__':
    main()
