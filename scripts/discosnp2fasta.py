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
#        any number of haplotypes and alleles per site; the loci are filtered by
#        default on their repeat score (RPT <= .5) and LONG flag (see info_fields.txt)
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


def data_samples(args):
    """The samples (Gx) present in the input: the VCF header (-H) or the genotype fields of the first bubble (-i)."""
    if args.infile:
        for header_a, _, _, _ in read_bubbles(os.path.join(args.path, args.infile)):
            return [field.split('_', 1)[0] for field in header_a.split('|') if GENOTYPE_PATTERN.match(field)]
        return []
    with open(os.path.join(args.path, args.haplotypes) + '.vcf') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
            if not line.startswith('#'):
                break
    sys.exit('ERROR: no #CHROM header line in {0}.vcf'.format(os.path.join(args.path, args.haplotypes)))


def match_lookup(include, samples):
    """Keeps the lookup entries of the samples present in the data; warns about the others.

    include: {Gx: taxon} of the lookup table, or None (no -l): every sample of the data, under its own name.
    When the samples of the data already carry the taxon names (VCF renamed, e.g. with bcftools reheader),
    they are not renamed: the lookup table only selects and orders them.
    """
    if include is None:
        return {sample: sample for sample in samples}
    taxa = set(include.values())
    by_id = sum(1 for sample in samples if sample in include)
    by_name = sum(1 for sample in samples if sample in taxa)
    if by_name > by_id:
        sys.stderr.write('NOTE: the samples of the data are already named ({0} of {1} match taxon names of the lookup '
                         'table): they are not renamed\n'.format(by_name, len(samples)))
        include = {taxon: taxon for taxon in include.values()}
    absent = [sample for sample in include if sample not in samples]
    if absent:
        sys.stderr.write('WARNING: {0} sample(s) of the lookup table are not in the data and are left out of the alignment: '
                         '{1}\n'.format(len(absent), ', '.join(x if x == include[x] else '{0} ({1})'.format(x, include[x])
                                                                for x in absent)))
    unlisted = [sample for sample in samples if sample not in include]
    if unlisted:
        sys.stderr.write('WARNING: {0} sample(s) of the data are not in the lookup table and are left out of the alignment: '
                         '{1}\n'.format(len(unlisted), ', '.join(unlisted)))
    kept = {sample: taxon for sample, taxon in include.items() if sample in samples}
    if not kept:
        sys.exit('ERROR: no sample of the lookup table is in the data (data samples: {0})'.format(', '.join(samples[:10])))
    return kept


class Alignment:
    """Per taxon lists of sequence chunks (one chunk per locus), joined only at the end: the consensus
    (haplotypes collapsed, IUPAC codes at the heterozygous sites) and the two haplotypes.

    The chunks of samples with the same genotype are the same string object, so a
    locus costs three references per taxon, whatever its length.
    """

    def __init__(self, taxa):
        self.taxa = taxa
        self.chunks = {x: ([], [], []) for x in taxa}
        self.partitions = []
        self.length = 0

    def add_locus(self, sequences):
        """sequences: {taxon: (consensus, haplotype_1, haplotype_2)} for every taxon."""
        n = len(next(iter(sequences.values()))[0])
        for taxon in self.taxa:
            for k, chunk in enumerate(sequences[taxon]):
                self.chunks[taxon][k].append(chunk)
        self.partitions.append('CHARSET p{0}={1}-{2};'.format(len(self.partitions) + 1, self.length + 1, self.length + n))
        self.length += n

    def rows(self, kind):
        """[(name, sequence)] in the order of the lookup table: 'consensus' (one row per taxon) or
        'haplotypes' (two rows per taxon, <taxon>_1 and <taxon>_2)."""
        rows = []
        for taxon in self.taxa:
            if kind == 'consensus':
                rows.append((taxon, ''.join(self.chunks[taxon][0])))
            else:
                rows.append((taxon + '_1', ''.join(self.chunks[taxon][1])))
                rows.append((taxon + '_2', ''.join(self.chunks[taxon][2])))
        return rows


def write_pseudo(handles, loc_info, taxon, locus, consensus, haplotypes, haplotype_names):
    """One locus of a taxon in the .pseudo.fasta files: the consensus, and the two haplotypes
    (haplotype_names: the name of the locus or bubble path of each haplotype, for -loc_info)."""
    if handles.get('consensus'):
        if loc_info:
            handles['consensus'].write('>{0}_c_{1}\n{2}\n'.format(taxon, locus, consensus))
        else:
            handles['consensus'].write('>{0}_{1}\n{2}\n'.format(taxon, locus, consensus))
    if handles.get('haplotypes'):
        for k in (0, 1):
            if loc_info:
                handles['haplotypes'].write('>{0}_{1}_{2}\n{3}\n'.format(taxon, k + 1, haplotype_names[k], haplotypes[k]))
            else:
                handles['haplotypes'].write('>{0}_{1}\n{2}\n'.format(taxon, k + 1, haplotypes[k]))


def end_pseudo_locus(handles):
    for handle in handles.values():
        if handle:
            handle.write('\n')


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


def bubbles_to_alignment(args, include, handles, stats):
    taxa = list(include.values())
    alignment = Alignment(taxa)
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
        het = ''.join(a if a == b else ambiguity(a, b) for a, b in zip(seq_a, seq_b))
        locus = header_a.split('|')[0]
        # a haplotype is named after the path it comes from (-loc_info)
        path_names = {'0': locus, '1': header_b.split('|')[0]}
        sequences = {}
        for sample, (g1, g2) in genotypes.items():
            taxon = include.get(sample)
            if taxon is None:
                continue
            if g1 == '.' or g2 == '.':
                consensus, haplotypes = missing, (missing, missing)
            else:
                consensus = paths[g1] if g1 == g2 else het
                haplotypes = (paths.get(g1, missing), paths.get(g2, missing))
            sequences[taxon] = (consensus,) + haplotypes
            write_pseudo(handles, args.loc_info, taxon, locus, consensus, haplotypes,
                         (path_names.get(g1, path_names['0']), path_names.get(g2, path_names['1'])))
        if not sequences:
            stats['no_sample'] += 1
            continue
        # taxa of the lookup table absent from this locus are missing data
        for taxon in taxa:
            if taxon not in sequences:
                sequences[taxon] = (missing, missing, missing)
        end_pseudo_locus(handles)
        alignment.add_locus(sequences)
    return alignment


########
# input 2: the outputs of disco_haplotypes.py (<prefix>.vcf and <prefix>_loci.fa)

def read_loci_fasta(loci_fasta):
    """Yields (locus name, consensus sequence), streaming. The header is
    '>locus_N|length_...|n_sites_...|positions_...': the name is its first field."""
    with open(loci_fasta) as handle:
        name = None
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:].split('|')[0]
            elif name is not None:
                yield name, line
                name = None


def read_vcf_loci(vcf_file):
    """Yields (locus name, samples, [(position, [alleles], info, [sample fields], {INFO flags})]) grouped by locus."""
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
            flags = {x for x in fields[7].split(';') if x and '=' not in x}
            sites.append((int(fields[1]), [fields[3]] + fields[4].split(','), info, fields[9:], flags))
    if sites:
        yield locus, samples, sites


def haplotypes_to_alignment(args, include, handles, stats):
    taxa = list(include.values())
    alignment = Alignment(taxa)
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
        # filter by rank: the median rank of the sites of the locus (Rk: best rank of the bubbles of a site);
        # the lowest one would let a single low ranked site discard a locus of many sites
        ranks = sorted(float(info['Rk']) for _, _, info, _, _ in sites if info.get('Rk', '.') != '.')
        median_rank = None
        if ranks:
            half = len(ranks) // 2
            median_rank = ranks[half] if len(ranks) % 2 else (ranks[half - 1] + ranks[half]) / 2
        if median_rank is not None and median_rank < args.minimum_rank:
            stats['rank'] += 1
            continue
        # filter on the repeat score: the highest RPT of the sites of the locus (RPT = 1 for the repeat
        # bubbles and the bubbles of giant components, GC)
        scores = [float(info['RPT']) for _, _, info, _, _ in sites if 'RPT' in info]
        if not scores:
            stats['no_rpt'] += 1
        elif max(scores) > args.maximum_repeat_score:
            stats['rpt'] += 1
            continue
        # filter the loci longer than one ddRAD fragment (LONG flag, run_discoSnpRad.sh -H): chimeras
        if args.skip_long and any('LONG' in flags for _, _, _, _, flags in sites):
            stats['long'] += 1
            continue
        # filter on number of polymorphic sites per locus
        if len(sites) < args.minimum_polymorphisms:
            stats['poly'] += 1
            continue
        # genotype of every sample at every site: (allele index, allele index, phased)
        genotypes = []
        for column in range(len(samples)):
            genotype = []
            for _, _, _, fields, _ in sites:
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
                for (position, alleles, _, _, _), (a, b, phased) in zip(sites, genotype):
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
                cache[genotype] = (''.join(ambiguity(x, y) for x, y in zip(hap_0, hap_1)), hap_0, hap_1)
            sequences[taxon] = cache[genotype]
            write_pseudo(handles, args.loc_info, taxon, locus, cache[genotype][0], cache[genotype][1:], (locus, locus))
        if not sequences:
            stats['no_sample'] += 1
            continue
        for taxon in taxa:
            if taxon not in sequences:
                sequences[taxon] = (missing, missing, missing)
        end_pseudo_locus(handles)
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
    parser.add_argument('-l', '--lookup', help='Gx to sample lookup table (tab separated); in project directory. Optional: '
                        'without it the samples keep their names; with -H, a VCF whose samples already carry the names '
                        'of the table (bcftools reheader) is not renamed, the table only selects and orders the samples',
                        type=str, required=False, default=None)
    parser.add_argument('-path', help='Path to project directory; default ./', type=str, required=False, default='./')
    parser.add_argument('-cons', help='both (default): write the consensus files (<out>.fas/.nex/.phy/.pseudo.fasta, haplotypes '
                        'collapsed, IUPAC codes at the heterozygous sites) and the haplotype files (<out>_haplotypes.*, '
                        'two rows per sample: <sample>_1 and <sample>_2); True: consensus files only; False: haplotype files only',
                        type=str, required=False, default='both')
    parser.add_argument('-loc_info', help='Include locus name in header; default False', type=str, required=False, default='False')
    parser.add_argument('-min_rank', '--minimum_rank', help='minimum rank (Rk with -H: the median of the sites of a locus); '
                        'parallog metric, default .4 (decimal)', type=float, required=False, default=.4)
    parser.add_argument('-max_rpt', '--maximum_repeat_score', help='-H only: skip the loci with a site of repeat score RPT above this '
                        '(repeats, giant components), default .5 (decimal); 1: no filter', type=float, required=False, default=.5)
    parser.add_argument('-skip_long', help='-H only: skip the loci flagged LONG (longer than one ddRAD fragment: chimeras); '
                        'default True', type=str, required=False, default='True')
    parser.add_argument('-max_miss', '--maximum_missingness', help='maximum missing data per locus, default .5 (decimal)', type=float, required=False, default=.5)
    parser.add_argument('-min_poly', '--minimum_polymorphisms', help='minimum number of SNPs per locus, default 3 (integer)', type=int, required=False, default=3)
    args = parser.parse_args()
    if args.cons.lower() == 'both':
        kinds = ('consensus', 'haplotypes')
    else:
        kinds = ('consensus',) if str2bool(args.cons) else ('haplotypes',)
    file_names = {'consensus': args.outfile, 'haplotypes': args.outfile + '_haplotypes'}
    args.loc_info = str2bool(args.loc_info)
    args.skip_long = str2bool(args.skip_long)

    lookup = read_lookup(os.path.join(args.path, args.lookup)) if args.lookup else None
    include = match_lookup(lookup, data_samples(args))
    stats = defaultdict(int)
    handles = {kind: open(os.path.join(args.path, file_names[kind] + '.pseudo.fasta'), 'w') for kind in kinds}
    try:
        if args.infile:
            alignment = bubbles_to_alignment(args, include, handles, stats)
        else:
            alignment = haplotypes_to_alignment(args, include, handles, stats)
    finally:
        for handle in handles.values():
            handle.close()
    with open(os.path.join(args.path, args.outfile + '.partitions'), 'w') as partition_handle:
        partition_handle.write(''.join(x + '\n' for x in alignment.partitions))

    print('********')
    print('Number of loci available: {0}'.format(stats['available']))
    for key, text in (('indel', 'INDEL bubbles skipped'), ('rank', 'loci below the minimum rank'),
                      ('rpt', 'loci above the maximum repeat score (RPT)'), ('long', 'loci flagged LONG'),
                      ('no_rpt', 'loci without RPT (older disco_haplotypes version): not filtered on it'),
                      ('poly', 'loci with too few polymorphisms'), ('miss', 'loci with too much missing data'),
                      ('no_sample', 'loci without any sample of the lookup table')):
        if stats[key]:
            print('  {0}: {1}'.format(text, stats[key]))
    print('Number of loci extracted: {0}'.format(len(alignment.partitions)))
    if not alignment.partitions:
        print('No locus passed the filters: no fasta, nexus or phylip file written')
        return

    # write out fasta, nexus and phylip formats: the consensus and/or the haplotypes
    for kind in kinds:
        rows = alignment.rows(kind)
        make_fasta(rows, args.path, file_names[kind])
        make_nexus(rows, alignment.partitions, args.path, file_names[kind])
        make_phylip(rows, args.path, file_names[kind])
        print('{0} files: {1}.fas, .nex, .phy, .pseudo.fasta ({2} rows)'.format(kind.capitalize(), file_names[kind], len(rows)))


if __name__ == '__main__':
    main()
