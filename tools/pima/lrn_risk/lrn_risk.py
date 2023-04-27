#!/usr/bin/env python

import argparse


BLACKLIST_HEADER = ['Blacklisted Gene', 'Reason', 'Risk Category']
VFDB_HEADER = ['Gene', 'Contig', '% Identity', '% Coverage', 'E-Value', 'Annotation', 'Comparison to Publicly Available Genomes']


def get_species_from_gtdb(f):
    # get GTDB species
    # assumes there is one genome in the GTDB-Tk output file
    with open(f, 'r') as fh:
        for line in fh:
            if line.find('user_genome') < 0:
                items = line.split('\t')
                tax = items[1].strip()
                tax = tax.split(';')[-1].strip()
                # split on GTDB species tag
                tax = tax.split('s__')[1].strip()
                if len(tax) == 0:
                    tax = '(Unknown Species)'
    return tax


def get_blast_genes(f):
    # reads genes detected via BLAST
    # BLAST header is as follows:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen
    d = {}
    with open(f, 'r') as fh:
        for line in fh:
            items = line.split('\t')
            gene = items[0]
            # contig = items[1]
            # pid = items[2]
            alen = items[3]
            # e = items[-4]
            qlen = items[-1]
            # calculate query coverage by dividing alignment length by query length
            qcov = round(float(alen) / float(qlen) * 100.0, 2)
            if gene not in d.keys():
                d[gene] = []
            d[gene].append('%s\t%s' % (line, str(qcov)))
    return d


def get_blacklist(v, b):
    # identify high-risk isolates based on blacklisted genes
    # blacklisted genes file contains two columns:
    # column 0=the gene name as it appears in the gene database
    # column 1=the reason why the gene was blacklisted, which will be reported
    # e.g., 'ANTHRAX TOXIN'
    bdict = {}
    with open(b, 'r') as fh:
        for line in fh:
            items = line.split('\t')
            gene = items[0].strip()
            val = items[1].strip()
            bdict[gene] = val
    blacklist_present = {}
    for key in v.keys():
        if key in bdict.keys():
            val = bdict[key]
            blacklist_present[key] = val
    return blacklist_present


def gene_dist(f, blast, gtdb):
    # get within-species prevalence of genes
    # for virulence factors (VFs): uses VFDB VFs detected via ABRicate's VFDB db
    # for AMR genes: uses AMR genes detected via ABRicate's ResFinder db
    # for VFs and AMR genes: genes were detected via ABRicate XXX
    # minimum nucleotide identity and coverage values >=80%
    # total of 61,161 genomes queried
    # takes VFDB or AMR gene distribution file as input (f)
    # BLAST file of VFDB or AMR genes (blast)
    # GTDB species (gtdb)
    # create dictionaries based on gene distribution
    d = {}
    annd = {}
    gtdbd = {}
    with open(f, 'r') as fh:
        for line in fh:
            items = line.split('\t')
            tax = items[0].strip()
            tax = tax.split('s__')[1].strip()
            if len(tax) == 0:
                tax = '(Unknown Species)'
            gene = items[1].strip()
            ann = items[-1].strip()
            denom = items[3].strip()
            d['%s___%s' % (tax, gene)] = line
            annd[gene] = ann
            gtdbd[tax] = denom
    # parse BLAST results
    finallines = []
    for key in blast.keys():
        blastval = blast[key]
        for bv in blastval:
            testkey = '%s___%s' % (gtdb, key)
            if testkey in d.keys() and gtdb != '(Unknown Species)':
                taxval = d[testkey]
                items = taxval.split('\t')
                tax = items[0].strip()
                tax = tax.split('s__')[1].strip()
                if len(tax) == 0:
                    tax = '(Unknown Species)'
                gene = items[1].strip()
                pres = items[2].strip()
                denom = items[3].strip()
                perc = items[4].strip()
                perc = str(round(float(perc), 2))
                ann = items[-1].strip()
                freetext = 'Gene {0} has been detected in {1}% of {2} genomes ({3} of {4} genomes queried)'.format(gene, perc, tax, pres, denom)
            elif gtdb != '(Unknown Species)':
                ann = annd[key]
                denom = gtdbd[gtdb]
                freetext = 'WARNING: Gene {0} ({1}) has never been detected in species {2} (n={3} genomes queried)! Interpret with caution!'.format(key, ann, gtdb, denom)
            else:
                ann = annd[key]
                freetext = 'WARNING: Genome belongs to an undescribed species. Interpret with caution!'
            finalline = '%s\t%s\t%s' % (bv, ann, freetext)
            finallines.append(finalline)
    return finallines


def output_blacklist(blacklist, blacklist_output_file):
    # takes detected blacklisted genes as input (blacklist)
    # blacklist results
    with open(blacklist_output_file, 'w') as fh:
        fh.write('%s\n' % '\t'.join(BLACKLIST_HEADER))
        if len(blacklist.keys()) == 0:
            # print this if no blacklisted genes are detected
            fh.write('(No blacklisted genes detected)\tNA\tNot high risk\n')
        else:
            # print this if blacklisted genes are detected
            # print a table with one row per detected blacklisted gene
            for key in blacklist.keys():
                val = blacklist[key]
                fh.write('%s\t%s\tHIGH RISK\n' % (key, val))


def output_vfdb(vfdist, vfdb_output_file):
    # takes distribution of virulence factors as input (vfdist)
    # VFDB results
    with open(vfdb_output_file, 'w') as fh:
        fh.write('%s\n' % '\t'.join(VFDB_HEADER))
        if len(vfdist) == 0:
            # print this if no VFs detected
            fh.write('%s\n' % '\t'.join(['(No VFs Detected)'] * 7))
        else:
            # print table of VFs if VFs detected
            for vline in vfdist:
                # blast_header=['Gene', 'Contig', 'Percent (%) Nucleotide Identity', 'Alignment Length', 'Mismatches', 'Gaps', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-Value', 'Bit Score',  'Identical Matches', 'Query Length']
                # lc_header=['Query Coverage', 'Annotation', 'Comparison to Publicly Available Genomes']
                items = vline.split('\t')
                vgene = items[0].strip()
                vcontig = items[1].strip()
                vid = items[2].strip()
                vcov = items[-3].strip()
                veval = items[-7].strip()
                vann = items[-2].strip()
                vnotes = items[-1].strip()
                vfinal = [vgene, vcontig, vid, vcov, veval, vann, vnotes]
                vfinal = '\t'.join(vfinal).strip()
                fh.write('%s\n' % vfinal)


def output_amr(amrdist, amr_output_file):
    # takes distribution of AMR genes as input (amrdist)
    # AMR results
    with open(amr_output_file, 'w') as fh:
        fh.write('%s\n' % '\t'.join(VFDB_HEADER))
        if len(amrdist) == 0:
            # print this if no AMR genes detected
            fh.write('%s\n' % '\t'.join(['(No AMR Genes Detected)'] * 7))
        else:
            # print this if AMR genes detected
            for aline in amrdist:
                # blast_header=['Gene', 'Contig', 'Percent (%) Nucleotide Identity', 'Alignment Length', 'Mismatches', 'Gaps', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-Value', 'Bit Score',  'Identical Matches', 'Query Length']
                # lc_header=['Query Coverage', 'Annotation', 'Comparison to Publicly Available Genomes']
                items = aline.split('\t')
                agene = items[0].strip()
                acontig = items[1].strip()
                aid = items[2].strip()
                acov = items[-3].strip()
                aeval = items[-7].strip()
                aann = items[-2].strip()
                anotes = items[-1].strip()
                afinal = [agene, acontig, aid, acov, aeval, aann, anotes]
                afinal = '\t'.join(afinal).strip()
                fh.write('%s\n' % afinal)


# lrnrisk_prototype arguments
parser = argparse.ArgumentParser()

parser.add_argument('--gtdb_file', action='store', dest='gtdb_file', help='Path to gtdbtk tsv file')
parser.add_argument('--virulence_factors_file', action='store', dest='virulence_factors_file', help='Path to tsv virulence factors file')
parser.add_argument('--amr_determinants_file', action='store', dest='amr_determinants_file', help='Path to AMR determinants tsv file')
parser.add_argument('--blacklist_file', action='store', dest='blacklist_file', help='Path to blacklisted high-risk virulence factors tsv file')
parser.add_argument('--vf_distribution_file', action='store', dest='vf_distribution_file', help='Path to virulence factor distribution tsv file')
parser.add_argument('--amr_distribution_file', action='store', dest='amr_distribution_file', help='Path to AMR determinant distribution tsv file')
parser.add_argument('--blacklist_output_file', action='store', dest='blacklist_output_file', help='Path to blacklist output file')
parser.add_argument('--vfdb_output_file', action='store', dest='vfdb_output_file', help='Path to vfdb output file')
parser.add_argument('--amr_output_file', action='store', dest='amr_output_file', help='Path to amr output file')

# parse arguments and run pipeline
args = parser.parse_args()

# print_output(blacklist, vf_distribution, amr_distribution, args.output, species)
virulence_genes = get_blast_genes(args.virulence_factors_file)
species = get_species_from_gtdb(args.gtdb_file)

blacklist = get_blacklist(virulence_genes, args.blacklist_file)
output_blacklist(blacklist, args.blacklist_output_file)

vf_distribution = gene_dist(args.vf_distribution_file, virulence_genes, species)
output_vfdb(vf_distribution, args.vfdb_output_file)

amr_genes = get_blast_genes(args.amr_determinants_file)
amr_distribution = gene_dist(args.amr_distribution_file, amr_genes, species)
output_amr(amr_distribution, args.amr_output_file)
