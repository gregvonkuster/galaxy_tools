#!/usr/bin/env python

import argparse
import csv
import os
import subprocess
import sys
import tempfile

import Bio.SeqIO
import numpy
import pandas
import matplotlib.pyplot as pyplot


def get_amr_in_feature_hits(amr_feature_hits):
    for k in amr_feature_hits.keys():
        if k.lower().find('amr') >= 0:
            return amr_feature_hits[k]
    return None


def load_fasta(fasta_file):
    sequence = pandas.Series(dtype=object)
    for contig in Bio.SeqIO.parse(fasta_file, 'fasta'):
        sequence[contig.id] = contig
    return sequence


def run_command(cmd):
    try:
        tmp_name = tempfile.NamedTemporaryFile(dir=".").name
        tmp_stderr = open(tmp_name, 'wb')
        proc = subprocess.Popen(args=cmd, shell=True, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
        tmp_stderr.close()
        if returncode != 0:
            # Get stderr, allowing for case where it's very large.
            tmp_stderr = open(tmp_name, 'rb')
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read(buffsize)
                    if not stderr or len(stderr) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            os.remove(tmp_name)
            stop_err(stderr)
    except Exception as e:
        stop_err('Command:\n%s\n\nended with error:\n%s\n\n' % (cmd, str(e)))


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def draw_amr_matrix(amr_feature_hits_files, amr_deletions_file, amr_mutations_file, amr_mutation_regions_file, amr_gene_drug_file, reference, reference_size, region_mutations_output_file, mutations_dir, output_dir):
    ofh = open('process_log', 'w')

    # Read amr_feature_hits_files.
    amr_feature_hits = pandas.Series(dtype=object)
    for amr_feature_hits_file in amr_feature_hits_files:
        feature_name = os.path.basename(amr_feature_hits_file)
        # Make sure the file is not empty.
        if os.path.isfile(amr_feature_hits_file) and os.path.getsize(amr_feature_hits_file) > 0:
            best_hits = pandas.read_csv(filepath_or_buffer=amr_feature_hits_file, sep='\t', header=None)
            ofh.write("\nFeature file %s will be processed\n" % os.path.basename(amr_feature_hits_file))
        else:
            ofh.write("\nEmpty feature file %s will NOT be processed\n" % os.path.basename(amr_feature_hits_file))
            best_hits = None
        amr_feature_hits[feature_name] = best_hits

    amr_hits = get_amr_in_feature_hits(amr_feature_hits)
    ofh.write("\namr_hits:\n%s\n" % str(amr_hits))
    if amr_hits is not None:
        amr_to_draw = pandas.DataFrame(columns=['gene', 'drug'])
        ofh.write("\namr_to_draw:\n%s\n" % str(amr_to_draw))
        # Read amr_drug_gene_file.
        amr_gene_drug = pandas.read_csv(amr_gene_drug_file, index_col=None, sep='\t', quoting=csv.QUOTE_NONE, header=None)
        ofh.write("\namr_gene_drug:\n%s\n" % str(amr_gene_drug))

        # Roll up AMR gene hits.
        ofh.write("\namr_hits.shape[0]: %s\n" % str(amr_hits.shape[0]))
        if amr_hits.shape[0] > 0:
            for gene_idx, gene in amr_hits.iterrows():
                ofh.write("gene_idx: %s\n" % str(gene_idx))
                ofh.write("gene: %s\n" % str(gene))
                gene_name = gene[3]
                ofh.write("gene_name: %s\n" % str(gene_name))
                ofh.write("amr_gene_drug[0]: %s\n" % str(amr_gene_drug[0]))
                drugs = amr_gene_drug.loc[amr_gene_drug[0] == gene_name, :][1]
                ofh.write("drugs: %s\n" % str(drugs))
                for drug in drugs:
                    amr_to_draw = amr_to_draw.append(pandas.Series([gene_name, drug], name=amr_to_draw.shape[0], index=amr_to_draw.columns))
                    ofh.write("\amr_to_draw: %s\n" % str(amr_to_draw))

        ofh.write("\namr_mutations_file si None: %s\n" % str(amr_mutations_file == 'None'))
        if amr_mutations_file not in [None, 'None'] and os.path.getsize(amr_mutations_file) > 0:
            amr_mutations = pandas.Series(dtype=object)
            if amr_mutation_regions_file is not None:
                mutation_regions = pandas.read_csv(amr_mutation_regions_file, header=0, sep='\t', index_col=False)
                if mutation_regions.shape[1] != 7:
                    ofh.write("\nMutation regions should be a six column file.\n")
                elif mutation_regions.shape[0] == 0:
                    ofh.write("\nNo rows in mutation regions file.\n")
                else:
                    # Make sure the positions in the BED file fall within the chromosomes provided in the reference sequence.
                    for mutation_region in range(mutation_regions.shape[0]):
                        mutation_region = mutation_regions.iloc[mutation_region, :]
                        if not (mutation_region[0] in reference):
                            ofh.write("\nMutation region :%s not found in reference genome.\n" % ' '.join(mutation_region.astype(str)))
                            continue
                        if not isinstance(mutation_region[1], numpy.int64):
                            ofh.write("\nNon-integer found in mutation region start (column 2): %s.\n" % str(mutation_region[1]))
                            break
                        elif not isinstance(mutation_region[2], numpy.int64):
                            ofh.write("\nNon-integer found in mutation region start (column 3): %s.\n" % str(mutation_region[2]))
                            break
                        if mutation_region[1] <= 0 or mutation_region[2] <= 0:
                            ofh.write("\nMutation region %s starts before the reference sequence.\n" % ' '.join(mutation_region.astype(str)))
                        if mutation_region[1] > len(reference[mutation_region[0]].seq) or mutation_region[2] > len(reference[mutation_region[0]].seq):
                            ofh.write("\nMutation region %s ends after the reference sequence.\n" % ' '.join(mutation_region.astype(str)))
                    for region_i in range(mutation_regions.shape[0]):
                        region = mutation_regions.iloc[region_i, :]
                        if not region.get('type', default='No Type') in ['snp', 'small-indel', 'any']:
                            continue
                        ofh.write("\nFinding AMR mutations for %s.\n" % str(region['name']))
                        region_dir = os.path.join(mutations_dir, 'region_' + str(region_i))
                        os.mkdir(region_dir)
                        region_bed = os.path.join(region_dir, 'region.bed')
                        mutation_regions.loc[[region_i], ].to_csv(path_or_buf=region_bed, sep='\t', header=False, index=False)
                        cmd = "bedtools intersect -nonamecheck -wb -a '%s' -b '%s' | awk '{BEGIN{getline < \"%s\" ;printf $0\"\t\";getline < \"%s\"; getline < \"%s\";print $0}{print}' > %s" % (region_bed, amr_mutations_file, amr_mutation_regions_file, amr_mutations_file, amr_mutations_file, region_mutations_output_file)
                        run_command(cmd)
                        try:
                            region_mutations = pandas.read_csv(region_mutations_output_file, sep='\t', header=0, index_col=False)
                        except Exception:
                            region_mutations = pandas.DataFrame()
                        if region_mutations.shape[0] == 0:
                            continue
                        # Figure out what kind of mutations are in this region.
                        region_mutation_types = pandas.Series(['snp'] * region_mutations.shape[0], name='TYPE', index=region_mutations.index)
                        region_mutation_types[region_mutations['REF'].str.len() != region_mutations['ALT'].str.len()] = 'small-indel'
                        region_mutation_drugs = pandas.Series(region['drug'] * region_mutations.shape[0], name='DRUG', index=region_mutations.index)
                        region_notes = pandas.Series(region['note'] * region_mutations.shape[0], name='NOTE', index=region_mutations.index)
                        region_mutations = pandas.concat([region_mutations, region_mutation_types, region_mutation_drugs, region_notes], axis=1)
                        region_mutations = region_mutations[['#CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'DRUG', 'NOTE']]
                        amr_mutations[region['name']] = region_mutations
            else:
                ofh.write("\nMutation region BED not received.\n")
            # Roll up potentially resistance conferring mutations.
            for mutation_region, mutation_hits in amr_mutations.iteritems():
                for mutation_idx, mutation_hit in mutation_hits.iterrows():
                    mutation_name = mutation_region + ' ' + mutation_hit['REF'] + '->' + mutation_hit['ALT']
                    amr_to_draw = amr_to_draw.append(pandas.Series([mutation_name, mutation_hit['DRUG']], name=amr_to_draw.shape[0], index=amr_to_draw.columns))

        if amr_deletions_file not in [None, 'None'] and os.path.getsize(amr_deletions_file) > 0:
            # TODO: So far, no samples have produced deletions, but we do have all the pices in place
            # within the workflow to receive the amr_deletions_file here, although it is currently
            # always empty...
            # Roll up deletions that might confer resistance.
            try:
                amr_deletions = pandas.read_csv(filepath_or_buffer=amr_deletions_file, sep='\t', header=None)
            except Exception:
                amr_deletions = pandas.DataFrame()
            if amr_deletions.shape[0] > 0:
                amr_deletions.columns = ['contig', 'start', 'stop', 'name', 'type', 'drug', 'note']
                amr_deletions = amr_deletions.loc[amr_deletions['type'].isin(['large-deletion', 'any']), :]
                for deletion_idx, deleted_gene in amr_deletions.iterrows():
                    amr_to_draw = amr_to_draw.append(pandas.Series(['\u0394' + deleted_gene[3], deleted_gene[5]], name=amr_to_draw.shape[0], index=amr_to_draw.columns))

        if amr_to_draw.shape[0] > 1:
            ofh.write("\nDrawing AMR matrix...\n")
            present_genes = amr_to_draw['gene'].unique()
            present_drugs = amr_to_draw['drug'].unique()
            amr_matrix = pandas.DataFrame(0, index=present_genes, columns=present_drugs)
            for hit_idx, hit in amr_to_draw.iterrows():
                amr_matrix.loc[hit[0], hit[1]] = 1
            amr_matrix_png = os.path.join(output_dir, 'amr_matrix.png')
            int_matrix = amr_matrix[amr_matrix.columns].astype(int)
            figure, axis = pyplot.subplots()
            axis.invert_yaxis()
            axis.set_yticks(numpy.arange(0.5, len(amr_matrix.index)), minor=False)
            axis.set_yticklabels(int_matrix.index.values)
            axis.set_xticks(numpy.arange(0.5, len(amr_matrix.columns)), minor=False)
            axis.set_xticklabels(amr_matrix.columns.values, rotation=90)
            axis.xaxis.tick_top()
            axis.xaxis.set_label_position('top')
            pyplot.tight_layout()
            pyplot.savefig(amr_matrix_png, dpi=300)
        else:
            ofh.write("\nEmpty AMR matrix, nothing to draw...\n")
    ofh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--amr_feature_hits_dir', action='store', dest='amr_feature_hits_dir', help='Directory of tabular files containing feature hits')
    parser.add_argument('--amr_deletions_file', action='store', dest='amr_deletions_file', default=None, help='AMR deletions BED file')
    parser.add_argument('--amr_mutations_file', action='store', dest='amr_mutations_file', default=None, help='AMR mutations TSV file')
    parser.add_argument('--amr_mutation_regions_file', action='store', dest='amr_mutation_regions_file', default=None, help='AMR mutation regions BED file')
    parser.add_argument('--amr_gene_drug_file', action='store', dest='amr_gene_drug_file', help='AMR_gene_drugs tsv file')
    parser.add_argument('--reference_genome', action='store', dest='reference_genome', help='Reference genome fasta file')
    parser.add_argument('--region_mutations_output_file', action='store', dest='region_mutations_output_file', default=None, help='Region mutations TSV output file')
    parser.add_argument('--mutations_dir', action='store', dest='mutations_dir', help='Mutations directory')
    parser.add_argument('--output_dir', action='store', dest='output_dir', help='Output directory')

    args = parser.parse_args()

    # Get the collection of feature hits files.  The collection
    # will be sorted alphabetically and will contain 2 files
    # named something like AMR_CDS_311_2022_12_20.fasta and
    # Incompatibility_Groups_2023_01_01.fasta.
    amr_feature_hits_files = []
    for file_name in sorted(os.listdir(args.amr_feature_hits_dir)):
        file_path = os.path.abspath(os.path.join(args.amr_feature_hits_dir, file_name))
        amr_feature_hits_files.append(file_path)

    # Load the reference genome into memory.
    reference = load_fasta(args.reference_genome)
    reference_size = 0
    for i in reference:
        reference_size += len(i.seq)

    draw_amr_matrix(amr_feature_hits_files, args.amr_deletions_file, args.amr_mutations_file, args.amr_mutation_regions_file, args.amr_gene_drug_file, reference, reference_size, args.region_mutations_output_file, args.mutations_dir, args.output_dir)
