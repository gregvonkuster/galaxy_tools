#!/usr/bin/env python

import argparse
import csv
import os

import numpy
import pandas
import matplotlib.pyplot as pyplot


def get_amr_in_feature_hits(amr_feature_hits):
    for k in amr_feature_hits.keys():
        if k.lower().find('amr') >= 0:
            return amr_feature_hits[k]
    return None


def draw_amr_matrix(amr_feature_hits_files, amr_deletions_file, amr_mutations_file, amr_gene_drug_file, output_dir):
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
        ofh.write("\namr_hits.shape[0]:%s\n" % str(amr_hits.shape[0]))
        if amr_hits.shape[0] > 0:
            for gene_idx, gene in amr_hits.iterrows():
                ofh.write("gene_idx:%s\n" % str(gene_idx))
                ofh.write("gene:%s\n" % str(gene))
                gene_name = gene[3]
                ofh.write("gene_name: %s\n" % str(gene_name))
                ofh.write("amr_gene_drug[0]: %s\n" % str(amr_gene_drug[0]))
                drugs = amr_gene_drug.loc[amr_gene_drug[0] == gene_name, :][1]
                ofh.write("drugs:%s\n" % str(drugs))
                for drug in drugs:
                    amr_to_draw = amr_to_draw.append(pandas.Series([gene_name, drug], name=amr_to_draw.shape[0], index=amr_to_draw.columns))
                    ofh.write("\amr_to_draw:%s\n" % str(amr_to_draw))

        if amr_mutations_file is not None:
            # TODO: So far, no samples have produced mutations, so we haven't been able
            # to produce a populated VarScan VCF file of mutations - https://github.com/appliedbinf/pima_md/blob/main/pima.py#L2923.
            # The call_amr_mutations Galaxy tool will currently produce this VarScan VCF file, but we need a sample that
            # will produce a populated file.  After we find one, we'll need to figure out how to implement this loop
            # https://github.com/appliedbinf/pima_md/blob/main/pima.py#L2925 in a Galaxy tool so that the VarScan VCF
            # file will be converted to the TSV amr_mutations_file that thsi tool expects. 
            # Roll up potentially resistance conferring mutations.
            for mutation_region, mutation_hits in amr_mutations.iteritems():
                for mutation_idx, mutation_hit in mutation_hits.iterrows():
                    mutation_name = mutation_region + ' ' + mutation_hit['REF'] + '->' + mutation_hit['ALT']
                    amr_to_draw = amr_to_draw.append(pandas.Series([mutation_name, mutation_hit['DRUG']], name=amr_to_draw.shape[0], index=amr_to_draw.columns))

        if amr_deletions_file is not None:
            # TODO: So far, no samples have produced deletions, but we do have all the pices in place
            # within the workflow to receive the amr_deletions_file here, although it is currently
            # always empty...
            # Roll up deletions that might confer resistance.
            amr_deletions = pandas.read_csv(filepath_or_buffer=amr_deletions_file, sep='\t', header=None)
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
    parser.add_argument('--amr_gene_drug_file', action='store', dest='amr_gene_drug_file', help='AMR_gene_drugs tsv file')
    parser.add_argument('--output_dir', action='store', dest='output_dir', help='Output directory')

    args = parser.parse_args()

    # Get thge collection of feature hits files.  The collection
    # will be sorted alphabetically and will contain 2 files
    # named something like AMR_CDS_311_2022_12_20.fasta and
    # Incompatibility_Groups_2023_01_01.fasta.
    amr_feature_hits_files = []
    for file_name in sorted(os.listdir(args.amr_feature_hits_dir)):
        file_path = os.path.abspath(os.path.join(args.amr_feature_hits_dir, file_name))
        amr_feature_hits_files.append(file_path)

    draw_amr_matrix(amr_feature_hits_files, args.amr_deletions_file, args.amr_mutations_file, args.amr_gene_drug_file, args.output_dir)
