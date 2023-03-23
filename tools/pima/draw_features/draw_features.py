#!/usr/bin/env python

import argparse
import os
import random

import matplotlib.pyplot as pyplot
import pandas
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord


AMR_COLOR = '#FED976'
INC_GROUPS_COLOR = '#0570B0'
FEATURE_COLORS = [AMR_COLOR, INC_GROUPS_COLOR]
FIGURE_WIDTH = 13


def get_random_color():
    number_of_colors = 16
    colors = ['#%s' % ' '.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(number_of_colors)]
    return random.choice(colors)


def draw_features(feature_hits_files, contigs, output_dir):
    ofh = open('process_log', 'w')
    # Read feature_hits_files.
    feature_hits = pandas.Series(dtype=object)
    feature_plots = pandas.Series(dtype=object)
    for feature_hits_file in feature_hits_files:
        feature_name = os.path.basename(feature_hits_file)
        # Make sure the file is not empty.
        if os.path.isfile(feature_hits_file) and os.path.getsize(feature_hits_file) > 0:
            best_hits = pandas.read_csv(filepath_or_buffer=feature_hits_file, sep='\t', header=None)
            ofh.write("\nFeature file %s will be processed\n" % os.path.basename(feature_hits_file))
        else:
            ofh.write("\nEmpty feature file %s will NOT be processed\n" % os.path.basename(feature_hits_file))
            best_hits = None
        feature_hits[feature_name] = best_hits

    # Draw one plot per contig for simplicity.
    ofh.write("\nProcessing contigs file: %s\n" % str(contigs))
    for contig in SeqIO.parse(contigs, 'fasta'):
        ofh.write("Processing contig: %s\n" % str(contig))
        contig_plot_png = os.path.join(output_dir, '%s.png' % str(contig.id))
        feature_sets_to_plot = pandas.Series(dtype=object)
        for feature_number in range(len(feature_hits)):
            feature_name = feature_hits.index.to_list()[feature_number]
            ofh.write("Processing feature name: %s\n" % str(feature_name))
            these_features = feature_hits[feature_name]
            if these_features is None or these_features.shape[0] == 0:
                # No features.
                continue
            contig_features = these_features.loc[these_features.iloc[:, 0] == contig.id, :]
            if contig_features is None or contig_features.shape[0] == 0:
                # No features.
                continue
            features_to_plot = []
            for i in range(contig_features.shape[0]):
                i = contig_features.iloc[i, :]
                if feature_number <= len(FEATURE_COLORS):
                    color = FEATURE_COLORS[feature_number]
                else:
                    color = get_random_color()
                features_to_plot += [GraphicFeature(start=i[1], end=i[2], label=i[3], strand=1 * i[5], color=color)]
            feature_sets_to_plot[feature_name] = features_to_plot
        ofh.write("Number of features to plot: %d\n" % len(feature_sets_to_plot))
        if len(feature_sets_to_plot) == 0:
            # No features.
            continue
        # Determine each plot height for later scaling
        expected_plot_heights = []
        for i in range(len(feature_sets_to_plot)):
            record = GraphicRecord(sequence_length=len(contig), features=feature_sets_to_plot[i])
            if i == len(feature_sets_to_plot) - 1:
                with_ruler = True
            else:
                with_ruler = False
            plot, _ = record.plot(figure_width=FIGURE_WIDTH, with_ruler=with_ruler)
            expected_plot_heights += [plot.figure.get_size_inches()[1]]
        plot_height_sum = sum(expected_plot_heights)
        # Make a figure with separate plots for each feature class.
        plots = pyplot.subplots(nrows=len(feature_sets_to_plot),
                                ncols=1,
                                sharex=True,
                                figsize=(FIGURE_WIDTH, plot_height_sum * .66666),
                                gridspec_kw={"height_ratios": expected_plot_heights})
        figure = plots[0]
        plots = plots[1]
        if len(feature_sets_to_plot) == 1:
            plots = [plots]
        # Add each feature class's plot with the pre-determined height.
        for i in range(len(feature_sets_to_plot)):
            record = GraphicRecord(sequence_length=len(contig), features=feature_sets_to_plot[i])
            if i == len(feature_sets_to_plot) - 1:
                with_ruler = True
            else:
                with_ruler = False
            plot, _ = record.plot(ax=plots[i], with_ruler=with_ruler, figure_width=FIGURE_WIDTH)
            ymin, ymax = plot.figure.axes[0].get_ylim()
            if i == 0:
                plot.text(x=0, y=ymax, s=contig.id)
        figure.tight_layout()
        ofh.write("Saving PNG plot file: %s\n" % str(contig_plot_png))
        figure.savefig(contig_plot_png)
        feature_plots[contig.id] = contig_plot_png
    ofh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--feature_hits_dir', action='store', dest='feature_hits_dir', help='Directory of tabular files containing feature hits')
    parser.add_argument('--contigs', action='store', dest='contigs', help='Fasta file of contigs')
    parser.add_argument('--output_dir', action='store', dest='output_dir', help='Output directory')

    args = parser.parse_args()

    # Get thge collection of feature hits files.  The collection
    # will be sorted alphabetically and will contain 2 files
    # named something like AMR_CDS_311_2022_12_20.fasta and
    # Incompatibility_Groups_2023_01_01.fasta.
    feature_hits_files = []
    for file_name in sorted(os.listdir(args.feature_hits_dir)):
        file_path = os.path.abspath(os.path.join(args.feature_hits_dir, file_name))
        feature_hits_files.append(file_path)

    draw_features(feature_hits_files, args.contigs, args.output_dir)
