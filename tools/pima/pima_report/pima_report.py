#!/usr/bin/env python

import argparse
import os
import pandas
import pypandoc
import re
import subprocess
import sys

from Bio import SeqIO
from datetime import date
from mdutils.mdutils import MdUtils
# FIXME: TableOfContents doesn't work.
# from mdutils.tools import TableOfContents

CDC_ADVISORY = 'The analysis and report presented here should be treated as preliminary.  Please contact the CDC/BDRD with any results regarding _Bacillus anthracis_.'


class PimaReport:

    def __init__(self, analysis_name=None, amr_deletions_file=None, amr_matrix_files=None, assembler_version=None,
                 assembly_fasta_file=None, assembly_name=None, bedtools_version=None, blastn_version=None,
                 circos_files=None, compute_sequence_length_file=None, contig_coverage_file=None, dbkey=None,
                 dnadiff_snps_file=None, dnadiff_version=None, errors_file=None, feature_bed_files=None,
                 feature_png_files=None, flye_assembly_info_file=None, genome_insertions_file=None, gzipped=None,
                 illumina_forward_read_file=None, illumina_reverse_read_file=None, kraken2_report_file=None,
                 kraken2_version=None, lrn_risk_amr_file=None, lrn_risk_blacklist_file=None, lrn_risk_vf_file=None,
                 minimap2_version=None, mutation_regions_bed_file=None, mutation_regions_tsv_files=None,
                 ont_file=None, pima_css=None, plasmids_file=None, quast_report_file=None, read_type=None,
                 reference_insertions_file=None, samtools_version=None, varscan_version=None):
        self.ofh = open("process_log.txt", "w")

        self.ofh.write("amr_deletions_file: %s\n" % str(amr_deletions_file))
        self.ofh.write("amr_matrix_files: %s\n" % str(amr_matrix_files))
        self.ofh.write("analysis_name: %s\n" % str(analysis_name))
        self.ofh.write("assembler_version: %s\n" % str(assembler_version))
        self.ofh.write("assembly_fasta_file: %s\n" % str(assembly_fasta_file))
        self.ofh.write("assembly_name: %s\n" % str(assembly_name))
        self.ofh.write("bedtools_version: %s\n" % str(bedtools_version))
        self.ofh.write("blastn_version: %s\n" % str(blastn_version))
        self.ofh.write("circos_files: %s\n" % str(circos_files))
        self.ofh.write("compute_sequence_length_file: %s\n" % str(compute_sequence_length_file))
        self.ofh.write("contig_coverage_file: %s\n" % str(contig_coverage_file))
        self.ofh.write("dbkey: %s\n" % str(dbkey))
        self.ofh.write("dnadiff_snps_file: %s\n" % str(dnadiff_snps_file))
        self.ofh.write("dnadiff_version: %s\n" % str(dnadiff_version))
        self.ofh.write("errors_file: %s\n" % str(errors_file))
        self.ofh.write("feature_bed_files: %s\n" % str(feature_bed_files))
        self.ofh.write("feature_png_files: %s\n" % str(feature_png_files))
        self.ofh.write("flye_assembly_info_file: %s\n" % str(flye_assembly_info_file))
        self.ofh.write("gzipped: %s\n" % str(gzipped))
        self.ofh.write("genome_insertions_file: %s\n" % str(genome_insertions_file))
        self.ofh.write("illumina_forward_read_file: %s\n" % str(illumina_forward_read_file))
        self.ofh.write("illumina_reverse_read_file: %s\n" % str(illumina_reverse_read_file))
        self.ofh.write("kraken2_report_file: %s\n" % str(kraken2_report_file))
        self.ofh.write("kraken2_version: %s\n" % str(kraken2_version))
        self.ofh.write("lrn_risk_amr_file: %s\n" % str(lrn_risk_amr_file))
        self.ofh.write("lrn_risk_blacklist_file: %s\n" % str(lrn_risk_blacklist_file))
        self.ofh.write("lrn_risk_vf_file: %s\n" % str(lrn_risk_vf_file))
        self.ofh.write("minimap2_version: %s\n" % str(minimap2_version))
        self.ofh.write("mutation_regions_bed_file: %s\n" % str(mutation_regions_bed_file))
        self.ofh.write("mutation_regions_tsv_files: %s\n" % str(mutation_regions_tsv_files))
        self.ofh.write("ont_file: %s\n" % str(ont_file))
        self.ofh.write("pima_css: %s\n" % str(pima_css))
        self.ofh.write("plasmids_file: %s\n" % str(plasmids_file))
        self.ofh.write("quast_report_file: %s\n" % str(quast_report_file))
        self.ofh.write("read_type: %s\n" % str(read_type))
        self.ofh.write("reference_insertions_file: %s\n" % str(reference_insertions_file))
        self.ofh.write("samtools_version: %s\n" % str(samtools_version))
        self.ofh.write("varscan_version: %s\n" % str(varscan_version))

        # General
        self.doc = None
        self.report_md = 'pima_report.md'

        # Inputs
        self.amr_deletions_file = amr_deletions_file
        self.amr_matrix_files = amr_matrix_files
        self.analysis_name = analysis_name.split('_')[0]
        self.ofh.write("self.analysis_name: %s\n" % str(self.analysis_name))
        if assembler_version is None:
            self.assembler_version = 'assembler (version unknown)'
        else:
            if read_type == 'ont':
                # Assembler is flye.
                assembler_version = assembler_version.rstrip(' _assembly info_')
            else:
                # Assembler is spades.
                assembler_version = assembler_version.rstrip(' _contigs')
            self.assembler_version = re.sub('_', '.', assembler_version)
        self.assembly_fasta_file = assembly_fasta_file
        self.assembly_name = re.sub('_', '.', assembly_name.rstrip(' _consensus_'))
        if bedtools_version is None:
            self.bedtools_version = 'bedtools (version unknown)'
        else:
            self.bedtools_version = re.sub('_', '.', bedtools_version.rstrip(' _genome insertions'))
        if blastn_version is None:
            self.blastn_version = 'blastn (version unknown)'
        else:
            self.blastn_version = re.sub('_', '.', blastn_version.rstrip(' _features_'))
        self.circos_files = circos_files
        self.compute_sequence_length_file = compute_sequence_length_file
        self.contig_coverage_file = contig_coverage_file
        self.dbkey = dbkey
        self.dnadiff_snps_file = dnadiff_snps_file
        if dnadiff_version is None:
            self.dnadiff_version = 'dnadiff (version unknown)'
        else:
            self.dnadiff_version = re.sub('_', '.', dnadiff_version.rstrip(' _snps_'))
        self.errors_file = errors_file
        self.feature_bed_files = feature_bed_files
        self.feature_png_files = feature_png_files
        self.flye_assembly_info_file = flye_assembly_info_file
        self.gzipped = gzipped
        self.genome_insertions_file = genome_insertions_file
        self.illumina_forward_read_file = illumina_forward_read_file
        self.illumina_reverse_read_file = illumina_reverse_read_file
        self.kraken2_report_file = kraken2_report_file
        if kraken2_version is None:
            self.kraken2_version = 'kraken2 (version unknown)'
        else:
            self.kraken2_version = re.sub('_', '.', kraken2_version.rstrip(' _report_'))
        self.lrn_risk_amr_file = lrn_risk_amr_file
        self.lrn_risk_blacklist_file = lrn_risk_blacklist_file
        self.lrn_risk_vf_file = lrn_risk_vf_file
        if minimap2_version is None:
            self.minimap2_version = 'minimap2 (version unknown)'
        else:
            self.minimap2_version = re.sub('_', '.', minimap2_version)
        self.mutation_regions_bed_file = mutation_regions_bed_file
        self.mutation_regions_tsv_files = mutation_regions_tsv_files
        self.pima_css = pima_css
        self.plasmids_file = plasmids_file
        self.quast_report_file = quast_report_file
        self.read_type = read_type.upper()
        self.reference_insertions_file = reference_insertions_file
        self.reference_insertions_file = reference_insertions_file
        if samtools_version is None:
            self.samtools_version = 'samtools (version unknown)'
        else:
            self.samtools_version = re.sub('_', '.', samtools_version)
        if varscan_version is None:
            self.varscan_version = 'varscan (version unknown)'
        else:
            self.varscan_version = re.sub('_', '.', varscan_version)

        # Titles
        self.alignment_title = 'Comparison with reference'
        self.alignment_notes_title = 'Alignment notes'
        self.amr_matrix_title = 'AMR matrix'
        self.assembly_methods_title = 'Assembly'
        self.assembly_notes_title = 'Assembly notes'
        self.basecalling_title = 'Basecalling'
        self.basecalling_methods_title = 'Basecalling'
        self.contamination_methods_title = 'Contamination check'
        self.contig_alignment_title = 'Alignment vs. reference contigs'
        self.feature_title = 'Features found in the assembly'
        self.feature_methods_title = 'Feature annotation'
        self.feature_plot_title = 'Feature annotation plots'
        self.large_indel_title = 'Large insertions & deletions'
        self.lrn_risk_title = 'LRNRisk isolate classification'
        self.methods_title = 'Methods'
        self.mutation_errors_title = 'Errors finding mutations in the sample'
        self.mutation_title = 'Mutations found in the sample'
        self.mutation_methods_title = 'Mutation screening'
        self.plasmid_methods_title = 'Plasmid annotation'
        self.plasmid_title = 'Plasmid annotation'
        self.reference_genome_title = 'Reference genome'
        self.reference_methods_title = 'Reference comparison'
        self.snp_indel_title = 'SNPs and small indels'
        self.summary_title = 'Summary'

        # Methods
        self.methods = pandas.Series(dtype='float64')
        self.methods[self.contamination_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.assembly_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.reference_genome_title] = pandas.Series(dtype='float64')
        self.methods[self.reference_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.mutation_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.feature_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.plasmid_methods_title] = pandas.Series(dtype='float64')

        # Notes
        self.assembly_notes = pandas.Series(dtype=object)
        self.alignment_notes = pandas.Series(dtype=object)
        self.contig_alignment = pandas.Series(dtype=object)

        # Values
        self.assembly_size = 0
        self.contig_info = None
        self.feature_hits = pandas.Series(dtype='float64')
        self.ont_fast5 = None
        self.ont_file = ont_file
        self.ont_n50 = None
        self.ont_read_count = None
        # TODO: should the following be passed as a parameter?
        self.ont_coverage_min = 30
        # TODO: should the following be passed as a parameter?
        self.ont_n50_min = 2500

        if self.read_type == 'ONT':
            self.ont_raw_fastq = self.analysis_name
            self.ont_bases = 0
            self.illumina_bases = None
            self.illumina_fastq = None
            self.illumina_length_mean = None
            self.illumina_read_count = None
        else:
            self.illumina_fastq = self.analysis_name
            self.illumina_bases = 0
            self.illumina_length_mean = 0
            self.illumina_read_count = 0
            self.ont_bases = None
            self.ont_raw_fastq = None

        # Actions
        self.did_guppy_ont_fast5 = False
        self.did_qcat_ont_fastq = False
        self.ofh.write("self.read_type: %s\n" % str(self.read_type))
        if self.read_type == 'ONT':
            self.info_ont_fastq(self.ont_file)
        else:
            self.info_illumina_fastq([self.illumina_forward_read_file, self.illumina_reverse_read_file])
        self.load_contig_info()

    def run_command(self, command):
        self.ofh.write("\nXXXXXX In run_command, command:\n%s\n\n" % str(command))
        try:
            return re.split('\\n', subprocess.check_output(command, shell=True).decode('utf-8'))
        except Exception:
            message = 'Command %s failed: exiting...' % command
            sys.exit(message)

    def format_kmg(self, number, decimals=0):
        self.ofh.write("\nXXXXXX In format_kmg, number:\n%s\n" % str(number))
        self.ofh.write("XXXXXX In format_kmg, decimals:\n%s\n\n" % str(decimals))
        if number == 0:
            return '0'
        magnitude_powers = [10**9, 10**6, 10**3, 1]
        magnitude_units = ['G', 'M', 'K', '']
        for i in range(len(magnitude_units)):
            if number >= magnitude_powers[i]:
                magnitude_power = magnitude_powers[i]
                magnitude_unit = magnitude_units[i]
                return ('{:0.' + str(decimals) + 'f}').format(number / magnitude_power) + magnitude_unit

    def load_contig_info(self):
        self.contig_info = pandas.Series(dtype=object)
        self.contig_info[self.read_type] = pandas.read_csv(self.contig_coverage_file, header=None, index_col=None, sep='\t').sort_values(1, axis=0, ascending=False)
        self.contig_info[self.read_type].columns = ['contig', 'size', 'coverage']
        mean_coverage = (self.contig_info[self.read_type].iloc[:, 1] * self.contig_info[self.read_type].iloc[:, 2]).sum() / self.contig_info[self.read_type].iloc[:, 1].sum()
        if mean_coverage <= self.ont_coverage_min:
            warning = '%s mean coverage ({:.0f}X) is less than the recommended minimum ({:.0f}X).'.format(mean_coverage, self.ont_coverage_min) % self.read_type
            self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))
        # Report if some contigs have low coverage.
        low_coverage = self.contig_info[self.read_type].loc[self.contig_info[self.read_type]['coverage'] < self.ont_coverage_min, :]
        if low_coverage.shape[0] >= 0:
            for contig_i in range(low_coverage.shape[0]):
                warning = '%s coverage of {:s} ({:.0f}X) is less than the recommended minimum ({:.0f}X).'.format(low_coverage.iloc[contig_i, 0], low_coverage.iloc[contig_i, 2], self.ont_coverage_min) % self.read_type
                self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))
        # See if some contigs have anonymously low coverage.
        fold_coverage = self.contig_info[self.read_type]['coverage'] / mean_coverage
        low_coverage = self.contig_info[self.read_type].loc[fold_coverage < 1 / 5, :]
        if low_coverage.shape[0] >= 0:
            for contig_i in range(low_coverage.shape[0]):
                warning = '%s coverage of {:s} ({:.0f}X) is less than 1/5 the mean coverage ({:.0f}X).'.format(low_coverage.iloc[contig_i, 0], low_coverage.iloc[contig_i, 2], mean_coverage) % self.read_type
                self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))

    def load_fasta(self, fasta):
        sequence = pandas.Series(dtype=object)
        for contig in SeqIO.parse(fasta, 'fasta'):
            sequence[contig.id] = contig
        return sequence

    def load_assembly(self):
        self.assembly = self.load_fasta(self.assembly_fasta_file)
        self.num_assembly_contigs = len(self.assembly)
        self.assembly_size = self.format_kmg(sum([len(x) for x in self.assembly]), decimals=1)

    def info_illumina_fastq(self, illumina_read_files):
        self.ofh.write("\nXXXXXX In info_illumina_fastq\n\n")
        if self.gzipped:
            opener = 'gunzip -c'
        else:
            opener = 'cat'
        for fastq_file in illumina_read_files:
            command = ' '.join([opener,
                                fastq_file,
                                '| awk \'{getline;s += length($1);getline;getline;}END{print s/(NR/4)"\t"(NR/4)"\t"s}\''])
            output = self.run_command(command)
            self.ofh.write("output:\n%s\n" % str(output))
            self.ofh.write("re.split('\\t', self.run_command(command)[0]:\n%s\n" % str(re.split('\\t', self.run_command(command)[0])))
            values = []
            for i in re.split('\\t', self.run_command(command)[0]):
                if i == '':
                    values.append(float('nan'))
                else:
                    values.append(float(i))
            self.ofh.write("values:\n%s\n" % str(values))
            self.ofh.write("values[0]:\n%s\n" % str(values[0]))
            self.illumina_length_mean += values[0]
            self.ofh.write("values[1]:\n%s\n" % str(values[1]))
            self.illumina_read_count += int(values[1])
            self.ofh.write("values[2]:\n%s\n" % str(values[2]))
            self.illumina_bases += int(values[2])
        self.illumina_length_mean /= 2
        self.illumina_bases = self.format_kmg(self.illumina_bases, decimals=1)

    def start_doc(self):
        header_text = 'Analysis of ' + self.analysis_name
        self.doc = MdUtils(file_name=self.report_md, title=header_text)

    def add_table_of_contents(self):
        self.doc.create_marker(text_marker="TableOfContents")
        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_run_information(self):
        self.ofh.write("\nXXXXXX In add_run_information\n\n")
        self.doc.new_line()
        self.doc.new_header(1, 'Run information')
        # Tables in md.utils are implemented as a wrapping function.
        Table_list = [
            "Category",
            "Information",
            "Date",
            date.today(),
            "ONT FAST5",
            self.wordwrap_markdown(self.ont_fast5),
            "ONT FASTQ",
            self.wordwrap_markdown(self.ont_raw_fastq),
            "Illumina FASTQ",
            self.wordwrap_markdown(self.illumina_fastq),
            "Assembly",
            self.wordwrap_markdown(self.assembly_name),
            "Reference",
            self.wordwrap_markdown(self.dbkey),
        ]
        self.doc.new_table(columns=2, rows=7, text=Table_list, text_align='left')
        self.doc.new_line()
        # FIXME: the following doesn't work.
        # self.add_table_of_contents()
        self.doc.new_line()

    def add_ont_library_information(self):
        self.ofh.write("\nXXXXXX In add_ont_library_information\n\n")
        if self.ont_n50 is None:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'ONT library statistics')
        Table_List = [
            "Category",
            "Quantity",
            "ONT N50",
            '{:,}'.format(self.ont_n50),
            "ONT reads",
            '{:,}'.format(self.ont_read_count),
            "ONT bases",
            '{:s}'.format(self.ont_bases),
            "Illumina FASTQ",
            "N/A",
            "Assembly",
            self.wordwrap_markdown(self.assembly_name),
            "Reference",
            self.wordwrap_markdown(self.dbkey),
        ]
        self.doc.new_table(columns=2, rows=7, text=Table_List, text_align='left')
        self.doc.new_line()

    def add_illumina_library_information(self):
        self.ofh.write("\nXXXXXX In add_illumina_library_information\n\n")
        if self.illumina_length_mean is None:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'Illumina library statistics')
        Table_List = [
            "Illumina Info.",
            "Quantity",
            'Illumina mean length',
            '{:.1f}'.format(self.illumina_length_mean),
            'Illumina reads',
            '{:,}'.format(self.illumina_read_count),
            'Illumina bases',
            '{:s}'.format(self.illumina_bases)
        ]
        self.doc.new_table(columns=2, rows=4, text=Table_List, text_align='left')

    def evaluate_assembly(self):
        assembly_info = pandas.read_csv(self.compute_sequence_length_file, sep='\t', header=None)
        assembly_info.columns = ['contig', 'length']
        self.contig_sizes = assembly_info
        # Take a look at the number of contigs, their sizes,
        # and circularity.  Warn if things don't look good.
        if assembly_info.shape[0] > 4:
            warning = 'Assembly produced {:d} contigs, more than ususally expected; assembly may be fragmented'.format(assembly_info.shape[0])
            self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))
        small_contigs = assembly_info.loc[assembly_info['length'] <= 3000, :]
        if small_contigs.shape[0] > 0:
            warning = 'Assembly produced {:d} small contigs ({:s}); assembly may include spurious sequences.'.format(small_contigs.shape[0], ', '.join(small_contigs['contig']))
            self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))

    def add_assembly_information(self):
        self.ofh.write("\nXXXXXX In add_assembly_information\n\n")
        if self.assembly_fasta_file is None:
            return
        self.load_assembly()
        self.doc.new_line()
        self.doc.new_header(2, 'Assembly statistics')
        Table_List = [
            "Category",
            "Information",
            "Contigs",
            str(self.num_assembly_contigs),
            "Assembly size",
            str(self.assembly_size),
        ]
        self.doc.new_table(columns=2, rows=3, text=Table_List, text_align='left')

    def info_ont_fastq(self, fastq_file):
        self.ofh.write("\nXXXXXX In info_ont_fastq, fastq_file:\n%s\n\n" % str(fastq_file))
        opener = 'cat'
        if self.gzipped:
            opener = 'gunzip -c'
        else:
            opener = 'cat'
        command = ' '.join([opener,
                            fastq_file,
                            '| awk \'{getline;print length($0);s += length($1);getline;getline;}END{print "+"s}\'',
                            '| sort -gr',
                            '| awk \'BEGIN{bp = 0;f = 0}',
                            '{if(NR == 1){sub(/+/, "", $1);s=$1}else{bp += $1;if(bp > s / 2 && f == 0){n50 = $1;f = 1}}}',
                            'END{printf "%d\\t%d\\t%d\\n", n50, (NR - 1), s;exit}\''])
        result = list(re.split('\\t', self.run_command(command)[0]))
        if result[1] == '0':
            warning = 'No ONT reads found'
            self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))
        self.ont_n50, self.ont_read_count, ont_raw_bases = [int(i) for i in result]
        command = ' '.join([opener,
                            fastq_file,
                            '| awk \'{getline;print length($0);getline;getline;}\''])
        result = self.run_command(command)
        result = list(filter(lambda x: x != '', result))
        self.ont_bases = self.format_kmg(ont_raw_bases, decimals=1)
        if self.ont_n50 <= self.ont_n50_min:
            warning = 'ONT N50 (%s) is less than the recommended minimum (%s)' % (str(self.ont_n50), str(self.ont_n50_min))
            self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))

    def wordwrap_markdown(self, string):
        if string:
            if len(string) < 35:
                return string
            else:
                if '/' in string:
                    adjust = string.split('/')
                    out = ''
                    max = 35
                    for i in adjust:
                        out = out + '/' + i
                        if len(out) > max:
                            out += '<br>'
                            max += 35
                    return out
                else:
                    out = [string[i:i + 35] for i in range(0, len(string), 50)]
                    return '<br>'.join(out)
        else:
            return string

    def add_contig_info(self):
        self.ofh.write("\nXXXXXX In add_contig_info\n\n")
        if self.contig_info is None or self.read_type not in self.contig_info.index:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'Assembly coverage by ' + self.read_type)
        Table_List = ["Contig", "Length (bp)", "Coverage (X)"]
        formatted = self.contig_info[self.read_type].copy()
        formatted.iloc[:, 1] = formatted.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
        for i in range(self.contig_info[self.read_type].shape[0]):
            Table_List = Table_List + formatted.iloc[i, :].values.tolist()
        row_count = int(len(Table_List) / 3)
        self.doc.new_table(columns=3, rows=row_count, text=Table_List, text_align='left')

    def add_assembly_notes(self):
        self.ofh.write("\nXXXXXX In add_assembly_notes\n\n")
        if len(self.assembly_notes) == 0:
            return
        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()
        self.doc.new_header(2, self.assembly_notes_title)
        for note in self.assembly_notes:
            self.doc.new_line(note)

    def add_contamination(self):
        self.ofh.write("\nXXXXXX In add_contamination\n\n")
        if self.kraken2_report_file is None:
            return
        # Read in the Kraken fractions and pull out the useful parts
        kraken_fracs = pandas.read_csv(self.kraken2_report_file, delimiter='\t', header=None)
        kraken_fracs.index = kraken_fracs.iloc[:, 4].values
        kraken_fracs = kraken_fracs.loc[kraken_fracs.iloc[:, 3].str.match('[UG]1?'), :]
        kraken_fracs = kraken_fracs.loc[(kraken_fracs.iloc[:, 0] >= 1) | (kraken_fracs.iloc[:, 3] == 'U'), :]
        kraken_fracs = kraken_fracs.iloc[:, [0, 1, 3, 5]]
        kraken_fracs.columns = ['Fraction', 'Reads', 'Level', 'Taxa']
        kraken_fracs['Fraction'] = (kraken_fracs['Fraction'] / 100).round(4)
        kraken_fracs.sort_values(by='Fraction', inplace=True, ascending=False)
        kraken_fracs['Taxa'] = kraken_fracs['Taxa'].str.lstrip()
        self.doc.new_line()
        self.doc.new_header(2, 'Contamination check')
        self.doc.new_line(self.read_type + ' classifications')
        self.doc.new_line()
        Table_List = ["Percent of Reads", "Reads", "Level", "Label"]
        for index, row in kraken_fracs.iterrows():
            Table_List = Table_List + row.tolist()
        row_count = int(len(Table_List) / 4)
        self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')
        if self.contamination_methods_title not in self.methods:
            self.methods[self.contamination_methods_title] = ''
        method = '%s was used to assign the raw reads into taxa.' % self.kraken2_version.rstrip('report')
        self.methods[self.contamination_methods_title] = self.methods[self.contamination_methods_title].append(pandas.Series(method))

    def add_alignment(self):
        self.ofh.write("\nXXXXXX In add_alignment\n\n")
        if self.quast_report_file is not None:
            # Process quast values.
            quast_report = pandas.read_csv(self.quast_report_file, header=0, index_col=0, sep='\t')
            quast_mismatches = int(float(quast_report.loc['# mismatches per 100 kbp', :][0]) * (float(quast_report.loc['Total length (>= 0 bp)', :][0]) / 100000.))
            quast_indels = int(float(quast_report.loc['# indels per 100 kbp', :][0]) * (float(quast_report.loc['Total length (>= 0 bp)', :][0]) / 100000.))
            self.doc.new_line()
            self.doc.new_header(level=2, title=self.alignment_title)
            self.doc.new_line()
            self.doc.new_header(level=3, title=self.snp_indel_title)
            Table_1 = [
                "Category",
                "Quantity",
                'SNPs',
                '{:,}'.format(quast_mismatches),
                'Small indels',
                '{:,}'.format(quast_indels)
            ]
        self.doc.new_table(columns=2, rows=3, text=Table_1, text_align='left')
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()
        # TODO: self.alignment_notes is not currently populated.
        if len(self.alignment_notes) > 0:
            self.doc.new_header(level=3, title=self.alignment_notes_title)
            for note in self.alignment_notes:
                self.doc.new_line(note)
        if len(self.circos_files) > 0:
            # Add circos PNG files.
            for circos_file in self.circos_files:
                contig = os.path.basename(circos_file)
                contig_title = 'Alignment to %s' % contig
                self.doc.new_line()
                self.doc.new_header(level=3, title=contig_title)
                self.doc.new_line('Blue color indicates query sequences aligned to the reference sequence, which is shown in yellow')
                self.doc.new_line(self.doc.new_inline_image(text='contig_title', path=os.path.abspath(circos_file)))
                self.doc.new_line('<div style="page-break-after: always;"></div>')
                self.doc.new_line()
        if self.dbkey == 'ref_genome':
            headers = ["* Chromosome - NC_007530.2 Bacillus anthracis str. 'Ames Ancestor', complete sequence",
                       "* pXO1 - NC_007322.2 Bacillus anthracis str. 'Ames Ancestor' plasmid pXO1, complete sequence",
                       "* pXO2 - NC_007323.3 Bacillus anthracis str. 'Ames Ancestor' plasmid pXO2, complete sequence"]
            method = '\n'.join(headers)
            self.methods[self.reference_genome_title] = self.methods[self.reference_genome_title].append(pandas.Series(method))
        method = 'The genome assembly was aligned against the reference sequence using %s.' % self.dnadiff_version
        self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(pandas.Series(method))

    def add_features(self):
        self.ofh.write("\nXXXXXX In add_features\n\n")
        if len(self.feature_bed_files) == 0:
            return
        for bbf in self.feature_bed_files:
            if os.path.getsize(bbf) > 0:
                best = pandas.read_csv(filepath_or_buffer=bbf, sep='\t', header=None)
                self.feature_hits[os.path.basename(bbf)] = best
        if len(self.feature_hits) == 0:
            return
        self.ofh.write("self.feature_hits: %s\n" % str(self.feature_hits))
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.feature_title)
        for feature_name in self.feature_hits.index.tolist():
            self.ofh.write("feature_name: %s\n" % str(feature_name))
            features = self.feature_hits[feature_name].copy()
            self.ofh.write("features: %s\n" % str(features))
            if features.shape[0] == 0:
                continue
            features.iloc[:, 1] = features.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            features.iloc[:, 2] = features.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
            self.doc.new_line()
            self.doc.new_header(level=3, title=feature_name)
            if (features.shape[0] == 0):
                continue
            for contig in pandas.unique(features.iloc[:, 0]):
                self.ofh.write("contig: %s\n" % str(contig))
                self.doc.new_line(contig)
                contig_features = features.loc[(features.iloc[:, 0] == contig), :]
                self.ofh.write("contig_features: %s\n" % str(contig_features))
                Table_List = ['Start', 'Stop', 'Feature', 'Identity (%)', 'Strand']
                for i in range(contig_features.shape[0]):
                    self.ofh.write("i: %s\n" % str(i))
                    feature = contig_features.iloc[i, :].copy(deep=True)
                    self.ofh.write("feature: %s\n" % str(feature))
                    feature[4] = '{:.3f}'.format(feature[4])
                    self.ofh.write("feature[1:].values.tolist(): %s\n" % str(feature[1:].values.tolist()))
                    Table_List = Table_List + feature[1:].values.tolist()
                self.ofh.write("Table_List: %s\n" % str(Table_List))
                row_count = int(len(Table_List) / 5)
                self.ofh.write("row_count: %s\n" % str(row_count))
                self.doc.new_line()
                self.ofh.write("Before new_table, len(Table_List):: %s\n" % str(len(Table_List)))
                self.doc.new_table(columns=5, rows=row_count, text=Table_List, text_align='left')
        blastn_version = 'The genome assembly was queried for features using %s.' % self.blastn_version
        bedtools_version = 'Feature hits were clustered using %s and the highest scoring hit for each cluster was reported.' % self.bedtools_version
        method = '%s  %s' % (blastn_version, bedtools_version)
        self.methods[self.feature_methods_title] = self.methods[self.feature_methods_title].append(pandas.Series(method))

    def add_feature_plots(self):
        self.ofh.write("\nXXXXXX In add_feature_plots\n\n")
        if len(self.feature_png_files) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title='Feature Plots')
        self.doc.new_paragraph('Only contigs with features are shown')
        for feature_png_file in self.feature_png_files:
            self.doc.new_line(self.doc.new_inline_image(text='Analysis', path=os.path.abspath(feature_png_file)))

    def add_mutations(self):
        self.ofh.write("\nXXXXXX In add_mutations\n\n")
        if len(self.mutation_regions_tsv_files) == 0:
            return
        try:
            mutation_regions = pandas.read_csv(self.mutation_regions_bed_file, sep='\t', header=0, index_col=False)
        except Exception:
            # Likely an empty file.
            return
        amr_mutations = pandas.Series(dtype=object)
        for region_i in range(mutation_regions.shape[0]):
            region = mutation_regions.iloc[region_i, :]
            region_name = str(region['name'])
            self.ofh.write("Processing mutations for region %s\n" % region_name)
            region_mutations_tsv_name = '%s_mutations.tsv' % region_name
            if region_mutations_tsv_name not in self.mutation_regions_tsv_files:
                continue
            region_mutations_tsv = self.mutation_regions_tsv_files[region_mutations_tsv_name]
            try:
                region_mutations = pandas.read_csv(region_mutations_tsv, sep='\t', header=0, index_col=False)
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
        if (amr_mutations.shape[0] > 0):
            # Report the mutations.
            self.doc.new_line()
            self.doc.new_header(level=2, title=self.mutation_title)
            for region_name in amr_mutations.index.tolist():
                region_mutations = amr_mutations[region_name].copy()
                self.doc.new_line()
                self.doc.new_header(level=3, title=region_name)
                if (region_mutations.shape[0] == 0):
                    self.doc.append('None')
                    continue
                region_mutations.iloc[:, 1] = region_mutations.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
                Table_List = ['Reference contig', 'Position', 'Reference', 'Alternate', 'Drug', 'Note']
                for i in range(region_mutations.shape[0]):
                    Table_List = Table_List + region_mutations.iloc[i, [0, 1, 3, 4, 5, 6]].values.tolist()
                row_count = int(len(Table_List) / 6)
                self.doc.new_table(columns=6, rows=row_count, text=Table_List, text_align='left')
        if os.path.getsize(self.errors_file) > 0:
            # Report the errors encountered when attempting
            # to find mutations in the sample.
            self.doc.new_line()
            self.doc.new_header(level=2, title=self.mutation_errors_title)
            with open(self.errors_file, 'r') as efh:
                for i, line in enumerate(efh):
                    line = line.strip()
                    if line:
                        self.doc.new_line('* %s' % line)
        method = '%s reads were mapped to the reference sequence using %s.' % (self.read_type, self.minimap2_version)
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pandas.Series(method))
        method = 'Mutations were identified using %s and %s.' % (self.samtools_version, self.varscan_version)
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pandas.Series(method))

    def add_amr_matrix(self):
        self.ofh.write("\nXXXXXX In add_amr_matrix\n\n")
        # Make sure that we have an AMR matrix to plot
        if len(self.amr_matrix_files) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.amr_matrix_title)
        self.doc.new_line('AMR genes and mutations with their corresponding drugs')
        for amr_matrix_file in self.amr_matrix_files:
            self.doc.new_line(self.doc.new_inline_image(text='AMR genes and mutations with their corresponding drugs',
                                                        path=os.path.abspath(amr_matrix_file)))

    def add_large_indels(self):
        self.ofh.write("\nXXXXXX In add_large_indels\n\n")
        large_indels = pandas.Series(dtype='float64')
        # Pull in insertions.
        try:
            reference_insertions = pandas.read_csv(filepath_or_buffer=self.reference_insertions_file, sep='\t', header=None)
        except Exception:
            reference_insertions = pandas.DataFrame()
        try:
            genome_insertions = pandas.read_csv(filepath_or_buffer=self.genome_insertions_file, sep='\t', header=None)
        except Exception:
            genome_insertions = pandas.DataFrame()
        large_indels['Reference insertions'] = reference_insertions
        large_indels['Query insertions'] = genome_insertions
        # Pull in deletions.
        try:
            amr_deletions = pandas.read_csv(filepath_or_buffer=self.amr_deletion_file, sep='\t', header=None)
        except Exception:
            amr_deletions = pandas.DataFrame()
        if amr_deletions.shape[0] > 0:
            amr_deletions.columns = ['contig', 'start', 'stop', 'name', 'type', 'drug', 'note']
            amr_deletions = amr_deletions.loc[amr_deletions['type'].isin(['large-deletion', 'any']), :]
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.large_indel_title)
        self.doc.new_line('This section is informative only when your isolates were identified as *Bacillus anthracis* strains')
        for genome in ['Reference insertions', 'Query insertions']:
            genome_indels = large_indels[genome].copy()
            self.doc.new_line()
            self.doc.new_header(level=3, title=genome)
            if (genome_indels.shape[0] == 0):
                continue
            genome_indels.iloc[:, 1] = genome_indels.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            genome_indels.iloc[:, 2] = genome_indels.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
            genome_indels.iloc[:, 3] = genome_indels.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
            Table_List = [
                'Reference contig', 'Start', 'Stop', 'Size (bp)'
            ]
            for i in range(genome_indels.shape[0]):
                Table_List = Table_List + genome_indels.iloc[i, :].values.tolist()
            row_count = int(len(Table_List) / 4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')
        method = 'Large insertions or deletions were found as the complement of aligned regions using %s.' % self.bedtools_version
        self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(pandas.Series(method))
        self.doc.new_line()
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_lrn_risk_info(self):
        self.ofh.write("\nXXXXXX In add_lrn_risk_info\n\n")
        if self.lrn_risk_amr_file is None and self.lrn_risk_blacklist_file is None and self.lrn_risk_vf_file is None:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.lrn_risk_title)
        # Process self.lrn_risk_amr_file.
        try:
            lrn_risk_amr = pandas.read_csv(filepath_or_buffer=self.lrn_risk_amr_file, sep='\t', header=0)
        except Exception:
            lrn_risk_amr = pandas.DataFrame()
        if lrn_risk_amr.shape[0] > 0:
            self.doc.new_line()
            self.doc.new_header(level=2, title="AMR Determinant Distribution")
            self.doc.new_line()
            Table_List = ["Gene", "Contig", "% Identity", "% Coverage", "E-Value", "Annotation", "Comparison to Publicly Available Genomes"]
            for index, row in lrn_risk_amr.iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List) / 7)
            self.doc.new_table(columns=7, rows=row_count, text=Table_List, text_align='left')
        # Process self.lrn_risk_blacklist_file.
        try:
            lrn_risk_blacklist = pandas.read_csv(filepath_or_buffer=self.lrn_risk_blacklist_file, sep='\t', header=0)
        except Exception:
            lrn_risk_blacklist = pandas.DataFrame()
        if lrn_risk_blacklist.shape[0] > 0:
            self.doc.new_line()
            self.doc.new_header(level=2, title="Blacklisted High-risk Virulence Factors")
            self.doc.new_line()
            Table_List = ["Blacklisted Gene", "Reason", "Risk Category"]
            for index, row in lrn_risk_blacklist.iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List) / 3)
            self.doc.new_table(columns=3, rows=row_count, text=Table_List, text_align='left')
        # Process self.lrn_risk_vf_file.
        try:
            lrn_risk_vf = pandas.read_csv(filepath_or_buffer=self.lrn_risk_vf_file, sep='\t', header=0)
        except Exception:
            lrn_risk_vf = pandas.DataFrame()
        if lrn_risk_vf.shape[0] > 0:
            self.doc.new_line()
            self.doc.new_header(level=2, title="Virulence Factor Distribution")
            self.doc.new_line()
            Table_List = ["Gene", "Contig", "% Identity", "% Coverage", "E-Value", "Annotation", "Comparison to Publicly Available Genomes"]
            for index, row in lrn_risk_vf.iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List) / 7)
            self.doc.new_table(columns=7, rows=row_count, text=Table_List, text_align='left')
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_plasmids(self):
        try:
            plasmids = pandas.read_csv(filepath_or_buffer=self.plasmids_file, sep='\t', header=0)
        except Exception:
            return
        plasmids = plasmids.copy()
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.plasmid_title)
        if (plasmids.shape[0] == 0):
            self.doc.new_line('None')
            return
        plasmids.iloc[:, 3] = plasmids.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 4] = plasmids.iloc[:, 4].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 5] = plasmids.iloc[:, 5].apply(lambda x: '{:,}'.format(x))
        Table_List = ['Genome contig', 'Plasmid hit', 'Plasmid acc.', 'Contig size', 'Aliged', 'Plasmid size']
        for i in range(plasmids.shape[0]):
            Table_List = Table_List + plasmids.iloc[i, 0:6].values.tolist()
        row_count = int(len(Table_List) / 6)
        self.doc.new_table(columns=6, rows=row_count, text=Table_List, text_align='left')
        method = 'The plasmid reference database was queried against the genome assembly using %s.' % self.minimap2_version
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))
        method = 'The resulting BAM was converted to a PSL using a custom version of sam2psl.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))
        method = 'Plasmid-to-genome hits were resolved using the pChunks algorithm.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))

    def add_methods(self):
        self.ofh.write("\nXXXXXX In add_methods\n\n")
        if len(self.methods) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.methods_title)
        for methods_section in self.methods.index.tolist():
            if self.methods[methods_section] is None or len(self.methods[methods_section]) == 0:
                continue
            self.doc.new_line()
            self.doc.new_header(level=3, title=methods_section)
            self.doc.new_paragraph(' '.join(self.methods[methods_section]))
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()

    def add_summary(self):
        self.ofh.write("\nXXXXXX In add_summary\n\n")
        # Add summary title
        self.doc.new_header(level=1, title=self.summary_title)
        # First section of Summary
        self.doc.new_header(level=1, title='CDC Advisory')
        self.doc.new_paragraph(CDC_ADVISORY)
        self.doc.new_line()
        self.add_run_information()
        self.add_ont_library_information()
        methods = []
        if self.did_guppy_ont_fast5:
            methods += ['ONT reads were basecalled using guppy']
        if self.did_qcat_ont_fastq:
            methods += ['ONT reads were demultiplexed and trimmed using qcat']
        self.methods[self.basecalling_methods_title] = pandas.Series(methods)
        self.add_illumina_library_information()
        self.add_assembly_information()
        self.add_contig_info()
        self.evaluate_assembly()
        if self.assembler_version is not None:
            if self.read_type == 'ONT':
                method = 'ONT reads were assembled using %s' % self.assembler_version
                self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pandas.Series(method))
                # Pull in the assembly summary and look at the coverage.
                assembly_info = pandas.read_csv(self.flye_assembly_info_file, header=0, index_col=0, sep='\t')
                # Look for non-circular contigs.
                open_contigs = assembly_info.loc[assembly_info['circ.'] == 'N', :]
                if open_contigs.shape[0] > 0:
                    open_contig_ids = open_contigs.index.values
                    warning = 'Flye reported {:d} open contigs ({:s}); assembly may be incomplete.'.format(open_contigs.shape[0], ', '.join(open_contig_ids))
                    self.assembly_notes = self.assembly_notes.append(pandas.Series(warning))
                else:
                    method = 'Illumina reads were assembled using %s' % self.assembler_version
                method = 'The genome assembly was polished using ONT reads and medaka.'
                self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pandas.Series(method))
        self.add_assembly_notes()

    def make_tex(self):
        self.doc.new_table_of_contents(table_title='detailed run information', depth=2, marker="tableofcontents")
        text = self.doc.file_data_text
        text = text.replace("##--[", "")
        text = text.replace("]--##", "")
        self.doc.file_data_text = text
        self.doc.create_md_file()

    def make_report(self):
        self.ofh.write("\nXXXXXX In make_report\n\n")
        self.start_doc()
        self.add_summary()
        self.add_contamination()
        self.add_alignment()
        self.add_features()
        self.add_feature_plots()
        self.add_mutations()
        self.add_large_indels()
        self.add_plasmids()
        self.add_amr_matrix()
        self.add_lrn_risk_info()
        # self.add_snps()
        self.add_methods()
        self.make_tex()
        # It took me quite a long time to find out that the value of the -t
        # (implied) argument in the following command must be 'html' instead of
        # the more logical 'pdf'.  see the answer from snsn in this thread:
        # https://github.com/jessicategner/pypandoc/issues/186
        self.ofh.write("\nXXXXX In make_report, calling pypandoc.convert_file...\n\n")
        pypandoc.convert_file(self.report_md,
                              'html',
                              extra_args=['--pdf-engine=weasyprint', '-V', '-css=%s' % self.pima_css],
                              outputfile='pima_report.pdf')
        self.ofh.close()


parser = argparse.ArgumentParser()

parser.add_argument('--amr_deletions_file', action='store', dest='amr_deletions_file', help='AMR deletions BED file')
parser.add_argument('--amr_matrix_png_dir', action='store', dest='amr_matrix_png_dir', help='Directory of AMR matrix PNG files')
parser.add_argument('--analysis_name', action='store', dest='analysis_name', help='Sample identifier')
parser.add_argument('--assembler_version', action='store', dest='assembler_version', default=None, help='Assembler version string')
parser.add_argument('--assembly_fasta_file', action='store', dest='assembly_fasta_file', help='Assembly fasta file')
parser.add_argument('--assembly_name', action='store', dest='assembly_name', help='Assembly identifier')
parser.add_argument('--bedtools_version', action='store', dest='bedtools_version', default=None, help='Bedtools version string')
parser.add_argument('--blastn_version', action='store', dest='blastn_version', default=None, help='Blastn version string')
parser.add_argument('--circos_png_dir', action='store', dest='circos_png_dir', help='Directory of circos PNG files')
parser.add_argument('--compute_sequence_length_file', action='store', dest='compute_sequence_length_file', help='Comnpute sequence length tabular file')
parser.add_argument('--contig_coverage_file', action='store', dest='contig_coverage_file', help='Contig coverage TSV file')
parser.add_argument('--dbkey', action='store', dest='dbkey', help='Reference genome identifier')
parser.add_argument('--dnadiff_snps_file', action='store', dest='dnadiff_snps_file', help='DNAdiff snps tabular file')
parser.add_argument('--dnadiff_version', action='store', dest='dnadiff_version', default=None, help='DNAdiff version string')
parser.add_argument('--errors_file', action='store', dest='errors_file', default=None, help='AMR mutations errors encountered txt file')
parser.add_argument('--feature_bed_dir', action='store', dest='feature_bed_dir', help='Directory of best feature hits bed files')
parser.add_argument('--feature_png_dir', action='store', dest='feature_png_dir', help='Directory of best feature hits png files')
parser.add_argument('--flye_assembly_info_file', action='store', dest='flye_assembly_info_file', default=None, help='Flye assembly info tabular file')
parser.add_argument('--genome_insertions_file', action='store', dest='genome_insertions_file', help='Genome insertions BED file')
parser.add_argument('--gzipped', action='store_true', dest='gzipped', default=False, help='Sample(s) is/are gzipped')
parser.add_argument('--illumina_forward_read_file', action='store', dest='illumina_forward_read_file', help='Illumina forward read file')
parser.add_argument('--illumina_reverse_read_file', action='store', dest='illumina_reverse_read_file', help='Illumina reverse read file')
parser.add_argument('--kraken2_report_file', action='store', dest='kraken2_report_file', default=None, help='kraken2 report file')
parser.add_argument('--kraken2_version', action='store', dest='kraken2_version', default=None, help='kraken2 version string')
parser.add_argument('--lrn_risk_amr_file', action='store', dest='lrn_risk_amr_file', default=None, help='LRN RISK AMR TSV file')
parser.add_argument('--lrn_risk_blacklist_file', action='store', dest='lrn_risk_blacklist_file', default=None, help='LRN RISK blacklist TSV file')
parser.add_argument('--lrn_risk_vf_file', action='store', dest='lrn_risk_vf_file', default=None, help='LRN RISK virulence factors TSV file')
parser.add_argument('--minimap2_version', action='store', dest='minimap2_version', default=None, help='minimap2 version string')
parser.add_argument('--mutation_regions_bed_file', action='store', dest='mutation_regions_bed_file', help='AMR mutation regions BRD file')
parser.add_argument('--mutation_regions_dir', action='store', dest='mutation_regions_dir', help='Directory of mutation regions TSV files')
parser.add_argument('--ont_file', action='store', dest='ont_file', help='ONT single read file')
parser.add_argument('--pima_css', action='store', dest='pima_css', help='PIMA css stypesheet')
parser.add_argument('--plasmids_file', action='store', dest='plasmids_file', default=None, help='pChunks plasmids TSV file')
parser.add_argument('--quast_report_file', action='store', dest='quast_report_file', help='Quast report tabular file')
parser.add_argument('--read_type', action='store', dest='read_type', help='Sample read type (ONT or Illumina)')
parser.add_argument('--reference_insertions_file', action='store', dest='reference_insertions_file', help='Reference insertions BED file')
parser.add_argument('--samtools_version', action='store', dest='samtools_version', default=None, help='Samtools version string')
parser.add_argument('--varscan_version', action='store', dest='varscan_version', default=None, help='Varscan version string')

args = parser.parse_args()

# Prepare the AMR matrix PNG files.
amr_matrix_files = []
for file_name in sorted(os.listdir(args.amr_matrix_png_dir)):
    file_path = os.path.abspath(os.path.join(args.amr_matrix_png_dir, file_name))
    amr_matrix_files.append(file_path)
# Prepare the circos PNG files.
circos_files = []
for file_name in sorted(os.listdir(args.circos_png_dir)):
    file_path = os.path.abspath(os.path.join(args.circos_png_dir, file_name))
    circos_files.append(file_path)
# Prepare the features BED files.
feature_bed_files = []
for file_name in sorted(os.listdir(args.feature_bed_dir)):
    file_path = os.path.abspath(os.path.join(args.feature_bed_dir, file_name))
    feature_bed_files.append(file_path)
# Prepare the features PNG files.
feature_png_files = []
for file_name in sorted(os.listdir(args.feature_png_dir)):
    file_path = os.path.abspath(os.path.join(args.feature_png_dir, file_name))
    feature_png_files.append(file_path)
# Prepare the mutation regions TSV files.
mutation_regions_files = []
for file_name in sorted(os.listdir(args.mutation_regions_dir)):
    file_path = os.path.abspath(os.path.join(args.feature_png_dir, file_name))
    mutation_regions_files.append(file_path)

markdown_report = PimaReport(args.analysis_name,
                             args.amr_deletions_file,
                             amr_matrix_files,
                             args.assembler_version,
                             args.assembly_fasta_file,
                             args.assembly_name,
                             args.bedtools_version,
                             args.blastn_version,
                             circos_files,
                             args.compute_sequence_length_file,
                             args.contig_coverage_file,
                             args.dbkey,
                             args.dnadiff_snps_file,
                             args.dnadiff_version,
                             args.errors_file,
                             feature_bed_files,
                             feature_png_files,
                             args.flye_assembly_info_file,
                             args.genome_insertions_file,
                             args.gzipped,
                             args.illumina_forward_read_file,
                             args.illumina_reverse_read_file,
                             args.kraken2_report_file,
                             args.kraken2_version,
                             args.lrn_risk_amr_file,
                             args.lrn_risk_blacklist_file,
                             args.lrn_risk_vf_file,
                             args.minimap2_version,
                             args.mutation_regions_bed_file,
                             mutation_regions_files,
                             args.ont_file,
                             args.pima_css,
                             args.plasmids_file,
                             args.quast_report_file,
                             args.read_type,
                             args.reference_insertions_file,
                             args.samtools_version,
                             args.varscan_version)
markdown_report.make_report()
