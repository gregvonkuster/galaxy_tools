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
from mdutils.tools import TableOfContents

CDC_ADVISORY = 'The analysis and report presented here should be treated as preliminary.  Please contact the CDC/BDRD with any results regarding _Bacillus anthracis_.'


class PimaReport:

    def __init__(self, analysis_name, assembly_fasta_file, assembly_name, contig_coverage_file, dbkey, illumina_fastq_file, gzipped, pima_css):
        self.ofh = open("process_log.txt", "w")

        self.ofh.write("analysis_name: %s\n" % str(analysis_name))
        self.ofh.write("assembly_fasta_file: %s\n" % str(assembly_fasta_file))
        self.ofh.write("assembly_name: %s\n" % str(assembly_name))
        self.ofh.write("contig_coverage_file: %s\n" % str(contig_coverage_file))
        self.ofh.write("dbkey: %s\n" % str(dbkey))
        self.ofh.write("illumina_fastq_file: %s\n" % str(illumina_fastq_file))
        self.ofh.write("gzipped: %s\n" % str(gzipped))
        self.ofh.write("pima_css: %s\n" % str(pima_css))

        # General
        self.doc = None
        self.report_md = 'pima_report.md'

        # Inputs
        self.analysis_name = analysis_name
        self.assembly_fasta_file = assembly_fasta_file
        self.assembly_name = assembly_name
        self.contig_coverage_file = contig_coverage_file
        self.dbkey = dbkey
        self.illumina_fastq_file = illumina_fastq_file
        self.gzipped = gzipped
        self.read_type = 'Illumina'
        self.ont_bases = None
        self.ont_n50 = None
        self.ont_read_count = None
        self.pima_css = pima_css

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
        self.methods_title = 'Methods'
        self.mutation_title = 'Mutations found in the sample'
        self.mutation_methods_title = 'Mutation screening'
        self.plasmid_methods_title = 'Plasmid annotation'
        self.reference_methods_title = 'Reference comparison'
        self.snp_indel_title = 'SNPs and small indels'
        #self.summary_title = 'Summary'
        self.summary_title = 'Analysis of %s' % analysis_name

        # Methods
        self.methods = pandas.Series(dtype='float64')
        self.methods[self.contamination_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.assembly_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.reference_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.mutation_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.feature_methods_title] = pandas.Series(dtype='float64')
        self.methods[self.plasmid_methods_title] = pandas.Series(dtype='float64')

        # Contamination
        self.kraken_fracs = pandas.Series(dtype=object)

        # Notes
        self.assembly_notes = pandas.Series(dtype=object)
        self.alignment_notes = pandas.Series(dtype=object)
        self.contig_alignment = pandas.Series(dtype=object)

        # Values
        self.assembly_size = 0
        self.contig_info = None
        self.did_flye_ont_fastq = False
        self.did_medaka_ont_assembly = False
        self.feature_hits = pandas.Series(dtype=object)
        self.feature_plots = pandas.Series(dtype=object)
        self.illumina_length_mean = 0
        self.illumina_read_count = 0
        self.illumina_bases = 0
        self.mean_coverage = 0
        self.num_assembly_contigs = 0

        # Actions
        self.did_guppy_ont_fast5 = False
        self.did_qcat_ont_fastq = False
        self.info_illumina_fastq()
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
        self.mean_coverage = (self.contig_info[self.read_type].iloc[:, 1] * self.contig_info[self.read_type].iloc[:, 2]).sum() / self.contig_info[self.read_type].iloc[:, 1].sum()

    def load_fasta(self, fasta):
        sequence = pandas.Series(dtype=object)
        for contig in SeqIO.parse(fasta, 'fasta'):
            sequence[contig.id] = contig
        return sequence

    def load_assembly(self):
        self.assembly = self.load_fasta(self.assembly_fasta_file)
        self.num_assembly_contigs = len(self.assembly)
        for i in self.assembly:
            self.assembly_size += len(i.seq)
        self.assembly_size = self.format_kmg(self.assembly_size, decimals=1)

    def info_illumina_fastq(self):
        self.ofh.write("\nXXXXXX In info_illumina_fastq\n\n")
        if self.gzipped:
            opener = 'gunzip -c'
        else:
            opener = 'cat'
        command = ' '.join([opener,
                            self.illumina_fastq_file,
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
        # The original PIMA code inserts self.illumina_fastq into
        # a list for no apparent reason.  We don't do that here.
        # self.illumina_length_mean /= len(self.illumina_fastq)
        self.illumina_length_mean /= 1
        self.illumina_bases = self.format_kmg(self.illumina_bases, decimals=1)

    def start_doc(self):
        #header_text = 'Analysis of %s' % self.analysis_name
        self.doc = MdUtils(file_name=self.report_md, title='')

    def add_tableOfContents(self):
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
            "N/A",
            "ONT FASTQ",
            "N/A",
            "Illumina FASTQ",
            self.wordwrap_markdown(self.analysis_name),
            "Assembly",
            self.wordwrap_markdown(self.assembly_name),
            "Reference",
            self.wordwrap_markdown(self.dbkey),
        ]
        self.doc.new_table(columns=2, rows=7, text=Table_list, text_align='left')
        self.doc.new_line()
        self.add_tableOfContents()
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
            self.wordwrap_markdown(self.illumina_fastq_file),
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
            self.error_out('No ONT reads found')
        ont_n50, ont_read_count, ont_raw_bases = [int(i) for i in result]

        command = ' '.join([opener,
                            fastq_file,
                            '| awk \'{getline;print length($0);getline;getline;}\''])
        result = self.run_command(command)
        result = list(filter(lambda x: x != '', result))
        ont_read_lengths = [int(i) for i in result]

        return ([ont_n50, ont_read_count, ont_raw_bases, ont_read_lengths])

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
        if self.contig_info is None:
            return
        for method in ['ONT', 'Illumina']:
            if method not in self.contig_info.index:
                continue
            self.doc.new_line()
            self.doc.new_header(2, 'Assembly coverage by ' + method)
            Table_List = ["Contig", "Length (bp)", "Coverage (X)"]
            formatted = self.contig_info[method].copy()
            formatted.iloc[:, 1] = formatted.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            for i in range(self.contig_info[method].shape[0]):
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
        #for note in self.analysis.assembly_notes:
        #    self.doc.new_line(note)

    def add_contamination(self):
        self.ofh.write("\nXXXXXX In add_contamination\n\n")
        if len(self.kraken_fracs) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(2, 'Contamination check')
        for read_type, kraken_fracs in self.kraken_fracs.iteritems():
            self.doc.new_line(read_type + ' classifications')
            self.doc.new_line()
            Table_List = ["Percent of Reads", "Reads", "Level", "Label"]
            for index, row in kraken_fracs.iterrows():
                Table_List = Table_List + row.tolist()
            row_count = int(len(Table_List) / 4)
            self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')
            if self.contamination_methods_title not in self.methods:
                self.methods[self.contamination_methods_title] = ''
        method = 'Kraken2 was used to assign the raw reads into taxa.'
        self.methods[self.contamination_methods_title] = self.methods[self.contamination_methods_title].append(pandas.Series(method))

    def add_alignment(self):
        self.ofh.write("\nXXXXXX In add_alignment\n\n")
        if len(self.contig_alignment) > 0:
            alignments = self.contig_alignment
        else:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.alignment_title)
        self.doc.new_line()
        self.doc.new_header(level=3, title=self.snp_indel_title)
        Table_1 = [
            "Category",
            "Quantity",
            'SNPs',
            #'{:,}'.format(self.analysis.quast_mismatches),
            'NA'
            'Small indels',
            #'{:,}'.format(self.analysis.quast_indels)
            'NA'
        ]
        self.doc.new_table(columns=2, rows=3, text=Table_1, text_align='left')
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()
        if len(self.alignment_notes) > 0:
            self.doc.new_header(level=3, title=self.alignment_notes_title)
            for note in self.alignment_notes:
                self.doc.new_line(note)
        for contig in alignments.index.tolist():
            contig_title = 'Alignment to %s' % contig
            image_png = alignments[contig]
            self.doc.new_line()
            self.doc.new_header(level=3, title=contig_title)
            self.doc.new_line(self.doc.new_inline_image(text='contig_title', path=os.path.abspath(image_png)))
            self.doc.new_line('<div style="page-break-after: always;"></div>')
            self.doc.new_line()
        method = 'The genome assembly was aligned against the reference sequencing using dnadiff.'
        self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(pandas.Series(method))

    def add_features(self):
        self.ofh.write("\nXXXXXX In add_features\n\n")
        if len(self.feature_hits) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.feature_title)
        for feature_name in self.feature_hits.index.tolist():
            features = self.feature_hits[feature_name].copy()
            if features.shape[0] == 0:
                continue
            features.iloc[:, 1] = features.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
            features.iloc[:, 2] = features.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
            self.doc.new_line()
            self.doc.new_header(level=3, title=feature_name)
            if (features.shape[0] == 0):
                continue
            for contig in pandas.unique(features.iloc[:, 0]):
                self.doc.new_line(contig)
                contig_features = features.loc[(features.iloc[:, 0] == contig), :]
                Table_List = ['Start', 'Stop', 'Feature', 'Identity (%)', 'Strand']
                for i in range(contig_features.shape[0]):
                    feature = contig_features.iloc[i, :].copy(deep=True)
                    feature[4] = '{:.3f}'.format(feature[4])
                    Table_List = Table_List + feature[1:].values.tolist()
                row_count = int(len(Table_List) / 5)
                self.doc.new_line()
                self.doc.new_table(columns=5, rows=row_count, text=Table_List, text_align='left')
        blastn_version = 'The genome assembly was queried for features using blastn.'
        bedtools_version = 'Feature hits were clustered using bedtools and the highest scoring hit for each cluster was reported.'
        method = '%s  %s' % (blastn_version, bedtools_version)
        self.methods[self.feature_methods_title] = self.methods[self.feature_methods_title].append(pandas.Series(method))
        self.methods[self.feature_methods_title] = self.methods[self.feature_methods_title].append([method])

    def add_feature_plots(self):
        self.ofh.write("\nXXXXXX In add_feature_plots\n\n")
        if len(self.feature_plots) == 0:
            return
        self.doc.new_line()
        self.doc.new_header(level=2, title='Feature Plots')
        self.doc.new_paragraph('Only contigs with features are shown')
        for contig in self.feature_plots.index.tolist():
            image_png = self.feature_plots[contig]
            self.doc.new_line(self.doc.new_inline_image(text='Analysis', path=os.path.abspath(image_png)))

    def add_mutations(self):
        self.ofh.write("\nXXXXXX In add_mutations\n\n")
        # Make sure we looked for mutations
        if not getattr(self, 'did_call_amr_mutations', False):
            return
        mutations = self.amr_mutations
        self.doc.new_line()
        self.doc.new_header(level=2, title=self.mutation_title)
        for region_name in mutations.index.tolist():
            region_mutations = mutations[region_name].copy()
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
        method = '%s reads were mapped to the reference sequence using minimap2.' % self.read_type
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pandas.Series(method))
        method = 'Mutations were identified using samtools mpileup and varscan.'
        self.methods[self.mutation_methods_title] = self.methods[self.mutation_methods_title].append(pandas.Series(method))

    def add_amr_matrix(self):
        self.ofh.write("\nXXXXXX In add_amr_matrix\n\n")
        # Make sure that we have an AMR matrix to plot
        #if not getattr(self.analysis, 'did_draw_amr_matrix', False):
        #    return
        #amr_matrix_png = self.analysis.amr_matrix_png
        #self.doc.new_line()
        #self.doc.new_header(level=2, title=self.amr_matrix_title)
        #self.doc.new_line('AMR genes and mutations with their corresponding drugs.')
        #self.doc.new_line(
        #    self.doc.new_inline_image(
        #        text='AMR genes and mutations with their corresponding drugs',
        #        path=amr_matrix_png
        #    )
        #)
        pass

    def add_large_indels(self):
        self.ofh.write("\nXXXXXX In add_large_indels\n\n")
        # Make sure we looked for mutations
        #if len(self.analysis.large_indels) == 0:
        #    return
        #large_indels = self.analysis.large_indels
        #self.doc.new_line()
        #self.doc.new_header(level=2, title=self.large_indel_title)
        #for genome in ['Reference insertions', 'Query insertions']:
        #    genome_indels = large_indels[genome].copy()
        #    self.doc.new_line()
        #    self.doc.new_header(level=3, title=genome)
        #    if (genome_indels.shape[0] == 0):
        #        continue
        #    genome_indels.iloc[:, 1] = genome_indels.iloc[:, 1].apply(lambda x: '{:,}'.format(x))
        #    genome_indels.iloc[:, 2] = genome_indels.iloc[:, 2].apply(lambda x: '{:,}'.format(x))
        #    genome_indels.iloc[:, 3] = genome_indels.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
        #    Table_List = [
        #        'Reference contig', 'Start', 'Stop', 'Size (bp)'
        #    ]
        #    for i in range(genome_indels.shape[0]):
        #        Table_List = Table_List + genome_indels.iloc[i, :].values.tolist()
        #    row_count = int(len(Table_List) / 4)
        #    self.doc.new_table(columns=4, rows=row_count, text=Table_List, text_align='left')
        #method = 'Large insertions or deletions were found as the complement of aligned regions using bedtools.'
        #self.methods[self.reference_methods_title] = self.methods[self.reference_methods_title].append(
        #    pandas.Series(method))
        #self.doc.new_line()
        #self.doc.new_line('<div style="page-break-after: always;"></div>')
        #self.doc.new_line()
        pass

    def add_plasmids(self):
        self.ofh.write("\nXXXXXX In add_plasmids\n\n")
        """
        if not getattr(self.analysis, 'did_call_plasmids', False):
            return
        # Make sure we looked for mutations
        #plasmids = self.analysis.plasmids
        if plasmids is None:
            return
        plasmids = plasmids.copy()
        self.doc.new_line()
        #self.doc.new_header(level=2, title=self.analysis.plasmid_title)
        if (plasmids.shape[0] == 0):
            self.doc.new_line('None')
            return
        plasmids.iloc[:, 3] = plasmids.iloc[:, 3].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 4] = plasmids.iloc[:, 4].apply(lambda x: '{:,}'.format(x))
        plasmids.iloc[:, 5] = plasmids.iloc[:, 5].apply(lambda x: '{:,}'.format(x))
        Table_List = [
            'Genome contig',
            'Plasmid hit',
            'Plasmid acc.',
            'Contig size',
            'Aliged',
            'Plasmid size'
        ]
        for i in range(plasmids.shape[0]):
            Table_List = Table_List + plasmids.iloc[i, 0:6].values.tolist()
        row_count = int(len(Table_List) / 6)
        self.doc.new_table(columns=6, rows=row_count, text=Table_List, text_align='left')
        method = 'The plasmid reference database was queried against the genome assembly using minimap2.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))
        method = 'The resulting SAM was converted to a PSL using a custom version of sam2psl.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))
        method = 'Plasmid-to-genome hits were resolved using the pChunks algorithm.'
        self.methods[self.plasmid_methods_title] = self.methods[self.plasmid_methods_title].append(pandas.Series(method))
        """
        pass

    def add_methods(self):
        self.ofh.write("\nXXXXXX In add_methods\n\n")
        self.doc.new_line('<div style="page-break-after: always;"></div>')
        self.doc.new_line()
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
        self.add_assembly_notes()
        if self.did_flye_ont_fastq:
            method = 'ONT reads were assembled using Flye.'
            self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pandas.Series(method))
        if self.did_medaka_ont_assembly:
            method = 'the genome assembly was polished using ont reads and medaka.'
            self.methods[self.assembly_methods_title] = self.methods[self.assembly_methods_title].append(pandas.series(method))

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

parser.add_argument('--analysis_name', action='store', dest='analysis_name', help='Sample identifier')
parser.add_argument('--assembly_fasta_file', action='store', dest='assembly_fasta_file', help='Assembly fasta file')
parser.add_argument('--assembly_name', action='store', dest='assembly_name', help='Assembly identifier')
parser.add_argument('--contig_coverage_file', action='store', dest='contig_coverage_file', help='Contig coverage TSV file')
parser.add_argument('--dbkey', action='store', dest='dbkey', help='Reference genome')
parser.add_argument('--illumina_fastq_file', action='store', dest='illumina_fastq_file', help='Input sample')
parser.add_argument('--gzipped', action='store_true', dest='gzipped', default=False, help='Input sample is gzipped')
parser.add_argument('--pima_css', action='store', dest='pima_css', help='PIMM css stypesheet')

args = parser.parse_args()

markdown_report = PimaReport(args.analysis_name, args.assembly_fasta_file, args.assembly_name, args.contig_coverage_file, args.dbkey, args.illumina_fastq_file, args.gzipped, args.pima_css)
markdown_report.make_report()
