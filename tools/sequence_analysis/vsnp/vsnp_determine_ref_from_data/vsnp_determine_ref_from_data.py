#!/usr/bin/env python

__version__ = "2.0.0"

import argparse
import gzip
import os
import re
import shutil
import sys
from collections import OrderedDict
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class GetReference:

    def __init__(self, read1, read2):
        self.fastq_list = []
        for read in [read1, read2]:
            if read is not None:
                self.fastq_list.append(read)
        if read2 is not None:
            self.paired = True
        else:
            self.paired = False
        self.sample_name = re.sub('[_.].*', '', os.path.basename(read1))
        # TODO: move these dicts to external yaml or tabular files.
        self.bovis_dict = self.get_bovis_dict()
        self.brucella_dict = self.get_brucella_dict()
        self.oligo_dict = self.get_oligo_dict()
        self.para_dict = self.get_para_dict()

    def get_bovis_dict(self):
        bovis = {}
        # Mycobacterium tuberculosis H37Rv tb1
        bovis["11101111"] = "NC_000962"
        bovis["11101101"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb2
        bovis["01100111"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb3
        bovis["01101011"] = "NC_000962"
        bovis["11101011"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb4a
        bovis["01101111"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb4b
        bovis["01101101"] = "NC_000962"
        bovis["11101101"] = "NC_000962"
        bovis["01101111"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb5
        bovis["11111111"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb6
        bovis["11001111"] = "NC_000962"
        # Mycobacterium tuberculosis H37Rv tb7
        bovis["10101110"] = "NC_000962"
        # Mycobacterium bovis AF2122/97
        bovis["11001110"] = "AF2122"
        bovis["11011110"] = "AF2122"
        bovis["11001100"] = "AF2122"
        return bovis

    def get_brucella_dict(self):
        brucella = {}
        # Unexpected findings
        brucella["1111111111111111"] = "odd"
        # Brucella abortus bv. 1 str. 9-941
        brucella["0111111111111111"] = "NC_006932"
        brucella["1101111111111111"] = "NC_006932"
        # Brucella abortus strain BER
        brucella["1011111111111111"] = "NZ_CP007682"
        # Brucella melitensis bv. 1 str. 16M
        brucella["1110111111111101"] = "NC_003317"
        brucella["0000010101101101"] = "NC_003317"
        # Brucella melitensis BwIM_SOM_36b
        brucella["1110111111111100"] = "NZ_CP018508"
        brucella["0000010101101100"] = "NZ_CP018508"
        # Brucella melitensis ATCC 23457
        brucella["1110111111111011"] = "NC_012441"
        brucella["0000010101101001"] = "NC_012441"
        brucella["0100010101101001"] = "NC_012441"
        brucella["1110011111101011"] = "NC_012441"
        # Brucella melitensis bv. 3 str. Ether
        brucella["1110111111110111"] = "NZ_CP007760"
        brucella["1110011111100111"] = "NZ_CP007760"
        # Brucella suis 1330
        brucella["1111011111111111"] = "NC_017251"
        # Brucella suis ATCC 23445
        brucella["1111101111111111"] = "NC_010169"
        # Brucella suis bv. 3 str. 686
        brucella["1111110111111101"] = "NZ_CP007719"
        # FIXME: what is the NCBI ID?
        brucella["1111111011111111"] = "Brucella_ceti1"
        brucella["1111111001111111"] = "Brucella_ceti1"
        # Brucella ceti TE10759-12
        brucella["1111111101111111"] = "NC_022905"
        # FIXME: what is the NCBI ID?
        brucella["1111111110111101"] = "Brucella_suis4"
        # Brucella canis ATCC 23365
        brucella["1111111110011101"] = "NC_010103"
        # Brucella ovis ATCC 25840
        brucella["1111111111101111"] = "NC_009505"
        return brucella

    def get_oligo_dict(self):
        oligo = {}
        oligo["01_ab1"] = "AATTGTCGGATAGCCTGGCGATAACGACGC"
        oligo["02_ab3"] = "CACACGCGGGCCGGAACTGCCGCAAATGAC"
        oligo["03_ab5"] = "GCTGAAGCGGCAGACCGGCAGAACGAATAT"
        oligo["04_mel"] = "TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
        oligo["05_suis1"] = "TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
        oligo["06_suis2"] = "GGCAATCATGCGCAGGGCTTTGCATTCGTC"
        oligo["07_suis3"] = "CAAGGCAGATGCACATAATCCGGCGACCCG"
        oligo["08_ceti1"] = "GTGAATATAGGGTGAATTGATCTTCAGCCG"
        oligo["09_ceti2"] = "TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
        oligo["10_canis4"] = "CTGCTACATAAAGCACCCGGCGACCGAGTT"
        oligo["11_canis"] = "ATCGTTTTGCGGCATATCGCTGACCACAGC"
        oligo["12_ovis"] = "CACTCAATCTTCTCTACGGGCGTGGTATCC"
        oligo["13_ether2"] = "CGAAATCGTGGTGAAGGACGGGACCGAACC"
        oligo["14_63B1"] = "CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
        oligo["15_16M0"] = "TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
        oligo["16_mel1b"] = "TCTGTCCAAACCCCGTGACCGAACAATAGA"
        oligo["17_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
        oligo["18_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
        oligo["19_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
        oligo["20_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
        oligo["21_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
        oligo["22_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
        oligo["23_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
        oligo["24_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
        oligo["25_para"] = "CCTTTCTTGAAGGGTGTTCG"
        oligo["26_para_sheep"] = "CGTGGTGGCGACGGCGGCGGGCCTGTCTAT"
        oligo["27_para_cattle"] = "TCTCCTCGGTCGGTGATTCGGGGGCGCGGT"
        return oligo

    def get_para_dict(self):
        para = {}
        # Mycobacterium avium subsp. paratuberculosis strain Telford
        para["110"] = "CP033688"
        # Mycobacterium avium subsp. paratuberculosis K-10
        para["101"] = "NC_002944"
        return para

    def output_ref_and_metrics(self, tool_data_path, output_ref, output_metrics, gzipped):
        count_summary = {}

        for v1 in self.oligo_dict.values():
            returned_value, count = self.get_seq_counts(v1, gzipped)
            for key, v2 in self.oligo_dict.items():
                if returned_value == v2:
                    count_summary.update({key: count})
                    # FIXME: can this be moved out of the loop?
                    count_summary = OrderedDict(sorted(count_summary.items()))

        count_list = []
        for v in count_summary.values():
            count_list.append(v)
        brucella_sum = sum(count_list[:16])
        bovis_sum = sum(count_list[16:24])
        para_sum = sum(count_list[24:])

        # Binary dictionary
        binary_dictionary = {}
        for k, v in count_summary.items():
            if v > 1:
                binary_dictionary.update({k: 1})
            else:
                binary_dictionary.update({k: 0})
        binary_dictionary = OrderedDict(sorted(binary_dictionary.items()))

        binary_list = []
        for v in binary_dictionary.values():
            binary_list.append(v)

        brucella_binary = binary_list[:16]
        brucella_string = ''.join(str(e) for e in brucella_binary)

        bovis_binary = binary_list[16:24]
        bovis_string = ''.join(str(e) for e in bovis_binary)

        para_binary = binary_list[24:]
        para_string = ''.join(str(e) for e in para_binary)

        with open(output_metrics, "w") as fh:
            fh.write("Brucella counts: ")
            for i in count_list[:16]:
                fh.write("%d," % i)
            fh.write("\nTB counts: ")
            for i in count_list[16:24]:
                fh.write("%d," % i)
            fh.write("\nPara counts: ")
            for i in count_list[24:]:
                fh.write("%d," % i)

            if brucella_sum > 3:
                group = "Brucella"
                if brucella_string in self.brucella_dict:
                    dbkey = self.brucella_dict[brucella_string]
                else:
                    dbkey = ""
            elif bovis_sum > 3:
                group = "TB"
                if bovis_string in self.bovis_dict:
                    dbkey = self.bovis_dict[bovis_string]
                else:
                    dbkey = ""
            elif para_sum >= 1:
                group = "paraTB"
                if para_string in self.para_dict:
                    dbkey = self.para_dict[para_string]
                else:
                    dbkey = ""
            else:
                group = ""
                dbkey = ""
            fh.write("\nGroup: %s, dbkey: %s\n" % (group, dbkey))
        reference_path = os.path.join(tool_data_path, dbkey, 'seq', '%s.fa' % dbkey)
        if os.path.isfile(reference_path):
            shutil.copy(reference_path, output_ref)
        else:
            sys.exit("Reference not determined")

    def get_seq_counts(self, value, gzipped):
        count = 0
        for fastq in self.fastq_list:
            if gzipped:
                with gzip.open(fastq, 'rt') as fh:
                    for title, seq, qual in FastqGeneralIterator(fh):
                        count += seq.count(value)
            else:
                with open(fastq, 'rt') as fh:
                    for title, seq, qual in FastqGeneralIterator(fh):
                        count += seq.count(value)
        return(value, count)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-gz', '--gzipped', action='store', dest='gzipped', required=False, default=None, help='Input files are gzipped')
    parser.add_argument('-ol', '--output_ref', action='store', dest='output_ref', help='Output reference file')
    parser.add_argument('-om', '--output_metrics', action='store', dest='output_metrics', help='Output metrics file')
    parser.add_argument('-td', '--tool_data_path', action='store', dest='tool_data_path', help='Location of cached genomes')

    args = parser.parse_args()

    gzipped = args.gzipped is not None
    reference = GetReference(args.read1, args.read2)
    reference.output_ref_and_metrics(args.tool_data_path, args.output_ref, args.output_metrics, gzipped)
