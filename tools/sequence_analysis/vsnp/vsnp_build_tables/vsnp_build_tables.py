#!/usr/bin/env python

import argparse
import os
import pandas
import pandas.io.formats.excel
import re

INPUT_JSON_DIR = 'input_json_dir'
INPUT_NEWICK_DIR = 'input_newick_dir'
# Maximum columns allowed in a LibreOffice
# spreadsheet is 1024.  Excel allows for
# 16,384 columns, but we'll set the lower
# number as the maximum since Galaxy is
# mostly run on Linux.
MAXCOLS = 1023
OUTPUT_EXCEL_DIR = 'output_excel_dir'


def excel_formatter(json_df, excel_file_name, group, gbk=None):
    pandas.io.formats.excel.header_style = None
    table_df = pandas.read_json(json_df, orient='split')
    # TODO: add support for a genbank file.
    if gbk:
        # table_df = annotate_table(table_df, group, gbk)
        pass
    else:
        table_df = table_df.append(pandas.Series(name='no annotations'))
    writer = pandas.ExcelWriter(excel_file_name, engine='xlsxwriter')
    table_df.to_excel(writer, sheet_name='Sheet1')
    writer_book = writer.book
    ws = writer.sheets['Sheet1']
    format_a = writer_book.add_format({'bg_color': '#58FA82'})
    format_g = writer_book.add_format({'bg_color': '#F7FE2E'})
    format_c = writer_book.add_format({'bg_color': '#0000FF'})
    format_t = writer_book.add_format({'bg_color': '#FF0000'})
    format_normal = writer_book.add_format({'bg_color': '#FDFEFE'})
    formatlowqual = writer_book.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    format_ambigous = writer_book.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    format_n = writer_book.add_format({'bg_color': '#E2CFDD'})
    rows, cols = table_df.shape
    ws.set_column(0, 0, 30)
    ws.set_column(1, cols, 2.1)
    ws.freeze_panes(2, 1)
    format_annotation = writer_book.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
    # Set last row.
    ws.set_row(rows+1, cols+1, format_annotation)
    # Make sure that row/column locations don't overlap.
    ws.conditional_format(rows - 2, 1, rows - 1, cols, {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': format_normal})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'A', 'format': format_a})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'G', 'format': format_g})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'C', 'format': format_c})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'T', 'format': format_t})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'S', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'Y', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'R', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'W', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'K', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'M', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'N', 'format': format_n})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': '-', 'format': format_n})
    format_rotation = writer_book.add_format({})
    format_rotation.set_rotation(90)
    for columnnum, columnname in enumerate(list(table_df.columns)):
        ws.write(0, columnnum + 1, columnname, format_rotation)
    format_annotation = writer_book.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
    # Set last row.
    ws.set_row(rows, 400, format_annotation)
    writer.save()


def output_cascade_table(cascade_order, mqdf, group):
    cascade_order_mq = pandas.concat([cascade_order, mqdf], join='inner')
    output_table(cascade_order_mq, "cascade", group)


def output_excel(df, type_str, group, count=None):
    # Output the temporary json file that
    # is used by the excel_formatter.
    if count is None:
        json_file_name = "%s_%s_order_mq.json" % (group, type_str)
    else:
        json_file_name = "%s_%s_order_mq_%d.json" % (group, type_str, count)
    df.to_json(json_file_name, orient='split')
    # Grouping took place upstream, so we
    # have inputs that are collections.
    if count is None:
        excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_table.xlsx" % (group, type_str))
    else:
        excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_table_%d.xlsx" % (group, type_str, count))
    # Output the Excel file.
    excel_formatter(json_file_name, excel_file_name, group, gbk=None)


def output_sort_table(cascade_order, mqdf, group):
    sort_df = cascade_order.T
    sort_df['abs_value'] = sort_df.index
    sort_df[['chrom', 'pos']] = sort_df['abs_value'].str.split(':', expand=True)
    sort_df = sort_df.drop(['abs_value', 'chrom'], axis=1)
    sort_df.pos = sort_df.pos.astype(int)
    sort_df = sort_df.sort_values(by=['pos'])
    sort_df = sort_df.drop(['pos'], axis=1)
    sort_df = sort_df.T
    sort_order_mq = pandas.concat([sort_df, mqdf], join='inner')
    output_table(sort_order_mq, "sort", group)


def output_table(df, type_str, group):
    count = 0
    chunk_start = 0
    chunk_end = 0
    column_count = df.shape[1]
    if column_count >= MAXCOLS:
        # Here the number of columns is greater than
        # the maximum allowed by Excel, so multiple
        # outputs will be produced.
        while column_count >= MAXCOLS:
            count += 1
            chunk_end += MAXCOLS
            df_of_type = df.iloc[:, chunk_start:chunk_end]
            output_excel(df_of_type, type_str, group, count=count)
            chunk_start += MAXCOLS
            column_count -= MAXCOLS
        count += 1
        df_of_type = df.iloc[:, chunk_start:]
        output_excel(df_of_type, type_str, group, count=count)
    else:
        output_excel(df, type_str, group)


parser = argparse.ArgumentParser()

parser.add_argument('--input_avg_mq_json', action='store', dest='input_avg_mq_json', help='Average MQ json file')
parser.add_argument('--input_newick', action='store', dest='input_newick', required=False, default=None, help='Newick file')
parser.add_argument('--input_snps_json', action='store', dest='input_snps_json', required=False, default=None, help='SNPs json file')

args = parser.parse_args()

avg_mq_series = pandas.read_json(args.input_avg_mq_json, typ='series', orient='split')
# Map quality to dataframe.
mqdf = avg_mq_series.to_frame(name='MQ')
mqdf = mqdf.T
# The assumption here is that the list of files
# in both INPUT_NEWICK_DIR and INPUT_JSON_DIR are
# named such that they are properly matched if
# the directories contain more than 1 file.
if args.input_newick is not None:
    newick_files = [args.input_newick]
else:
    newick_files = []
    for file_name in os.listdir(INPUT_NEWICK_DIR):
        file_path = os.path.abspath(os.path.join(INPUT_NEWICK_DIR, file_name))
        newick_files.append(file_path)
if args.input_snps_json is not None:
    json_files = [args.input_snps_json]
else:
    json_files = []
    for file_name in os.listdir(INPUT_JSON_DIR):
        file_path = os.path.abspath(os.path.join(INPUT_JSON_DIR, file_name))
        json_files.append(file_path)

for i, newick_file in enumerate(newick_files):
    json_file = json_files[i]
    # Hopefully the newick file name will be something like
    # Mbovis-01D6_snps.fasta so the group is meaningful.
    newick_file_name_base = os.path.basename(newick_file)
    # Eliminate the extension.
    newick_file_name_base = os.path.splitext(newick_file_name_base)[0]
    # Get the group.
    group = newick_file_name_base.split("_")[0]
    snps_df = pandas.read_json(json_file, orient='split')
    with open(newick_file, 'r') as fh:
        for line in fh:
            line = re.sub('[:,]', '\n', line)
            line = re.sub('[)(]', '', line)
            line = re.sub('[0-9].*\.[0-9].*\n', '', line)
            line = re.sub('root\n', '', line)
    sample_order = line.split('\n')
    sample_order = list(filter(None, sample_order))
    sample_order.insert(0, 'root')
    tree_order = snps_df.loc[sample_order]
    # Count number of SNPs in each column.
    snp_per_column = []
    for column_header in tree_order:
        count = 0
        column = tree_order[column_header]
        for element in column:
            if element != column[0]:
                count = count + 1
        snp_per_column.append(count)
    row1 = pandas.Series(snp_per_column, tree_order.columns, name="snp_per_column")
    # Count number of SNPS from the
    # top of each column in the table.
    snp_from_top = []
    for column_header in tree_order:
        count = 0
        column = tree_order[column_header]
        # for each element in the column
        # skip the first element
        for element in column[1:]:
            if element == column[0]:
                count = count + 1
            else:
                break
        snp_from_top.append(count)
    row2 = pandas.Series(snp_from_top, tree_order.columns, name="snp_from_top")
    tree_order = tree_order.append([row1])
    tree_order = tree_order.append([row2])
    # In pandas=0.18.1 even this does not work:
    # abc = row1.to_frame()
    # abc = abc.T --> tree_order.shape (5, 18), abc.shape (1, 18)
    # tree_order.append(abc)
    # Continue to get error: "*** ValueError: all the input arrays must have same number of dimensions"
    tree_order = tree_order.T
    tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
    tree_order = tree_order.T
    # Remove snp_per_column and snp_from_top rows.
    cascade_order = tree_order[:-2]
    # Output the cascade table.
    output_cascade_table(cascade_order, mqdf, group)
    # Output the sorted table.
    output_sort_table(cascade_order, mqdf, group)