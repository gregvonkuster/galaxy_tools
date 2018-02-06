#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
        make_option(c("--chrom_bed_input"), action="store", dest="chrom_bed_input", defaul=NULL, help="Chromosome windows positions file"),
        make_option(c("--exclude_bed_input"), action="store", dest="exclude_bed_input", defaul=NULL, help="File(s) containing regions to exclude"),
        make_option(c("--chrom_len_file"), action="store", dest="chrom_len_file", default=NULL, help="Chromosome lengths file"),
        make_option(c("--ideaspre_input_config"), action="store", dest="ideaspre_input_config", help="Preprocessing input config file"),
        make_option(c("--output"), action="store", dest="output", help="Primary output dataset"),
        make_option(c("--output_files_path"), action="store", dest="output_files_path", help="Primary output dataset extra files path"),
        make_option(c("--chromosome_windows"), action="store", dest="chromosome_windows", default=NULL, help="Windows positions by chroms config file"),
        make_option(c("--window_size"), action="store", dest="window_size", type="integer", default=NULL, help="Window size in base pairs")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

tmp_dir = "tmp";
cbi_file = "chrom_bed_input.bed";

if (is.null(opt$chrom_bed_input)) {
    # Create a chromosome windows positions file
    # using the received chromosome lengths file
    # and the window size.
    cmd = paste("bedtools makewindows -g", opt$chrom_len_file, "-w", opt$window_size, ">", cbi_file, sep=" ");
    system(cmd);
} else {
    if (!is.null(opt$exclude_bed_input)) {
        # Copy the received chrom_bed_input
        # since we will alter it.
        file.copy(opt$chrom_bed_input, cbi_file);
    } else {
        cbi_file = opt$chrom_bed_input;
    }
}
# Exclude regions if necessary.
if (!is.null(opt$exclude_bed_input)) {
    exclude_bed_inputs = as.character(opt$exclude_bed_input);
    exclude_bed_files = strsplit(exclude_bed_inputs, ",");
    tmp_file = paste("tmp", cbi_file, sep="_");
    for (exclude_bed_file in exclude_bed_files) {
        cmd = paste("bedtools subtract -a", cbi_file, "-b", exclude_bed_file, ">", tmp_file, sep=" ");
        system(cmd);
        cmd = paste("mv", tmp_file, cbi_file, sep=" ");
        system(cmd);
    }
}
# Read the chromosome windows positions file
# to get the smallest window size in the file
# (i.e., the minimum of column 3 - column 2.
cbi = fread(cbi_file);
min_window_size = min(cbi[,3]-cbi[,2]);
# Read the ideaspre_input_config text file which has this format:
# "cell type name" "epigenetic factor name" "file path" "file name" "datatype"
ideaspre_input_config = as.matrix(read.table(opt$ideaspre_input_config));
# Process data to windows mean.
for (i in 1:dim(ideaspre_input_config)[1]) {
    file_path = ideaspre_input_config[i, 3]
    file_name = ideaspre_input_config[i, 4]
    datatype = ideaspre_input_config[i, 5]
    if (datatype == "bam") {
        cmd = paste("samtools index", file_path);
        system(cmd);
        bigwig_file_name = paste(file_name, "bw", sep=".");
        cmd = paste("bamCoverage --bam", file_path, "-o", bigwig_file_name, "--binSize", min_window_size);
        system(cmd);
    } else {
        bigwig_file_name = file_path;
    }
    bed_file_name = paste(file_name, "bed", sep=".");
    bed_file_path = paste("tmp", bed_file_name, sep="/");
    cmd = paste("bigWigAverageOverBed", bigwig_file_name, opt$chrom_bed_input, "stdout | cut -f5 >", bed_file_path);
    system(cmd);
    cmd = paste("gzip -f", bed_file_path);
    system(cmd);
}
# Create file1.txt.
cmd = paste("cut -d' '", opt$ideaspre_input_config, "-f1,2 > file1.txt", sep=" ");
system(cmd);
# Compress the bed files in the tmp directory.
tmp_gzipped_files = paste(tmp_dir, "*.bed.gz", sep="/");
# Create file2.txt.
cmd = paste("ls", tmp_gzipped_files, "> file2.txt", sep=" ");
system(cmd);
# Create IDEAS_input_config.txt  with the format required by IDEAS.
ideas_input_config = "IDEAS_input_config.txt"
cmd = paste("paste -d' ' file1.txt file2.txt >", ideas_input_config, sep=" " );
system(cmd);
# Move IDEAS_input_config.txt to the output directory.
to_path = paste(opt$output_files_path, ideas_input_config, sep="/");
file.rename(ideas_input_config, to_path);
# Archive the tmp directory.
cmd = "tar -cvf tmp.tar tmp";
system(cmd);
# Compress the archive.
cmd = "gzip tmp.tar";
system(cmd);
# Move the tmp archive to the output directory.
to_path = paste(opt$output_files_path, "tmp.tar.gz", sep="/");
file.rename("tmp.tar.gz", to_path);
# Handle file names for display in the primary dataset if necessary.
to_path = paste(opt$output_files_path, "chromosomes.bed", sep="/");
if (is.null(opt$chrom_bed_input)) {
    # Move cbi_file to the output directory,
    # naming it chromosomes.bed.
    file.rename(cbi_file, to_path);
} else {
    # Copy opt$chrom_bed_input to the output
    # directory, naming it chromosomes.bed.
    file.copy(opt$chrom_bed_input, to_path);
}
if (!is.null(opt$chromosome_windows)) {
    # Move chromosome_windows.txt to the output directory.
    to_path = paste(opt$output_files_path, opt$chromosome_windows, sep="/");
    file.rename(opt$chromosome_windows, to_path);
}
