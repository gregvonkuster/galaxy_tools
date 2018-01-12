#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--input_dir"), action="store", dest="input_dir", help="IDEAS para files directory"),
    make_option(c("--output_dir"), action="store", dest="output_dir", help="PDF output directory"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--in_training_mode"), action="store_true", dest="in_training_mode", default=FALSE, help="Flag for training mode")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

heatmap_path = paste(opt$script_dir, "create_heatmap.R", sep="/");
source(heatmap_path);

if (opt$in_training_mode) {
    ext = ".para0";
    pattern = "\\.para0$";
} else {
    ext = ".para";
    pattern = "\\.para$"
}
para_files = list.files(path=opt$input_dir, pattern=pattern, full.names=TRUE);
for (i in 1:length(para_files)) {
    para_file = para_files[i];
    para_file_base_name = strsplit(para_file, split="/")[[1]][2];
    output_file_base_name = gsub(ext, "", para_file_base_name);
    output_file_name = paste(output_file_base_name, "state", i, "pdf", sep=".");
    output_file_path = paste(opt$output_dir, output_file_name, sep="/");
    data_frame = read.table(para_file, comment="!", header=T);
    create_heatmap(data_frame, output_file_path);
}
