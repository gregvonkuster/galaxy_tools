(.venv) [galaxy@IDEAS ideas_preprocessor]$ cat ideas_preprocessor.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
        make_option(c("--bed_input"), action="store", dest="bed_input", defaul=NULL, help="Chromosome windows positions file"),
        make_option(c("--exclude_input"), action="store", dest="exclude_input", defaul=NULL, help="File(s) containing regions to exclude"),
        make_option(c("--bychr"), action="store_true", dest="bychr", defaul=FALSE, help="Separate files by chromosome"),
        make_option(c("--chrom_len_file"), action="store", dest="chrom_len_file", default=NULL, help="Chromosome lengths file"),
        make_option(c("--ideas_input_config"), action="store", dest="ideas_input_config", help="Preprocessing output config file"),
        make_option(c("--ideaspre_input_config"), action="store", dest="ideaspre_input_config", help="Preprocessing input config file"),
        make_option(c("--output"), action="store", dest="output", help="Primary output dataset"),
        make_option(c("--output_files_path"), action="store", dest="output_files_path", help="Primary output dataset extra files path"),
        make_option(c("--output_hid"), action="store", dest="output_hid", help="Primary output dataset hid"),
        make_option(c("--reads_per_bp"), action="store", dest="reads_per_bp", type="integer", help="Calculate the signal in each genomic window"),
        make_option(c("--restrict_to_chroms"), action="store", dest="restrict_to_chroms", default=NULL, help="Restrict processing to specified chromosomes"),
        make_option(c("--standardize_datasets"), action="store_true", dest="standardize_datasets", default=FALSE, help="Standardize all datasets"),
        make_option(c("--chromosome_windows"), action="store", dest="chromosome_windows", default=NULL, help="Windows positions by chroms config file"),
        make_option(c("--window_size"), action="store", dest="window_size", type="integer", default=NULL, help="Window size in base pairs")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

create_primary_html = function(output, output_hid, output_files_path) {
    files = list.files(path=output_files_path);
    s <- paste('<html><head></head><body>', sep="\n");
    s <- paste(s, '<h3>History item ', output_hid, ' files prepared for IDEAS</h3>\n', sep="");
    s <- paste(s, '<ul>\n', sep="");
    for (i in 1:length(files)) {
        s <- paste(s, '<li><a href="', files[i], '">', files[i], '</a></li>\n', sep="");
    }
    s <- paste(s, '</ul>\n</body>\n</html>', sep="");
    cat(s, file=output);
}

tmp_dir = "tmp";
output_tmp_dir = paste(opt$output_files_path, tmp_dir, sep="/");
dir.create(output_tmp_dir, showWarnings=FALSE);

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
        cmd = paste("bamCoverage --bam", file_path, "-o", bigwig_file_name, "--binSize", opt$window_size);
        system(cmd);
    } else {
        bigwig_file_name = file_path;
    }
    bed_file_name = paste(file_name, "bed", sep=".");
    bed_file_path = paste("tmp", bed_file_name, sep="/");
    cmd = paste("bigWigAverageOverBed", bigwig_file_name, opt$bed_input, "stdout | cut -f5 >", bed_file_path);
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
# Create the ideas_input_config with the format required by IDEAS.
cmd = paste("paste -d' ' file1.txt file2.txt >", opt$ideas_input_config, sep=" ");
system(cmd);
# Move the ideas_input_config to the output directory.
to_path = paste(opt$output_files_path, opt$ideas_input_config, sep="/");
file.rename(opt$ideas_input_config, to_path);
# Move the compressed bed files in the tmp
# directory to the output tmp directory.
tmp_files = list.files(path=tmp_dir);
for (i in 1:length(tmp_files)) {
    from_path = paste(tmp_dir, tmp_files[i], sep="/");
    to_path = paste(output_tmp_dir, tmp_files[i], sep="/");
    file.rename(from_path, to_path);
}
if (!is.null(opt$chromosome_windows)) {
    # Move the chromosome_windows to the output directory.
    to_path = paste(opt$output_files_path, opt$chromosome_windows, sep="/");
    file.rename(opt$chromosome_windows, to_path);
}
# Create the primary HTML dataset.
create_primary_html(opt$output, opt$output_hid, opt$output_files_path);
