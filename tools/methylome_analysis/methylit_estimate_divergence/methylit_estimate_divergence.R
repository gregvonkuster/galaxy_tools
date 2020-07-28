#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--bayesian"), action="store", dest="bayesian", help="Perform the estimations based on posterior estimations of methylation"),
    make_option(c("--column1"), action="store", dest="column1", type="integer" default=NULL, help="Column number of methylated read counts"),
    make_option(c("--column2"), action="store", dest="column2", type="integer" default=NULL, help="Column number of unmethylated read counts"),
    make_option(c("--high_coverage"), action="store", dest="high_coverage", type="integer" default=NULL, help="High coverage floor for samples"),
    make_option(c("--input_indiv_dir"), action="store", dest="input_indiv_dir", help="Directory of GRanges files"),
    make_option(c("--jd"), action="store", dest="jd", default=NULL, help="Add a column with values of J-information divergence"),
    make_option(c("--logbase"), action="store", dest="logbase", default=NULL, help="Logarithm base used to compute the J-information divergence"),
    make_option(c("--meth_level"), action="store", dest="meth_level", default=NULL, help="Perform the estimations based on posterior estimations of methylation levels"),
    make_option(c("--min_coverage1"), action="store", dest="min_coverage1", type="integer" help="Minimum coverage for cytosine site in sample 1"),
    make_option(c("--min_coverage2"), action="store", dest="min_coverage2", type="integer" default=NULL, help="Minimum coverage for cytosine site in sample 2"),
    make_option(c("--min_meth1"), action="store", dest="min_meth1", type="integer" help="Minimum read counts of methylated cytosine in sample 1"),
    make_option(c("--min_meth2"), action="store", dest="min_meth2", type="integer" default=NULL, help="Minimum read counts of methylated cytosine in sample 2"),
    make_option(c("--min_sitecov"), action="store", dest="min_sitecov", type="integer" help="Minimum total coverage for samples"),
    make_option(c("--min_umeth1"), action="store", dest="min_umeth1", type="integer" help="Minimum read counts of unmethylated cytosine in sample 1"),
    make_option(c("--min_umeth2"), action="store", dest="min_umeth2", type="integer" default=NULL, help="Minimum read counts of unmethylated cytosine in sample 2"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer" help="The number of cores to use"),
    make_option(c("--output"), action="store", dest="output", help="Output Granges converted to data frame file"),
    make_option(c("--percentile"), action="store", dest="percentile", type="float" help="Threshold to remove the outliers from each file and all files stacked"),
    make_option(c("--ref"), action="store", dest="ref", help="The reference individual")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Convert bayesian to boolean.
if (opt$bayesian == 'yes') {
    bayesian = TRUE;
} else {
    bayesian = FALSE;
}

# Convert columns to an integer or integer vector if not NULL.
if (!is.null(opt$column1) && !is.null(opt$column2)) {
    columns <- integer(opt$column1, opt$column2);
} else if (!is.null(opt$column1)) {
    columns <- opt$column1;
} else {
    columns <- NULL;
}

# Get the list of input GRanges files for the individual.
input_indiv_files <- list.files(path=opt$input_indiv_dir, full.names=TRUE);
num_input_indiv_files <- length(input_indiv_files);
# Convert the list to a list of GRanges files.
granges_list <- list();
for (i in 1:num_input_indiv_files) {
    input_indiv_file <- input_indiv_files[i];
    grange = load(input_indiv_file);
    granges_list[[i]] <- grange;
}

# Convert jd to a boolean.
if (opt$jd == 'yes') {
    jd = TRUE;
} else {
    jd = FALSE;
}

# Convert meth_level to a boolean.
if (opt$meth_level == 'yes') {
    meth_level = TRUE;
} else {
    meth_level = FALSE;
}

# Convert min_coverage to an integer or integer vector if not NULL.
if (!is.null(opt$min_coverage1) && !is.null(opt$min_coverage2)) {
    min_coverage <- integer(opt$min_coverage1, opt$min_coverage2);
} else if (!is.null(opt$min_coverage1)) {
    min_coverage <- opt$min_coverage1;
}

# Convert min_meth to an integer or integer vector if not NULL.
if (!is.null(opt$min_meth1) && !is.null(opt$min_meth2)) {
    min_meth <- integer(opt$min_meth1, opt$min_meth2);
} else if (!is.null(opt$min_meth1)) {
    min_meth <- opt$min_meth1;
}

# Convert min_umeth to an integer or integer vector if not NULL.
if (!is.null(opt$min_umeth1) && !is.null(opt$min_umeth2)) {
    min_umeth <- integer(opt$min_umeth1, opt$min_umeth2);
} else if (!is.null(opt$min_umeth1)) {
    min_umeth <- opt$min_umeth1;
}

ref <- load(opt$ref);

# Compute the information divergences of methylation levels.
inf_div_obj <- estimateDivergence(ref,
                                  indiv,
                                  bayesian=bayesian,
                                  columns=columns,
                                  min.coverage=min_coverage,
                                  min.meth=min_meth,
                                  min.umeth=min_umeth,
                                  min.sitecov=opt$min_sitecov,
                                  high.coverage=opt$high_coverage,
                                  percentile=opt$percentile,
                                  jd=$jd,
                                  num.cores=opt$num_cores,
                                  tasks=0L,
                                  meth.level=meth_level,
                                  logbase=opt$logbase,
                                  verbose=TRUE);

cat("inf_div_obj:\n");
inf_div_obj
cat("typeof(inf_div_obj:\n");
typeof(inf_div_obj

