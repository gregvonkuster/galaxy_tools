#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--bayesian"), action="store", dest="bayesian", help="Perform the estimations based on posterior estimations of methylation"),
    make_option(c("--col1"), action="store", dest="col1", type="integer", default=NULL, help="Column number of methylated read counts"),
    make_option(c("--col2"), action="store", dest="col2", type="integer", default=NULL, help="Column number of unmethylated read counts"),
    make_option(c("--high_coverage"), action="store", dest="high_coverage", type="integer", default=NULL, help="High coverage floor for samples"),
    make_option(c("--input_indiv_dir"), action="store", dest="input_indiv_dir", help="Directory of GRanges files"),
    make_option(c("--jd"), action="store", dest="jd", default=NULL, help="Add a column with values of J-information divergence"),
    make_option(c("--logbase"), action="store", dest="logbase", default=NULL, help="Logarithm base used to compute the J-information divergence"),
    make_option(c("--meth_level"), action="store", dest="meth_level", default=NULL, help="Perform the estimations based on posterior estimations of methylation levels"),
    make_option(c("--mcov1"), action="store", dest="mcov1", type="integer", help="Minimum coverage for cytosine site in sample 1"),
    make_option(c("--mcov2"), action="store", dest="mcov2", type="integer", default=NULL, help="Minimum coverage for cytosine site in sample 2"),
    make_option(c("--mmeth1"), action="store", dest="mmeth1", type="integer", help="Minimum read counts of methylated cytosine in sample 1"),
    make_option(c("--mmeth2"), action="store", dest="mmeth2", type="integer", default=NULL, help="Minimum read counts of methylated cytosine in sample 2"),
    make_option(c("--min_sitecov"), action="store", dest="min_sitecov", type="integer", help="Minimum total coverage for samples"),
    make_option(c("--mum1"), action="store", dest="mum1", type="integer", help="Minimum read counts of unmethylated cytosine in sample 1"),
    make_option(c("--mum2"), action="store", dest="mum2", type="integer", default=NULL, help="Minimum read counts of unmethylated cytosine in sample 2"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="The number of cores to use"),
    make_option(c("--output"), action="store", dest="output", help="Output Granges converted to data frame file"),
    make_option(c("--percentile"), action="store", dest="percentile", type="double", help="Threshold to remove the outliers from each file and all files stacked"),
    make_option(c("--ref"), action="store", dest="ref", help="A GRange file for the reference individual")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

cat("\nargs:\n", toString(args), "\n");

# Convert bayesian to boolean.
if (is.null(opt$bayesian)) {
    bayesian <- FALSE;
} else if (opt$bayesian=='yes') {
    bayesian <- TRUE;
} else {
    bayesian <- FALSE;
}

# Convert columns to an integer or integer vector if not NULL.
if (!is.null(opt$col1) & !is.null(opt$col2)) {
    columns <- integer(length=2);
    columns[[1]] <- opt$col1;
    columns[[2]] <- opt$col2;
} else if (!is.null(opt$col1)) {
    columns <- opt$col1;
} else {
    columns <- NULL;
}

# Read the reference individual into a data frame.
ref = readRDS(file=opt$ref);

# Read the list of input files for the individual into data frames.
input_indiv_files <- list.files(path=opt$input_indiv_dir, full.names=TRUE);
cat("\ninput_indiv_files: ", toString(input_indiv_files), "\n");
num_input_indiv_files <- length(input_indiv_files);
cat("num_input_indiv_files: ", num_input_indiv_files, "\n");

grange_list <- list();
for (i in 1:num_input_indiv_files) {
    cat("i: ", i, "\n");
    input_indiv_file <- normalizePath(input_indiv_files[[i]]);
    grange_list[[i]] <- readRDS(input_indiv_file);
}

# Convert jd to a boolean.
if (is.null(opt$jd)) {
    jd <- FALSE;
} else if (opt$jd=='yes') {
    jd <- TRUE;
} else {
    jd <- FALSE;
}

# Convert meth_level to a boolean.
if (is.null(opt$meth_level)) {
    meth_level <- FALSE;
} else if (opt$meth_level=='yes') {
    meth_level <- TRUE;
} else {
    meth_level <- FALSE;
}

# Convert min_coverage to an integer or integer vector if not NULL.
if (!is.null(opt$mcov1) && !is.null(opt$mcov2)) {
    min_coverage <- integer(length=2);
    min_coverage[[1]] <- opt$mcov1;
    min_coverage[[2]] <- opt$mcov2;
} else if (!is.null(opt$mcov1)) {
    min_coverage <- opt$mcov1;
}

cat("\nmin_coverage: ", min_coverage, "\n");

# Convert min_meth to an integer or integer vector if not NULL.
if (!is.null(opt$mmeth1) && !is.null(opt$mmeth2)) {
    min_meth <- integer(length=2);
    min_meth[[1]] <- opt$mmeth1;
    min_meth[[2]] <- opt$mmeth2;
} else if (!is.null(opt$mmeth1)) {
    min_meth <- opt$mmeth1;
}

cat("\nmin_meth: ", min_meth, "\n");

# Convert min_umeth to an integer or integer vector if not NULL.
if (!is.null(opt$mum1) && !is.null(opt$mum2)) {
    min_umeth <- integer(length=2);
    min_umeth[[1]] <- opt$mum1;
    min_umeth[[2]] <- opt$mum2;
} else if (!is.null(opt$mum1)) {
    min_umeth <- opt$mum1;
}

percentile <- formatC(opt$percentile, digits=3, format="f")


cat("\nminu_meth: ", min_umeth, "\n");
cat("\nbayesian: ", bayesian, "\n");
cat("\ncolumns: ", toString(columns), "\n");
cat("\nmin_coverage: ", min_coverage, "\n");
cat("\nmin_meth: ", min_meth, "\n");
cat("\nmin_umeth: ", min_umeth, "\n");
cat("\nopt$min_sitecov: ", opt$min_sitecov, "\n");
cat("\nopt$high_coverage: ", opt$high_coverage, "\n");
cat("\npercentile: ", percentile, "\n");
cat("\njd: ", jd, "\n");
cat("\nopt$num_cores: ", opt$num_cores, "\n");
cat("\nmeth_level: ", meth_level, "\n");
cat("\nopt$logbase: ", opt$logbase, "\n");

# Compute the information divergences of methylation levels.
inf_div_obj <- estimateDivergence(ref,
                                  grange_list,
                                  bayesian=bayesian,
                                  columns=columns,
                                  min.coverage=min_coverage,
                                  min.meth=min_meth,
                                  min.umeth=min_umeth,
                                  min.sitecov=opt$min_sitecov,
                                  high.coverage=opt$high_coverage,
                                  percentile=percentile,
                                  jd=jd,
                                  num.cores=opt$num_cores,
                                  tasks=0L,
                                  meth.level=meth_level,
                                  logbase=opt$logbase,
                                  verbose=TRUE);

cat("inf_div_obj:\n", toString(inf_div_obj, "\n");
cat("typeof(inf_div_obj:\n");
typeof(inf_div_obj

