#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("xtable"))

get_empirical_cumulative_probability_distributions_critical_values <- function(grange_list) {
    # FIXME: this function currently throws this exception:
    # error in evaluating the argument 'args' in selecting a method for function 'do.call': error in evaluating the
    # argument 'x' in selecting a method for function 'quantile': non-numeric argument to mathematical function
    # Calls: ... quantile -> .handleSimpleError -> h -> .handleSimpleError -> h
    critical_val <- do.call(rbind, lapply(grange_list, function(x) {
        hd.95 = quantile(x$hdiv, 0.95)
        tv.95 = quantile(abs(x$bay.TV), 0.95)
        return(c(tv = tv.95, hd = hd.95))
    }))
    return(critical_val);
}

get_cytosine_site_coverage <- function(grange_list, high_coverage) {
    # Output cytosine read counts for grange_list.
    covr <- lapply(grange_list, function(x) {
        cov1 <- x$c1 + x$t1;
        cov2 <- x$c2 + x$t2;
        cov <- apply(cbind(cov1, cov2), 1, max);
        return(cov);
    })

    do.call(rbind, lapply(covr, function(x) {
        q60 <- quantile(x, 0.6);
        q9999 <- quantile(x, 0.9999);
        idx1 <- which(x >= q60);
        idx2 <- which(x <= high_coverage);
        q95 <- quantile(x, 0.95);
        idx <- intersect(idx1, idx2);
        return(c(round(summary(x)),
                 q60,
                 quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
                 '#sites_ge_8' = sum(x >= 8),
                 'q60_le_high_coverage' = sum((x >= q60) & (x <= high_coverage)),
                 '#sites_gt_high_coverage' = sum(x > high_coverage)
                )
              )
         }
    ))
}

get_methylated_read_statistics <- function(grange_list, high_coverage) {
    # Descriptive statistics for methylated reads.
    methc <- lapply(grange_list, function(x) {
        r <- apply(cbind(x$c1, x$c2), 1, max)
        return(r)
    })

    do.call(rbind, lapply(methc, function(x) {
        q10 <- quantile(x, 0.1)
        q60 <- quantile(x, 0.6)
        q9999 <- quantile(x, 0.9999)
        idx1 <- which(x >= q60)
        idx2 <- which(x <= high_coverage)
        q95 <- quantile(x, 0.95)
        idx <- intersect(idx1, idx2)
        return(c(round(summary(x)),
                q10,
                q60,
                quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
                '#sites_ge_8' = sum(x >= 8),
                'q60_to_high_coverage' = sum((x >= q60) & (x <= high_coverage)),
                '#sites_gt_high_coverage' = sum(x > high_coverage)))
        })
    )
}

output_statistics <- function(output_file, csc_df, mrs_df) {
    sink(output_file);
    cat("<html><head></head><body>");
    cat("<h3>Cytosine Site Coverage</h3>");
    print(xtable(csc_df), type="html");
    cat("<h3>Methylated Reads Statistics</h3>");
    print(xtable(mrs_df), type="html");
    # cat("<h3>Critical Values from Empirical Cumulative Probability Distributions</h3>");
    # print(xtable(critical_vals_df), type="html");
    cat("</body></html>");
    sink();
}

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
    make_option(c("--output_crc"), action="store", dest="output_crc", help="Output cytosince read counts file"),
    make_option(c("--output_infdiv"), action="store", dest="output_infdiv", help="Output infDiv file"),
    make_option(c("--percentile"), action="store", dest="percentile", type="double", help="Threshold to remove the outliers from each file and all files stacked"),
    make_option(c("--ref"), action="store", dest="ref", help="A GRange file for the reference individual")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

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
num_input_indiv_files <- length(input_indiv_files);

grange_list <- list();
for (i in 1:num_input_indiv_files) {
    input_indiv_file <- input_indiv_files[[i]];
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

# Convert min_meth to an integer or integer vector if not NULL.
if (!is.null(opt$mmeth1) && !is.null(opt$mmeth2)) {
    min_meth <- integer(length=2);
    min_meth[[1]] <- opt$mmeth1;
    min_meth[[2]] <- opt$mmeth2;
} else if (!is.null(opt$mmeth1)) {
    min_meth <- opt$mmeth1;
}

# Convert min_umeth to an integer or integer vector if not NULL.
if (!is.null(opt$mum1) && !is.null(opt$mum2)) {
    min_umeth <- integer(length=2);
    min_umeth[[1]] <- opt$mum1;
    min_umeth[[2]] <- opt$mum2;
} else if (!is.null(opt$mum1)) {
    min_umeth <- opt$mum1;
}

percentile <- formatC(opt$percentile, digits=3, format="f")

############
# Debugging.
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
cat("ref: \n");
show(ref);
cat("\n\n");
cat("grange_list: \n");
for (i in 1:num_input_indiv_files) {
    grange <- grange_list[[i]];
    show(grange);
    cat("\n\n");
}
############

# Compute the information divergences of methylation levels.
infDiv <- estimateDivergence(ref,
                             grange_list,
                             Bayesian=bayesian,
                             columns=columns,
                             min.coverage=min_coverage,
                             min.meth=min_meth,
                             min.umeth=min_umeth,
                             min.sitecov=opt$min_sitecov,
                             high.coverage=opt$high_coverage,
                             percentile=percentile,
                             JD=jd,
                             num.cores=opt$num_cores,
                             meth.level=meth_level,
                             logbase=opt$logbase,
                             verbose=TRUE);

saveRDS(infDiv, file=opt$output_infdiv, compress=TRUE);

csc_df <- get_cytosine_site_coverage(infDiv, opt$high_coverage);

############
# Debugging.
cat("\ninfDiv names: ", toString(names(infDiv)), "\n");
cat("\ncytosine site coverage:\n");
show(csc_df);
cat("\n\n");
############

mrs_df <- get_methylated_read_statistics(infDiv, opt$high_coverage);

############
# Debugging.
cat("\nmethylated reads statistics:\n");
show(mrs_df);
cat("\n\n");
############

# FIXME: see the function above.
# critical_vals_df <- get_empirical_cumulative_probability_distributions_critical_values(grange_list);

############
# Debugging.
# cat("\ncritical values from empirical cumulative probability distributions:\n");
# show(critical_vals_df);
# cat("\n\n");
############

# output_statistics(opt$output_crc, csc_df, mrs_df, critical_vals_df);
output_statistics(opt$output_crc, csc_df, mrs_df);

