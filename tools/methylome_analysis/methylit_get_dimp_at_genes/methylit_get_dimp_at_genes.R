#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocGenerics"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("S4Vectors"))
suppressPackageStartupMessages(library("xtable"))

option_list <- list(
    make_option(c("--gene_annot"), action="store", dest="gene_annot", help="Data Manager installed gene annotation file"),
    make_option(c("--ignore_strand"), action="store", dest="ignore_strand", help="Ignore strand"),
    make_option(c("--input"), action="store", dest="input", help="File containing any of the following objects: pDMP, InfDiv, GRangesList, GRanges or a list of GRanges"),
    make_option(c("--output_rdata"), action="store", dest="output_rdata", help="Output rdata file"),
    make_option(c("--output_type"), action="store", dest="output_type", help="Output type"),
    make_option(c("--overlap_type"), action="store", dest="overlap_type", help="Overlap type"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

ignore_strand <- string_to_boolean(opt$ignore_strand, default=FALSE);
gene_annot <- readGFFAsGRanges(opt$gene_annot);

############
# Debugging.
cat("\nclass(gene_annot):\n");
class(gene_annot);
cat("\n\n");
cat("\nhead(gene_annot):\n");
head(gene_annot);
cat("\n\n");
############

genes <- gene_annot[ gene_annot$type == "gene", c( "gene_id", "biotype", "Name" ) ];
genes <- genes[ genes$biotype == "protein_coding", "gene_id" ];
seqlevels(genes) <- "chr13";

############
# Debugging.
cat("\ngenes:\n");
show(genes);
cat("\n\n");
cat("\nclass(genes):\n");
class(genes);
cat("\n\n");
############

gr <- readRDS(opt$input);
############
# Debugging.
cat("\n\ngr:\n");
gr
cat("\n\n");
cat("\nclass(gr):\n");
class(gr);
cat("\n\n");
############

num_granges <- length(gr);
grange_names <- names(gr);
dmp_pat_genes_list <- list();

############
# Debugging.
cat("\nnum_granges: ", num_granges, "\n");
cat("\nignore_strand: ", ignore_strand, "\n");
cat("\nopt$output_type: ", opt$output_type, "\n");
cat("\nopt$overlap_type: ", opt$overlap_type, "\n");
cat("\ngrange_names: ", toString(grange_names), "\n");
############

for (i in 1:num_granges) {
    name <- grange_names[[i]];
    grange <- gr[[name]];
    dmp_pat_genes_list[[i]]  <- getDIMPatGenes(GR=grange,
                                               GENES=genes,
                                               type=opt$overlap_type,
                                               ignore.strand=ignore_strand);
#  TODO: the output value is not correct here  output=opt$output_type);
}

names(dmp_pat_genes_list) <- grange_names;

############
# Debugging.
cat("\ndmp_pat_genes_list:\n");
dmp_pat_genes_list
cat("\n\n");
############

# Save the dmp_pat_genes_list.
saveRDS(dmp_pat_genes_list, file=opt$output_rdata, compress=TRUE);

