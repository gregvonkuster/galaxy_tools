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
    make_option(c("--output_grangelist"), action="store", dest="output_grangelist", help="Output grangelist file"),
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

## Tom commented out so that TAIR10/arabidopsis can work (TAIR10.gff3 has no "biotype" nor chr13:
#genes <- gene_annot[ gene_annot$type == "gene", c( "gene_id", "biotype", "Name" ) ];
#genes <- genes[ genes$biotype == "protein_coding", "gene_id" ];
#seqlevels(genes) <- "chr13";
## Tom added in so that object "genes" exists:/citation>
genes <- gene_annot;
gr <- readRDS(opt$input);

############
# Debugging.
cat("\ngenes:\n");
show(genes);
cat("\n\n");
cat("\nclass(genes):\n");
class(genes);
cat("\n\n");
cat("\n\ngr:\n");
gr
cat("\n\n");
cat("\nclass(gr):\n");
class(gr);
cat("\n\n");
cat("\nignore_strand: ", ignore_strand, "\n");
cat("\nopt$output_type: ", opt$output_type, "\n");
cat("\nopt$overlap_type: ", opt$overlap_type, "\n");
############

if (is(gr, "GRanges")) {
    # The getDIMPatGenes.GRanges function signature
    # does not include the output parameter.
    # The ret_val object will be a GRanges or a list
    # depending on the value of the opt$output_type
    # parameter.
    ret_val <- getDIMPatGenes(GR=gr, GENES=genes, type=opt$overlap_type, ignore.strand=ignore_strand);
} else {
    ret_val <- getDIMPatGenes(GR=gr, GENES=genes, type=opt$overlap_type, ignore.strand=ignore_strand, output=opt$output_type);
}

############
# Debugging.
cat("\nret_val:\n");
ret_val
cat("\nclass(ret_val):\n");
class(ret_val)
cat("\n\n");
############

# Save the GRanges.
saveRDS(ret_val, file=opt$output_grangelist, compress=TRUE);

