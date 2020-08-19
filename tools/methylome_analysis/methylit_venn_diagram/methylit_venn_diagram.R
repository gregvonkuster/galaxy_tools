#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VennDiagram"))

option_list <- list(
    make_option(c("--input"), action="store", dest="input", help="pDMP file"),
    make_option(c("--output"), action="store", dest="output", help="Output file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Here df should be a data frame with 2 columns,
# Category and HelingerDivergence
pdmp <- readRDS(opt$input);

############
# Debugging.
cat("\nInput pDMP:\n\n");
pdmp
cat("\n\n");
############

num_granges <- length(pdmp);
grange_names <- names(pdmp);

############
# Debugging.
cat("\nnum_granges: ", num_granges, "\n");
cat("\ngrange_names: ", grange_names, "\n");
cat("\n\n");
############

# Start PDF device driver.
dev.new(width=20, height=20);
pdf(file=opt$output, width=20, height=20);

if (num_granges == 3) {
    # The VennDiagram package has routines for drawing
    # 1 to 5 sets (https://rdrr.io/cran/VennDiagram/man/).
    # We'll currently support drawing Venn Diagrams
    # with 3 sets.
    # TODO: Add support for additional sized sets.
    name1 <- grange_names[[1]];
    name2 <- grange_names[[2]];
    name3 <- grange_names[[3]];
    n12 <- length(GenomicRanges::intersect(pdmp[[name1]], pdmp[[name2]]));
    n23 <- length(GenomicRanges::intersect(pdmp[[name1]], pdmp[[name3]]));
    n13 <- length(GenomicRanges::intersect(pdmp[[name2]], pdmp[[name3]]));
    n123 <- length(Reduce(GenomicRanges::intersect, list(pdmp[[name1]], pdmp[[name2]], pdmp[[name3]])));

    ############
    # Debugging.
    cat("\nname1: ", name1, "\n");
    cat("\nname2: ", name2, "\n");
    cat("\nname3: ", name3, "\n");
    cat("\nn12: ", n12, "\n");
    cat("\nn23: ", n23, "\n");
    cat("\nn13: ", n13, "\n");
    cat("\nn123: ", n123, "\n");
    cat("\n\n");
    ############

    grid.newpage();
    v = draw.triple.venn(area1 = length(pdmp[[name1]]),
                         area2 = length(pdmp[[name2]]),
                         area3 = length(pdmp[[name3]]),
                         n12 = n12,
                         n23 = n23,
                         n13 = n13,
                         n123 = n123,
                         category = grange_names,
                         lty = rep("blank", 3),
                         fill = c("blue", "yellow", "magenta"),
                         alpha = c(0.1, 0.2, 0.3),
                         cat.pos = c(-80, 90, 0),
                         cat.col = c("blue", "darkorange", "red"),
                         cat.dist = c( -0.1, -0.08, -0.26),
                         cex = rep(1.7, 7),
                         cat.cex = c( 1.5, 1.5, 1.5),
                         label.col = c( "blue", "darkorange", "darkorange", "red", "white", "red", "red"),
                         scaled = TRUE);
}
grid.draw(v);
dev.off();

