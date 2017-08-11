#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-i", "--input"), action="store", dest="input", help="Input .bed.gz file produced by prepMat"),
    make_option(c("-o", "--output"), action="store", dest="output", help="Output file"),
    make_option(c("-w", "--work_dir"), action="store", dest="work_dir", help="Working directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

data_table <- read.table(opt$input)
as.matrix(data_table)
status <- match(data_table[,3], missing)
data_table[,3] <- paste(opt$work_dir, data_table[,1], ".", data_table[,2], ".bed.gz", sep="")
data_table <- data_table[is.na(status)==TRUE,]
write.table(array(data_table, dim=c(length(data_table)/3, 3)), file=opt$output, quote=FALSE, row.names=FALSE, col.names=FALSE)
