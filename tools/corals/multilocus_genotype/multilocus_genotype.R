#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("RPostgres"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("vegan"))

option_list <- list(
    make_option(c("--database_connection_string"), action="store", dest="database_connection_string", help="Corals (stag) database connection string"),
    make_option(c("--input_affy_metadata"), action="store", dest="input_affy_metadata", help="Affymetrix 96 well plate input file"),
    make_option(c("--input_pop_info"), action="store", dest="input_pop_info", help="Population information input file"),
    make_option(c("--input_vcf"), action="store", dest="input_vcf", help="VCF input file"),
    make_option(c("--output_stag_db_report"), action="store", dest="output_stag_db_report", help="stag db report output file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_file_path = function(file_name) {
    file_path = paste("output_plots_dir", file_name, sep="/");
    return(file_path);
}

get_database_connection <- function(db_conn_string) {
    # Instantiate database connection.
    # The connection string has this format:
    # postgresql://user:password@host/dbname
    conn_items <- strsplit(db_conn_string, "://")[[1]];
    string_needed <- conn_items[2];
    items_needed <- strsplit(string_needed, "@")[[1]];
    user_pass_string <- items_needed[1];
    host_dbname_string <- items_needed[2];
    user_pass_items <- strsplit(user_pass_string, ":")[[1]];
    host_dbname_items <- strsplit(host_dbname_string, "/")[[1]];
    user <- user_pass_items[1];
    pass <- user_pass_items[2];
    host <- host_dbname_items[1];
    dbname <- host_dbname_items[2];
    # FIXME: is there a way to not hard-code the port?
    conn <- DBI::dbConnect(RPostgres::Postgres(), host=host, port='5432', dbname=dbname, user=user, password=pass);
    return (conn);
}

# Read in VCF input file.
vcf <- read.vcfR(opt$input_vcf);

#Missing GT in samples submitted
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
myMiss <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
myMiss <- (myMiss / nrow(vcf)) * 100;
miss <- data.frame(myMiss);

hets <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))} );
hets <- (hets / nrow(vcf)) * 100;
ht <- data.frame(hets);

refA <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))} );
refA <- (refA / nrow(vcf)) * 100;
rA <- data.frame(refA);

altB <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))} );
altB <- (altB / nrow(vcf)) * 100;
aB <- data.frame(altB);

# Convert VCF file into a genind for the Poppr package.
# TODO: probably should not hard-code 2 cores.
gl <- vcfR2genlight(vcf, n.cores=2);
gind <- new("genind", (as.matrix(gl)));

# Add population information to the genind object.
poptab <- read.table(opt$input_pop_info, check.names=FALSE, header=T, na.strings=c("", "NA"));
colnames(poptab) <- c("row_id", "affy_id", "affy_id2", "user_specimen_id", "user_specimen_id2", "region");
gind@pop <- as.factor(poptab$region);

# Convert genind object to a genclone object.
obj2 <- as.genclone(gind);

# Calculate the bitwise distance between individuals.
xdis <- bitwise.dist(obj2);

# Multilocus genotypes (threshold of 16%).
mlg.filter(obj2, distance=xdis) <- 0.016;
m <- mlg.table(obj2, background=TRUE, color=TRUE);

# Create table of MLGs.
id <- mlg.id(obj2);
dt <- data.table(id, keep.rownames=TRUE);
setnames(dt, c("id"), c("affy_id"));

# Read user's Affymetrix 96 well plate csv file.
pinfo <- read.csv(opt$input_affy_metadata, stringsAsFactors=FALSE);
pinfo <- pinfo$user_specimen_id;
pi <- data.table(pinfo);
setnames(pi, c("pinfo"), c("user_specimen_id"));

# Connect to database.
conn <- get_database_connection(opt$database_connection_string);

# Import the sample table.
mD <- tbl(conn, "sample");

# Select user_specimen_id and mlg columns.
smlg <- mD %>% select(user_specimen_id, coral_mlg_clonal_id, symbio_mlg_clonal_id);

# Convert to dataframe.
sm <- data.frame(smlg);
sm[sm==""] <- NA;

# Convert missing data into data table.
mi <-setDT(miss, keep.rownames=TRUE)[];
setnames(mi, c("rn"), c("affy_id"));
setnames(mi, c("myMiss"), c("percent_missing_data_coral"));
# Round missing data to two digits.
mi$percent_missing_data_coral <- round(mi$percent_missing_data_coral, digits=2);

# Convert heterozygosity data into data table.
ht <-setDT(ht, keep.rownames=TRUE)[];
setnames(ht, c("rn"), c("affy_id"));
setnames(ht, c("hets"), c("percent_mixed_coral"));
# Round missing data to two digits.
ht$percent_mixed<-round(ht$percent_mixed, digits=2);

# Convert refA data into data.table.
rA <-setDT(rA, keep.rownames=TRUE)[];
setnames(rA, c("rn"), c("affy_id"));
setnames(rA, c("refA"), c("percent_reference_coral"));
# round missing data to two digits.
rA$percent_reference<-round(rA$percent_reference, digits=2);

# Convert altB data into data table.
aB <-setDT(aB, keep.rownames=TRUE)[];
setnames(aB, c("rn"), c("affy_id"));
setnames(aB, c("altB"), c("percent_alternative_coral"));
# Round missing data to two digits.
aB$percent_alternative<-round(aB$percent_alternative, digits=2);

#convert mlg id to data.table format
dt <- data.table(id, keep.rownames=TRUE);
setnames(dt, c("id"), c("affy_id"));

# Transform.
df3 <- dt %>%
    group_by(row_number()) %>%
    dplyr::rename(group='row_number()') %>%
    unnest (user_specimen_id) %>%
    # Join with mlg table.
    left_join(sm %>%
              select("affy_id","coral_mlg_clonal_id"),
              by='affy_id');

# If found in database, group members on previous mlg id.
uniques <- unique(df3[c("group", "coral_mlg_clonal_id")]);
uniques <- uniques[!is.na(uniques$coral_mlg_clonal_id),];
na.mlg <- which(is.na(df3$coral_mlg_clonal_id));
na.group <- df3$group[na.mlg];
df3$coral_mlg_clonal_id[na.mlg] <- uniques$coral_mlg_clonal_id[match(na.group, uniques$group)];

# Determine if the sample mlg matched previous genotyped sample.
df4<- df3 %>%
    group_by(group) %>%
    mutate(DB_match = ifelse(is.na(coral_mlg_clonal_id),"no_match","match"));

# Create new mlg id for samples that did not match those in the database.
none <- unique(df4[c("group", "coral_mlg_clonal_id")]);
none <- none[is.na(none$coral_mlg_clonal_id),];
na.mlg2 <- which(is.na(df4$coral_mlg_clonal_id));
n.g <- df4$group[na.mlg2];
ct <- length(unique(n.g));

# List of new group ids, the sequence starts at the number of
# ids present in df4$coral_mlg_clonal_ids plus 1.  Not sure if
# the df4 file contains all ids.  If it doesn't then look below
# to change the seq() function. 
n.g_ids <- sprintf("HG%04d", seq((sum(!is.na(unique(df4["coral_mlg_clonal_id"]))) + 1), by=1, length=ct));
# This is a key for pairing group with new ids.
rat <- cbind(unique(n.g), n.g_ids);
# this for loop assigns the new id iteratively for all that have NA.
for (i in 1:length(na.mlg2)) {
    df4$coral_mlg_clonal_id[na.mlg2[i]] <- n.g_ids[match(df4$group[na.mlg2[i]], unique(n.g))];
}

# subset poptab for all samples.
subpop <- poptab[c(2, 4)];

# Merge data frames for final table.
report_user <- pi %>%
    left_join(subpop %>%
        select("affy_id", "user_specimen_id"),
        by='user_specimen_id') %>%
    left_join(df4 %>%
        select("affy_id", "coral_mlg_clonal_id", "DB_match"),
        by='affy_id') %>%
    left_join(mi %>%
        select("affy_id", "percent_missing_data_coral"),
        by='affy_id') %>%
    left_join(ht %>%
        select("affy_id", "percent_mixed_coral"),
        by='affy_id') %>%
    left_join(rA %>%
        select("affy_id", "percent_reference_coral"),
        by='affy_id') %>%
    left_join(aB %>%
        select("affy_id", "percent_alternative_coral"),
        by='affy_id') %>%
    mutate(DB_match = ifelse(is.na(DB_match), "failed", DB_match))%>%
    mutate(coral_mlg_clonal_id=ifelse(is.na(coral_mlg_clonal_id), "failed", coral_mlg_clonal_id))%>%
    ungroup() %>%
    select(-group);

write.csv(report_user, file=paste(opt$output_stag_db_report), quote=FALSE);

# Combine sample information for database.
report_db <- pinfo %>%
    left_join(report_user %>%
        select("user_specimen_id", "affy_id", "coral_mlg_clonal_id", "DB_match",
               "percent_missing_data_coral", "percent_mixed_coral", "percent_reference_coral",
               "percent_alternative_coral"),
        by='user_specimen_id')

# Rarifaction curve.
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("geno_rarifaction_curve.pdf");
pdf(file=file_path, width=10, height=7);
rarecurve(m, ylab="Number of expected MLGs", sample=min(rowSums(m)), border=NA, fill=NA, font=2, cex=1, col="blue");
dev.off();

# Genotype accumulation curve, sample is number of
# loci randomly selected for to make the curve.
dev.new(width=10, height=7);
file_path = get_file_path("geno_accumulation_curve.pdf");
pdf(file=file_path, width=10, height=7);
genotype_curve(gind, sample=5, quiet=TRUE);
dev.off();

# Create a phylogeny of samples based on distance matrices.
palette(brewer.pal(n=12, name='Set3'));
set.seed(999);
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("nj_phylogeny.pdf");
pdf(file=file_path, width=10, height=7);
# Organize branches by clade.
theTree <- obj2 %>%
    aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE) %>%
    ladderize();
plot.phylo(theTree, tip.color=cols[obj2$pop], label.offset=0.0125, cex=0.7, font=2, lwd=4);
# Add a scale bar showing 5% difference..
add.scale.bar(length=0.05, cex=0.65);
nodelabels(theTree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
dev.off();

