#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dbplyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("rworldmap"))
suppressPackageStartupMessages(library("RPostgres"))
suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("yarrr"))
theme_set(theme_bw())

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
    conn <- DBI::dbConnect(RPostgres::Postgres(), host=host, port="5432", dbname=dbname, user=user, password=pass);
    return (conn);
}

time_elapsed <- function(start_time) {
    cat("Elapsed time: ", proc.time() - start_time, "\n\n");
}

time_start <- function(msg) {
    start_time <- proc.time();
    cat(msg, "...\n");
    return(start_time);
}

# Read in VCF input file.
start_time <- time_start("Reading VCF input");
vcf <- read.vcfR(opt$input_vcf);
time_elapsed(start_time);

# Convert VCF file into a genind for the Poppr package.
start_time <- time_start("Converting VCF to a genind object");
genind_obj <- vcfR2genind(vcf);
time_elapsed(start_time);

# Add population information to the genind object.
population_info_data_table <- read.table(opt$input_pop_info, check.names=FALSE, header=F, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t");
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region");
genind_obj@pop <- as.factor(population_info_data_table$region);
strata(genind_obj) <- data.frame(pop(genind_obj));

# Convert genind object to a genclone object.
start_time <- time_start("Converting genind to genclone object");
genind_clone <- as.genclone(genind_obj);
time_elapsed(start_time);

# Calculate the bitwise distance between individuals.
start_time <- time_start("Calculating the bitwise distance between individuals");
bitwise_distance <- bitwise.dist(genind_clone);
time_elapsed(start_time);

# Multilocus genotypes (threshold of 3.2%).
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032;
m <- mlg.table(genind_clone, background=TRUE, color=TRUE);

# Create list of MLGs.
mlg_ids <- mlg.id(genind_clone);

# Read user's Affymetrix 96 well plate tabular file.
affy_metadata_data_frame <- read.table(opt$input_affy_metadata, header=FALSE, stringsAsFactors=FALSE, sep="\t", na.strings = c("", "NA"));
colnames(affy_metadata_data_frame) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                                        "region", "latitude", "longitude", "geographic_origin", "sample_location",
                                        "latitude_outplant", "longitude_outplant", "depth", "disease_resist",
                                        "bleach_resist", "mortality","tle", "spawning", "collector_last_name",
                                        "collector_first_name", "organization", "collection_date", "email", "seq_facility",
                                        "array_version", "public", "public_after_date", "sperm_motility", "healing_time",
                                        "dna_extraction_method", "dna_concentration", "registry_id");
affy_metadata_data_frame$user_specimen_id <- as.character(affy_metadata_data_frame$user_specimen_id);
user_specimen_ids <- as.character(affy_metadata_data_frame$user_specimen_id);
# The specimen_id_field_call_data_table looks like this:
# user_specimen_ids V2
# test_002          prolifera
# test_003          prolifera
specimen_id_field_call_data_table <- data.table(user_specimen_ids, affy_metadata_data_frame$field_call);
# Rename the user_specimen_ids column.
setnames(specimen_id_field_call_data_table, c("user_specimen_ids"), c("user_specimen_id"));
# Rename the V2 column.
setnames(specimen_id_field_call_data_table, c("V2"), c("field_call"));

# Connect to database.
conn <- get_database_connection(opt$database_connection_string);
# Import the sample table.
sample_table <- tbl(conn, "sample");
# Import the genotype table.
genotype_table <- tbl(conn, "genotype");
# Select columns from the sample table and the
# genotype table joined by genotype_id.
sample_table_columns <- sample_table %>% select(user_specimen_id, affy_id, genotype_id);
smlg <- sample_table_columns %>%
    left_join(genotype_table %>%
              select("id", "coral_mlg_clonal_id", "symbio_mlg_clonal_id"),
              by=c("genotype_id"="id"));
# Convert to dataframe.
smlg_data_frame <- data.frame(smlg);
# Name the columns.
colnames(smlg_data_frame) <- c("user_specimen_id", "affy_id", "genotype_id", "coral_mlg_clonal_id", "symbio_mlg_clonal_id");

# Missing GT in samples submitted.
start_time <- time_start("Discovering missing GT in samples");
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
missing_gt <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
missing_gt <- (missing_gt / nrow(vcf)) * 100;
missing_gt_data_frame <- data.frame(missing_gt);
# The specimen_id_field_call_data_table looks like this:
# rn                                 missing_gt
# a100000-4368120-060520-256_I07.CEL 0.06092608
# a100000-4368120-060520-256_K07.CEL 0.05077173
missing_gt_data_table <-setDT(missing_gt_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(missing_gt_data_table, c("rn"), c("affy_id"));
# Rename the missing_gt column.
setnames(missing_gt_data_table, c("missing_gt"), c("percent_missing_data_coral"));
# Round data to two digits.
missing_gt_data_table$percent_missing_data_coral <- round(missing_gt_data_table$percent_missing_data_coral, digits=2);
time_elapsed(start_time);

# Heterozygous alleles.
start_time <- time_start("Discovering heterozygous alleles");
heterozygous_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))});
heterozygous_alleles <- (heterozygous_alleles / nrow(vcf)) * 100;
heterozygous_alleles_data_frame <- data.frame(heterozygous_alleles);
# The heterozygous_alleles_data_table looks like this:
# rn                                 heterozygous_alleles
# a100000-4368120-060520-256_I07.CEL 73.94903
# a100000-4368120-060520-256_K07.CEL 74.40089
heterozygous_alleles_data_table <- setDT(heterozygous_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(heterozygous_alleles_data_table, c("rn"), c("affy_id"));
# Rename the heterozygous_alleles column.
setnames(heterozygous_alleles_data_table, c("heterozygous_alleles"), c("percent_heterozygous_coral"));
# Round data to two digits.
heterozygous_alleles_data_table$percent_heterozygous_coral <- round(heterozygous_alleles_data_table$percent_heterozygous_coral, digits=2);
time_elapsed(start_time);

# Reference alleles.
start_time <- time_start("Discovering reference alleles");
reference_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))});
reference_alleles <- (reference_alleles / nrow(vcf)) * 100;
reference_alleles_data_frame <- data.frame(reference_alleles);
# The reference_alleles_data_table looks like this:
# rn                                 reference_alleles
# a100000-4368120-060520-256_I07.CEL 11.60642
# a100000-4368120-060520-256_K07.CEL 11.45918
reference_alleles_data_table  <- setDT(reference_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(reference_alleles_data_table, c("rn"), c("affy_id"));
# Rename the reference_alleles column.
setnames(reference_alleles_data_table, c("reference_alleles"), c("percent_reference_coral"));
# Round data to two digits.
reference_alleles_data_table$percent_reference <- round(reference_alleles_data_table$percent_reference, digits=2);
time_elapsed(start_time);

# Alternative alleles
start_time <- time_start("Discovering alternative alleles");
alternative_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))});
alternative_alleles <- (alternative_alleles / nrow(vcf)) * 100;
alternative_alleles_data_frame <- data.frame(alternative_alleles);
# The alternative_alleles_data_table looks like this:
# rn                                 alternative_alleles
# a100000-4368120-060520-256_I07.CEL 14.38363
# a100000-4368120-060520-256_K07.CEL 14.08916
alternative_alleles_data_table <- setDT(alternative_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(alternative_alleles_data_table, c("rn"), c("affy_id"));
# Rename the alternative_alleles column.
setnames(alternative_alleles_data_table, c("alternative_alleles"), c("percent_alternative_coral"));
# Round data to two digits.
alternative_alleles_data_table$percent_alternative <- round(alternative_alleles_data_table$percent_alternative, digits=2);
time_elapsed(start_time);

# The alternative_alleles_data_table looks like this:
# mlg_ids
# a550962-4368120-060520-500_M23.CEL
# a550962-4368120-060520-256_A19.CEL
mlg_ids_data_table <- data.table(mlg_ids, keep.rownames=TRUE);
# Rename the mlg_ids column.
setnames(mlg_ids_data_table, c("mlg_ids"), c("affy_id"));

# sample_mlg_tibble looks like this:
# A tibble: 262 x 3
# Groups:   group [?]
# group affy_id                            coral_mlg_clonal_id
# <int> <chr>                              <chr>              
# 1     a550962-4368120-060520-500_M23.CEL NA                 
# 2     a550962-4368120-060520-256_A19.CEL HG0006
sample_mlg_tibble <- mlg_ids_data_table %>%
    group_by(row_number()) %>%
    dplyr::rename(group="row_number()") %>%
    unnest (affy_id) %>%
    # Join with mlg table.
    left_join(smlg_data_frame %>%
              select("affy_id","coral_mlg_clonal_id"),
              by="affy_id");

# If found in database, group members on previous mlg id.
uniques <- unique(sample_mlg_tibble[c("group", "coral_mlg_clonal_id")]);
uniques <- uniques[!is.na(uniques$coral_mlg_clonal_id),];
na.mlg <- which(is.na(sample_mlg_tibble$coral_mlg_clonal_id));
na.group <- sample_mlg_tibble$group[na.mlg];
sample_mlg_tibble$coral_mlg_clonal_id[na.mlg] <- uniques$coral_mlg_clonal_id[match(na.group, uniques$group)];

# Determine if the sample mlg matched previous genotyped sample.
# sample_mlg_match_tibble looks like this:
# A tibble: 262 x 4
# Groups:   group [230]
# group affy_id                            coral_mlg_clonal_id DB_match
# <int> <chr>                              <chr>               <chr>   
# 1     a550962-4368120-060520-500_M23.CEL NA                  no_match
# 2     a550962-4368120-060520-256_A19.CEL HG0006              match 
sample_mlg_match_tibble <- sample_mlg_tibble %>%
    group_by(group) %>%
    mutate(DB_match = ifelse(is.na(coral_mlg_clonal_id),"no_match", "match"));

# Create new mlg id for samples with no matches in the database.
none <- unique(sample_mlg_match_tibble[c("group", "coral_mlg_clonal_id")]);
none <- none[is.na(none$coral_mlg_clonal_id),];
na.mlg2 <- which(is.na(sample_mlg_match_tibble$coral_mlg_clonal_id));
n.g <- sample_mlg_match_tibble$group[na.mlg2];
ct <- length(unique(n.g));

# List of new group ids, the sequence starts at the number of
# ids present in sample_mlg_match_tibble$coral_mlg_clonal_ids
# plus 1.
# FIXME: Not sure if # the sample_mlg_match_tibble file
# contains all ids.  If it doesn't then look below to change
# the seq() function.
n.g_ids <- sprintf("HG%04d", seq((sum(!is.na(unique(sample_mlg_match_tibble["coral_mlg_clonal_id"]))) + 1), by=1, length=ct));

#####################################
# FIXME: the following code is commented because it is not used.
# Pair group with new ids.
# rat looks like this:
#             n.g_ids 
#  [1,] "1"   "HG0135"
#  [2,] "10"  "HG0136"
#rat <- cbind(unique(n.g), n.g_ids);
#####################################

# Assign the new id iteratively for all that have NA.
for (i in 1:length(na.mlg2)) {
    sample_mlg_match_tibble$coral_mlg_clonal_id[na.mlg2[i]] <- n.g_ids[match(sample_mlg_match_tibble$group[na.mlg2[i]], unique(n.g))];
}

# Subset population_info_data_table for all samples.
# affy_id_user_specimen_id_vector looks like this:
# affy_id                            user_specimen_id
# a100000-4368120-060520-256_I07.CEL 13704
# a100000-4368120-060520-256_K07.CEL 13706
affy_id_user_specimen_id_vector <- population_info_data_table[c(2, 3)];

# Merge data frames for final table.
start_time <- time_start("Merging data frames");
stag_db_report <- specimen_id_field_call_data_table %>%
    left_join(affy_id_user_specimen_id_vector %>%
        select("affy_id", "user_specimen_id"),
        by="user_specimen_id") %>%
    left_join(sample_mlg_match_tibble %>%
        select("affy_id", "coral_mlg_clonal_id", "DB_match"),
        by="affy_id") %>%
    left_join(missing_gt_data_table %>%
        select("affy_id", "percent_missing_data_coral"),
        by="affy_id") %>%
    left_join(heterozygous_alleles_data_table %>%
        select("affy_id", "percent_heterozygous_coral"),
        by="affy_id") %>%
    left_join(reference_alleles_data_table %>%
        select("affy_id", "percent_reference_coral"),
        by="affy_id") %>%
    left_join(alternative_alleles_data_table %>%
        select("affy_id", "percent_alternative_coral"),
        by="affy_id") %>%
    mutate(DB_match = ifelse(is.na(DB_match), "failed", DB_match))%>%
    mutate(coral_mlg_clonal_id = ifelse(is.na(coral_mlg_clonal_id), "failed", coral_mlg_clonal_id)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 40 & percent_alternative_coral <= 44.5, "A.palmata","other")) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 45.5 & percent_alternative_coral <= 50, "A.cervicornis", genetic_coral_species_call)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_heterozygous_coral > 40, "A.prolifera", genetic_coral_species_call)) %>%
    ungroup() %>%
    select(-group);
time_elapsed(start_time);

start_time <- time_start("Writing csv output");
write.csv(stag_db_report, file=opt$output_stag_db_report, quote=FALSE);
time_elapsed(start_time);

# Database sample table.
sample_db <- affy_metadata_data_frame %>%
  left_join(stag_db_report %>%
      select("user_specimen_id","affy_id", "percent_missing_data_coral", "percent_heterozygous_coral", "percent_reference_coral", "percent_alternative_coral"),
      by='user_specimen_id');

# Representative clone for genotype table.
start_time <- time_start("Creating representative clone for genotype table");
no_dup_genotypes_genind <- clonecorrect(genind_clone, strata = ~pop.genind_obj.);
id_rep <- mlg.id(no_dup_genotypes_genind);
id_data_table <- data.table(id_rep, keep.rownames=TRUE);
# Rename the id_rep column.
setnames(id_data_table, c("id_rep"), c("affy_id"));
time_elapsed(start_time);

# # Combine with previously genotyped samples in the database.
start_time <- time_start("Selecting from various database tables");
representative_mlg_tibble <- id_data_table %>%
    group_by(row_number()) %>%
    rename(group='row_number()') %>%
    unnest(affy_id) %>%
    left_join(stag_db_report %>%
              select("coral_mlg_clonal_id", "user_specimen_id", "affy_id"),
              by='affy_id') %>%
    mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_clonal_id), "", affy_id)) %>%
    ungroup() %>%
    select(-group);

# Database genotype table.
genotype_table_join <- sample_mlg_match_tibble %>%
    left_join(representative_mlg_tibble %>%
    select("affy_id", "coral_mlg_rep_sample_id", "user_specimen_id"),
    by='affy_id') %>%
    ungroup() %>%
    select(-group);

# Database taxonomy table.
taxonomy_table_join <- stag_db_report %>%
    select(genetic_coral_species_call, affy_id) %>%
    mutate(genus_name = ifelse(genetic_coral_species_call == genetic_coral_species_call[grep("^A.*", genetic_coral_species_call)], "Acropora", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.palmata", "palmata", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.cervicornis", "cervicornis", species_name)) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.prolifera", "prolifera", species_name));
time_elapsed(start_time);
# Table of alleles for the new samples subset to new plate data.
# Create vector indicating number of individuals desired from
# affy_id column of stag_db_report data table.
i <- ifelse(is.na(stag_db_report[1]), "", stag_db_report[[1]]);
i <- i[!apply(i== "", 1, all),];
sample_alleles_vector <- genind_clone[i, mlg.reset=FALSE, drop=FALSE];

#####################################
# FIXME: the following code is commented because it is not used.
# Convert to data frame.
#at_96 <- genind2df(sample_alleles_vector, sep="");
#at_96 <- at_96 %>%
#select(-pop);
# Allele string for Allele.table in database.
#uat_96 <- unite(at_96, alleles, 1:19696, sep=" ", remove=TRUE);
#uat_96 <- setDT(uat_96, keep.rownames=TRUE)[];
#setnames(uat_96, c("rn"), c("user_specimen_id"));
#####################################

# Create a phylogeny tree of samples based on distance matrices.
# cols looks like this:
#       blue1         red       green        pink      orange       blue2 
# "#0C5BB0FF" "#EE0011FF" "#15983DFF" "#EC579AFF" "#FA6B09FF" "#149BEDFF" 
#      green2      yellow   turquoise        poop 
# "#A1C720FF" "#FEC10BFF" "#16A08CFF" "#9A703EFF"
cols <- piratepal("basel");
set.seed(999);

# Start PDF device driver.
start_time <- time_start("Creating nj_phylogeny_tree.pdf");
dev.new(width=10, height=7);
file_path = get_file_path("nj_phylogeny_tree.pdf");
pdf(file=file_path, width=10, height=7);
# Organize branches by clade.
nj_phylogeny_tree <- sample_alleles_vector %>%
    aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE) %>%
    ladderize();
nj_phylogeny_tree$tip.label <- stag_db_report$user_specimen_id[match(nj_phylogeny_tree$tip.label, stag_db_report$affy_id)];
plot.phylo(nj_phylogeny_tree, tip.color=cols[sample_alleles_vector$pop], label.offset=0.0125, cex=0.3, font=2, lwd=4, align.tip.label=F, no.margin=T);
# Add a scale bar showing 5% difference.
add.scale.bar(0, 0.95, length=0.05, cex=0.65, lwd=3);
nodelabels(nj_phylogeny_tree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
legend("topright", legend=c(levels(sample_alleles_vector$pop)), text.col=cols, xpd=T, cex=0.8);
dev.off()
time_elapsed(start_time);

# Subset VCF to the user samples.
start_time <- time_start("Subsetting vcf to the user samples");
l <- length(i);
n <- ncol(vcf@gt);
s <- n - l;
svcf <- vcf[, s:n];
write.vcf(svcf, "subset.vcf.gz");
vcf.fn <- "subset.vcf.gz";
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only");
genofile <- snpgdsOpen(filename="test3.gds", readonly=FALSE);
gds_array <- read.gdsn(index.gdsn(genofile, "sample.id"));
# gds_array looks like this:
# [1] "a550962-4368120-060520-500_A03.CEL" "a550962-4368120-060520-500_A05.CEL"
# [3] "a550962-4368120-060520-500_A09.CEL" "a550962-4368120-060520-500_A11.CEL"
gds_data_frame <- data.frame(gds_array);
# gds_data_frame looks like this:
# gds_array
# a550962-4368120-060520-500_A03.CEL
# a550962-4368120-060520-500_A05.CEL
gds_data_table <- setDT(gds_data_frame, keep.rownames=FALSE)[];
# Rename the gds_array column.
setnames(gds_data_table, c("gds_array"), c("affy_id"));
# affy_id_region_list looks like this:
# affy_id                            region
# a100000-4368120-060520-256_I07.CEL USVI
# a100000-4368120-060520-256_K07.CEL USVI
affy_id_region_list <- population_info_data_table[c(2,4)];
gds_data_table_join <- gds_data_table %>%
    left_join(affy_id_region_list %>%
        select("affy_id", "region"),
        by='affy_id')%>%
        drop_na();
samp.annot <- data.frame(pop.group=c(gds_data_table_join$region));
add.gdsn(genofile, "sample.annot", samp.annot);
# population_code looks like this:
# [1] 18.361733   18.361733   18.361733   18.361733   18.361733   18.361733  
# [7] 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009
population_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"));
pop.group <- as.factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")));
# pop.group looks like this:
# [1] 18.361733   18.361733   18.361733   18.361733   18.361733   18.361733  
# [7] 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009
time_elapsed(start_time);

# Distance matrix calculation.
start_time <- time_start("Calculating distance matrix");
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE);
time_elapsed(start_time);

# Cluster analysis on the genome-wide IBS pairwise distance matrix.
start_time <- time_start("Clustering the genome-wide IBS pairwise distance matrix");
set.seed(100);
par(cex=0.6, cex.lab=1, cex.axis=1.5,cex.main=2);
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, autosome.only=FALSE));
time_elapsed(start_time);

# Default clustering.
start_time <- time_start("Creating ibs_default.pdf");
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("ibs_default.pdf");
pdf(file=file_path, width=10, height=7);
rv <- snpgdsCutTree(ibs.hc, col.list=cols, pch.list=15);
snpgdsDrawTree(rv, main="Color by Cluster", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(rv$samp.group), xpd=T, col=cols[1:nlevels(rv$samp.group)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# Color cluster by region.
start_time <- time_start("Creating ibs_region.pdf");
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("ibs_region.pdf");
pdf(file=file_path, width=10, height=7);
race <- as.factor(population_code);
rv2 <- snpgdsCutTree(ibs.hc, samp.group=race,col.list=cols, pch.list=15);
snpgdsDrawTree(rv2, main="Color by Region", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(race), xpd=T, col=cols[1:nlevels(race)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# close GDS file.
snpgdsClose(genofile);

# Sample MLG on a map.
start_time <- time_start("Creating mlg_map.pdf");
world_map_spatial_data_frame <- getMap(resolution="low");
world_map_data_frame <- as.data.frame(world_map_spatial_data_frame);
affy_metadata_data_frame$mlg <- stag_db_report$coral_mlg_clonal_id;
n <- nrow(affy_metadata_data_frame);
mxlat <- max(affy_metadata_data_frame$latitude, na.rm=TRUE);
mnlat <- min(affy_metadata_data_frame$latitude, na.rm=TRUE);
mxlong <- max(affy_metadata_data_frame$longitude, na.rm=TRUE);
mnlong <- min(affy_metadata_data_frame$longitude, na.rm=TRUE);
# FIXME: figure out a way to replace the sf calls below.
#world_map_prep <- ggplot(data=world_map_data_frame) + geom_sf() + coord_sf(xlim=c(mnlong-3, mxlong+3), ylim=c(mnlat-3, mxlat+3), expand=FALSE);
world_map_prep <- ggplot(data=world_map_data_frame, ylim=c(mnlat-3, mxlat+3), expand=FALSE);
colourCount = length(unique(affy_metadata_data_frame$mlg));
getPalette = colorRampPalette(piratepal("basel"));
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("mlg_map.pdf");
pdf(file=file_path, width=10, height=7);
world_map_plot <- world_map_prep +
      geom_point(data=affy_metadata_data_frame, aes(x=longitude, y=latitude, group=mlg, color=mlg), alpha=.7, size=3) +
      scale_color_manual(values=getPalette(colourCount)) +
      theme(legend.position="bottom") +
      guides(color=guide_legend(nrow=8,byrow=F));
dev.off()
time_elapsed(start_time);

# Missing data barplot.
start_time <- time_start("Creating missing_data.pdf");
population_info_data_table$miss <- stag_db_report$percent_missing_data_coral[match(missing_gt_data_frame$affy_id, stag_db_report$affy_id)];
test2 <- which(!is.na(population_info_data_table$miss));
miss96 <- population_info_data_table$miss[test2];
name96 <- population_info_data_table$user_specimen_id[test2];
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("missing_data.pdf");
pdf(file=file_path, width=10, height=7);
par(mar = c(8, 4, 4, 2));
x <- barplot(miss96, las=2, col=cols, ylim=c(0, 3), cex.axis=0.8, space=0.8, ylab="Missingness (%)", xaxt="n");
text(cex=0.6, x=x-0.25, y=-.05, name96, xpd=TRUE, srt=60, adj=1);
dev.off()
time_elapsed(start_time);

# Generate a pie chart for each sample with a genotype.
# Store the numerical and user_specimen_id values from
# stag_db_report for the charts (user_specimen_id names
# will be used to label each chart).
start_time <- time_start("Creating percent_breakdown.pdf");
dt1 <- data.table(stag_db_report);
dt1 <- stag_db_report[c(-2, -3, -4)];
dt1 <- na.omit(dt1);
# Translate to N (i.e., number of samples with a
# genotype) columns and 5 rows.
tdt1 <- t(dt1);
# Make another data table and transpose it the same as dt1 to
# get numerics. These will feed into the creation of N vectors.
dt2 <- data.table(stag_db_report);
dt2 <- stag_db_report[c(-1, -2, -3, -4)];
# Translate to N columns and 5 rows.
tdt2 <- t(dt2);
tdt1_matrix <- as.matrix(tdt1[-1,]);
# The number of columns is the number of samples with genotypes.
nc <- ncol(tdt1_matrix);
mode(tdt1_matrix) <- "numeric";
spy <- rowMeans(tdt1_matrix);
dev.new(width=10, height=7);
file_path = get_file_path("percent_breakdown.pdf");
pdf(file=file_path, width=10, height=7);
# Average pie of all samples.
labels <- paste(c("missing data", "mixed", "reference", "alternative"), " (", round(spy, 1), "%)", sep="");
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
main <- "Average breakdown of SNP assignments across all samples";
pie(spy, labels=labels, radius=0.60, col=col, main=main, cex.main=.75);
par(mfrow=c(3, 2));
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
for (i in 1:nc) {
    tmp_labels <- paste(labels, " (", round(tdt1_matrix[,i], 1), "%)", sep="");
    main <- paste("Breakdown of SNP assignments for", tdt1[1, i]);
    pie(tdt1_matrix[,i], labels=tmp_labels, radius=0.90, col=col, main=main, cex.main=.85, cex=0.75);
}
dev.off()
time_elapsed(start_time);
