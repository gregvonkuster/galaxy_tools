#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-b", "--build"), action="store", dest="build", help="Genome build"),
    make_option(c("-c", "--chrom_len_file"), action="store", dest="chrom_len_file", help="Chromosome length file"),
    make_option(c("-d", "--header"),  action="store", dest="header", default=NULL, help="Track header"),
    make_option(c("-e", "--state_name"),  action="store", dest="state_name", help="State name"),
    make_option(c("-f", "--hub_id"),  action="store", dest="hub_id", help="Not sure what this is"),
    make_option(c("-g", "--email"),  action="store", dest="email", help="User email address"),
    make_option(c("-i", "--cell_info"),  action="store", dest="cell_info", default=NULL, help="Not sure what this is"),
    make_option(c("-n", "--hub_name"),  action="store", dest="hub_name", default=NULL, help="Not sure what this is"),
    make_option(c("-p", "--input_dir_para"), action="store", dest="input_dir_para", help="Directory containing .para outputs from IDEAS"),
    make_option(c("-q", "--input_dir_state"), action="store", dest="input_dir_state", help="Directory containing .state outputs from IDEAS"),
    make_option(c("-u", "--output_track_db"),  action="store", dest="output_track_db", help="Output track db file"),
    make_option(c("-w", "--output_trackhub"),  action="store", dest="output_trackhub", help="Output hub file"),
    make_option(c("-x", "--output_trackhub_files_path"),  action="store", dest="output_trackhub_files_path", help="Output hub extra files path")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

create_color_scheme = function(input_dir_para) {
    # Create the color scheme.
    para_files <- list.files(path=input_dir_para, full.names=TRUE);
    mc = NULL;
    x = read.table(para_files[1], comment="!", nrows=1);
    l = as.integer(regexpr("\\*", as.matrix(x)));
    l = min(which(l>0)) - 2;
    x = as.matrix(read.table(para_files[1]));
    if (length(para_files) > 1) {
        for (i in 2:length(para_files)) {
            x = x + as.matrix(read.table(para_files[i]));
        }
    }
    p = x[,1] / sum(x[,1]);
    m = array(as.matrix(x[,1:l+1] / x[,1]), dim=c(dim(x)[1], l));
    state_color <- get_state_color(m, mc);
    return(state_color);
}

get_state_color = function(statemean, markcolor=NULL)
{   
    if (length(markcolor) == 0) {
        markcolor = rep("", dim(statemean)[2]);
        markcolor[order(apply(statemean, 2, sd), decreasing=T)] = hsv((1:dim(statemean)[2]-1)/dim(statemean)[2], 1, 1);
        markcolor = t(col2rgb(markcolor));
    }
    rg = apply(statemean, 1, range);
    mm = NULL;
    for(i in 1:dim(statemean)[1]) {
       mm = rbind(mm, (statemean[i,]-rg[1,i])+1e-10)/(rg[2,i]-rg[1,i]+1e-10));
    }
    mm = mm^6;
    if(dim(mm)[2] > 1) {
        mm = mm / (apply(mm, 1, sum) + 1e-10);
    }
    mycol = mm%*%markcolor;
    s = apply(statemean, 1, max);
    s = (s - min(s)) / (max(s) - min(s) + 1e-10);
    h = t(apply(mycol, 1, function(x) {rgb2hsv(x[1], x[2], x[3])}));
    h[,2] = h[,2] * s;
    h = apply(h, 1, function(x) {hsv(x[1], x[2], x[3])});
    rt = cbind(apply(t(col2rgb(h)), 1, function(x) {paste(x, collapse=",")}), h);
    return(rt);
}

create_track = function(input_dir_state, chrom_len_file, base_track_file_name, state_color, state_name, header) {
    state_files <- list.files(path=input_dir_state, full.names=TRUE);
    genome_size = read.table(chrom_len_file);
    g = NULL;
    for(i in state_files) {
        tg = as.matrix(fread(i));
        t = NULL;
        for(j in 1:dim(genome_size)[1]) {
            t = c(t, which(tg[,2]==as.character(genome_size[j, 1]) & as.numeric(tg[,4]) > as.numeric(genome_size[j, 2])));
        }   
        if (length(t) > 0) {
            tg = tg[-t,];
        }
        t = which(is.na(match(tg[,2], genome_size[,1]))==T);
        if (length(t)>0) {
            tg = tg[-t,];
        }
        g = rbind(g, tg);
    }
    uchr = sort(unique(as.character(g[,2])));
    g1 = NULL;
    for(i in uchr) {
        t = which(g[,2]==i);
        g1 = rbind(g1, g[t[order(as.integer(g[t, 3]))],]);
    }
    g = NULL;
    chr = as.character(g1[,2]);
    posst = as.numeric(g1[,3]);
    posed = as.numeric(g1[,4]);
    state = as.matrix(g1[,5:(dim(g1)[2]-1)]);
    if (length(state_name)==0) {
        state_name=0:max(state);
    }
    L = dim(g1)[1];
    n = dim(state)[2];
    if (length(header) > 0) {
        colnames(g1) = header;
    }
    cells = as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
    g1 = NULL;
    options(scipen=999);
    tt = which(chr[2:L]!=chr[2:L-1]);
    tt = c(tt, which(posst[2:L]!=posed[2:L-1]));
    tt = sort(unique(tt));
    for(i in 1:n) {
        tstate = state[,i];
        t = c(tt, which(tstate[2:L]!=tstate[2:L-1]));
        t = sort(unique(t));
        t0 = c(0, t) + 1;
        t = c(t, L);
        np = cbind(chr[t], posst[t0], posed[t], tstate[t]);
        x = cbind(np[,1:3], state_name[as.integer(np[,4])+1], 1000, ".", np[,2:3], state_color[as.numeric(np[,4])+1]);
        track_file_name_bed <- paste(base_track_file_name, i, "bed", sep=".")
        track_file_name_bigbed <- paste(base_track_file_name, i, "bigbed", sep=".")
        write.table(as.matrix(x), track_file_name_bed, quote=F, row.names=F, col.names=F);
        system(paste("bedToBigBed ", track_file_name_bed, chrom_len_file, " ", track_file_name_bigbed));
        system(paste("rm ", track_file_name_bed));
    }
    return(cells);
}

create_track_db = function(input_dir_state, chrom_len_file, header, tracks_dir, hub_id, hub_name, state_color, state_name, cell_info) {
    base_track_file_name <- paste(tracks_dir, hub_id, sep="");
    cells = create_track(input_dir_state, chrom_len_file, base_track_file_name, state_color, state_name, header);
    if (length(cell_info) == 0) {
        cell_info = cbind(cells, cells, cells, "#000000");
        cell_info = array(cell_info, dim=c(length(cells), 4));
    }
    cell_info = as.matrix(cell_info);
    track_db = NULL;
    for(i in 1:length(cells)) {
        ii = which(cells[i] == cell_info[,1]);
        if (length(ii) == 0) {
            next;
        }
        ii = ii[1];
        track_db = c(track_db, paste("track bigBed", i, sep=""));
        track_db = c(track_db, paste("priority", ii));
        track_db = c(track_db, "type bigBed 9 .");
        track_db = c(track_db, "itemRgb on");
        track_db = c(track_db, "maxItems 100000");
        track_db = c(track_db, paste("bigDataUrl ", targetURL, hubid,".",i,".bb",sep=""));
        track_db = c(track_db, paste("shortLabel", cell_info[ii, 2]));
        track_db = c(track_db, paste("longLabel", paste(hub_name, cell_info[ii, 3])));
        track_db = c(track_db, paste("color", paste(c(col2rgb(cell_info[ii, 4])), collapse=",")));
        track_db = c(track_db, "visibility dense");
        track_db = c(track_db, ""); 
    }
    return(track_db);
}

if (length(opt$hub_name) == 0) {
    hub_name <- opt$hub_id;
} else {
    hub_name <- opt$hub_name;
}

# Create the color scheme.
state_color <- create_color_scheme(opt$input_dir_para);

# Create the tracks.
tracks_dir <- paste(opt$output_trackhub_files_path, "/", "tracks", "/", sep="");
dir.create(tracks_dir, showWarnings=FALSE);
track_db <- create_track_db(opt$input_dir_state, opt$chrom_len_file, opt$header, tracks_dir, opt$hub_id, hub_name, state_color, opt$state_name, opt$cell_info);

# Create the primary HTML dataset.


# Create the trackDb.txt output.
track_db_file_path <- paste(tracks_dir, "/", paste("trackDb.txt", sep=""), sep="")
write.table(track_db, file=track_db_file_path, quote=F, row.names=F, col.names=F);

# Create the genomes.txt output.
contents <- c(paste("genome", opt$build), paste("trackDb ", opt$build, "/", "trackDb.txt", sep=""))
genomes_file_path <- paste(opt$output_trackhub_files_path, "/", "genomes.txt", sep="")
write.table(contents, file=genomes_file_path, quote=F, row.names=F, col.names=F);

# Create the hub.txt output.
contents <- c(paste("hub", opt$hub_id), paste("shortLabel", opt$hub_id), paste("longLabel", hub_name), paste("genomesFile genomes.txt", sep=""), paste("email", opt$email))
hub_file_path <- paste(opt$output_trackhub_files_path, "/", "hub.txt", sep="")
write.table(contents, file=hub_file_path, quote=F, row.names=F, col.names=F);
