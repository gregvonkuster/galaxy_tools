#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-b", "--build"), action="store", dest="build", help="Genome build"),
    make_option(c("-c", "--chrom_len_file"), action="store", dest="chrom_len_file", help="Chromosome length file"),
    make_option(c("-d", "--header"),  action="store", dest="header", default=NULL, help="Track header"),
    make_option(c("-e", "--state_name"),  action="store", dest="state_name", help="State name"),
    make_option(c("-f", "--hub_id"),  action="store", dest="hub_id", help="Not sure what this is"),
    make_option(c("-i", "--cell_info"),  action="store", dest="cell_info", default=NULL, help="Not sure what this is"),
    make_option(c("-n", "--hub_name"),  action="store", dest="hub_name", default=NULL, help="Not sure what this is"),
    make_option(c("-p", "--parameter_files"),  action="store", dest="parameter_files", help="Comma-separated list of IDEAS parameter files"),
    make_option(c("-s", "--state_files"),  action="store", dest="state_files", help="Comma-separated list of IDEAS state files"),
    make_option(c("-t", "--target_url"),  action="store", dest="target_url", help="target url for tracks, not sure what it is used for"),
    make_option(c("-u", "--output_track_db"),  action="store", dest="output_track_db", help="Output track db file"),
    make_option(c("-v", "--output_build"),  action="store", dest="output_build", help="Output genome build file"),
    make_option(c("-w", "--output_hub"),  action="store", dest="output_hub", help="Output hub file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

#ideas_state_files <- strsplit(opt$state_files, ",")

get_state_color<-function(statemean, markcolor=NULL)
{   
    if(length(markcolor)==0) {
        markcolor=rep("",dim(statemean)[2]);
        markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1);
        markcolor=t(col2rgb(markcolor));
    }
    rg=apply(statemean,1,range);
    mm=NULL;
    for(i in 1:dim(statemean)[1]) {
       mm=rbind(mm,(statemean[i,]-rg[1,i])/(rg[2,i]-rg[1,i]+1e-10));
    }
    mm = mm^6; 
    mm = mm / (apply(mm, 1, sum) + 1e-10);
    mycol = mm%*%markcolor;
    s = apply(statemean, 1, max);
    s = (s - min(s)) / (max(s) - min(s) + 1e-10);
    h = t(apply(mycol, 1, function(x){rgb2hsv(x[1], x[2], x[3])}));
    h[,2] = h[,2] * s;
    h = apply(h, 1, function(x){hsv(x[1], x[2], x[3])});
    rt = cbind(apply(t(col2rgb(h)), 1, function(x){paste(x, collapse=",")}), h);
    return(rt);
}

create_track<-function(state_files, chrom_len_file, outpref, state_color, header, state_name) {
    message("Reading state file: ", appendLF=FALSE);
    genomesz = read.table(chrom_len_file);
    g = NULL;
    for(i in state_files) {
        message(paste(i, " ", sep=""), appendLF=F);
        tg = as.matrix(fread(i));
        t = NULL;
        for(j in 1:dim(genomesz)[1]) {
            t = c(t, which(tg[,2]==as.character(genomesz[j, 1]) & as.numeric(tg[,4]) > as.numeric(genomesz[j, 2])));
        }   
        if(length(t) > 0) {
            tg = tg[-t,];
        }
        t = which(is.na(match(tg[,2], genomesz[,1]))==T);
        if(length(t)>0) {
            tg = tg[-t,];
        }
        print(c(dim(g),dim(tg)));
        g = rbind(g, tg);
    }
    message("Done");
    uchr = sort(unique(as.character(g[,2])));
    g1 = NULL;
    for(i in uchr){
        t = which(g[,2]==i);
        g1 = rbind(g1, g[t[order(as.integer(g[t, 3]))],]);
    }
    g = NULL;
    chr = as.character(g1[,2]);
    posst = as.numeric(g1[,3]);
    posed = as.numeric(g1[,4]);
    state = as.matrix(g1[,5:(dim(g1)[2]-1)]);
    if(length(state_name)==0) {
        state_name=0:max(state);
    }
    L = dim(g1)[1];
    n = dim(state)[2];
    if(length(header) > 0) {
        colnames(g1) = header;
    }
    cells = as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
    g1 = NULL;
    message("Generating tracks");
    options(scipen=999);
    tt = which(chr[2:L]!=chr[2:L-1]);
    tt = c(tt,which(posst[2:L]!=posed[2:L-1]));
    tt = sort(unique(tt));
    for(i in 1:n) {
        tstate = state[,i];
        t = c(tt,which(tstate[2:L]!=tstate[2:L-1]));
        t = sort(unique(t));
        t0 = c(0, t) + 1;
        t = c(t, L);
        np = cbind(chr[t], posst[t0], posed[t], tstate[t]);
        x = cbind(np[,1:3], state_name[as.integer(np[,4])+1], 1000, ".", np[,2:3], state_color[as.numeric(np[,4])+1]);
        write.table(as.matrix(x), paste(outpref, i, "bed1", sep="."), quote=F, row.names=F, col.names=F);
        print(x[1,]);
        system(paste("bedToBigBed ", outpref, ".", i, ".bed1 ", chrom_len_file, " ", outpref, ".", i, ".bb", sep=""));
        system(paste("rm ", paste(outpref, i, "bed1", sep=".")));
    }
    return(cells);
}

if(length(opt$hub_name) == 0) {
    hub_name = opt$hub_id;
}

# Create the tracks output directory.
tracks_output_dir = paste("tracks_", opt$hub_id, "/", sep="");
dir.create(tracks_output_dir, showWarnings=FALSE);

# Create the color scheme.
mc = NULL;
x = read.table(opt$parameter_files, comment="!", nrows=1);
l = as.integer(regexpr("\\*", as.matrix(x)));
l = min(which(l>0))-2;
x = as.matrix(read.table(opt$parameter_files));
if(length(opt$parameter_files) > 1) {
    for(i in 2:length(opt$parameter_files)) {
        x=x+as.matrix(read.table(paste(opt$parameter_files[i],".para",sep="")));
    }
}
p = x[,1] / sum(x[,1]);
m = array(as.matrix(x[,1:l+1] / x[,1]), dim=c(dim(x)[1], l));
state_color = get_state_color(m, mc);

# Create the tracks.
cells = create_track(opt$state_files, opt$chrom_len_file, paste(tracks_output_dir, opt$hub_id, sep=""), state_color, header=opt$header, state_name=opt$state_name);
cell_info = opt$cell_info
if(length(cell_info) == 0) {
    cell_info = cbind(cells, cells, cells, "#000000");
    cell_info = array(cell_info, dim=c(length(cells), 4));
}
cell_info = as.matrix(cell_info);
track_db = NULL;

for(i in 1:length(cells)) {
    ii = which(cells[i] == cell_info[,1]);
    if(length(ii) == 0) {
        next;
    }
    ii = ii[1];
    track_db = c(track_db, paste("track bigBed", i, sep=""));
    track_db = c(track_db, paste("priority", ii));
    track_db = c(track_db, "type bigBed 9 .");
    track_db = c(track_db, "itemRgb on");
    track_db = c(track_db, "maxItems 100000");
    track_db = c(track_db, paste("bigDataUrl ", opt$target_url, opt$hub_id, ".", i, ".bb", sep=""));
    track_db = c(track_db, paste("shortLabel", cell_info[ii, 2]));
    track_db = c(track_db, paste("longLabel", paste(opt$hub_name, cell_info[ii, 3])));
    track_db = c(track_db, paste("color", paste(c(col2rgb(cell_info[ii, 4])), collapse=",")));
    track_db = c(track_db, "visibility dense");
    track_db = c(track_db, ""); 
}

# Write the outputs.
write.table(track_db, opt$output_track_db, quote=F, row.names=F, col.names=F);
write.table(c(paste("genome", opt$build), opt$output_build), paste(tracks_output_dir, "genomes_", opt$hub_id, ".txt", sep=""), quote=F, row.names=F, col.names=F);
write.table(c(paste("hub", opt$hub_id), paste("shortLabel", opt$hub_id), paste("longLabel", opt$hub_name), paste("genomesFile genomes_", opt$hub_id, ".txt", sep=""), "email yzz2 at psu.edu"), opt$output_hub, quote=F, row.names=F, col.names=F);
