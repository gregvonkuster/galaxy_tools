#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("viridisLite"))

option_list <- list(
        make_option(c("--build"), action="store", dest="build", help="Genome build"),
        make_option(c("--chrom_len_file"), action="store", dest="chrom_len_file", help="Chromosome length file"),
        make_option(c("--email"),  action="store", dest="email", help="User email address"),
        make_option(c("--galaxy_url"),  action="store", dest="galaxy_url", help="Galaxy instance base URL"),
        make_option(c("--hub_name"),  action="store", dest="hub_name", default=NULL, help="Hub name without spaces"),
        make_option(c("--input_dir_para"), action="store", dest="input_dir_para", help="Directory containing .para outputs from IDEAS"),
        make_option(c("--input_dir_state"), action="store", dest="input_dir_state", help="Directory containing .state outputs from IDEAS"),
        make_option(c("--long_label"), action="store", dest="long_label", help="Hub long label"),
        make_option(c("--output_trackhub"),  action="store", dest="output_trackhub", help="Output hub file"),
        make_option(c("--output_trackhub_files_path"),  action="store", dest="output_trackhub_files_path", help="Output hub extra files path"),
        make_option(c("--output_trackhub_id"),  action="store", dest="output_trackhub_id", help="Encoded output_trackhub dataset id"),
        make_option(c("--short_label"), action="store", dest="short_label", help="Hub short label")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

create_primary_html = function(output_trackhub, tracks_dir, build) {
    track_files <- list.files(path=tracks_dir);
    s <- paste('<html><head></head><body>', sep="\n");
    s <- paste(s, '<h3>Contents of directory ~/myHub/', build, ' required by UCSC TrackHub</h3>\n', sep="");
    s <- paste(s, '<ul>\n', sep="")
    for (i in 1:length(track_files)) {
        s <- paste(s, '<li><a href="', 'myHub/', build, "/", track_files[i], '">', track_files[i], '</a></li>\n', sep="");
    }
    s <- paste(s, '</ul>\n</body>\n</html>', sep="");
    cat(s, file=output_trackhub);
}

create_track = function(input_dir_para, input_dir_state, chrom_len_file, base_track_file_name) {
    para_files <- list.files(path=input_dir_para, full.names=TRUE);
    cat("para_files: ", para_files, "\n");
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
    state_name = 0:max(state);
    L = dim(g1)[1];
    n = dim(state)[2];
    cells = as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
    g1 = NULL;
    options(scipen=999);
    tt = which(chr[2:L]!=chr[2:L-1]);
    tt = c(tt, which(posst[2:L]!=posed[2:L-1]));
    tt = sort(unique(tt));
    state_color <- get_state_color(para_files[1])
    for(i in 1:n) {
        cat("i: ", i, "\n");
        tstate = state[,i];
        t = c(tt, which(tstate[2:L]!=tstate[2:L-1]));
        t = sort(unique(t));
        t0 = c(0, t) + 1;
        t = c(t, L);
        np = cbind(chr[t], posst[t0], posed[t], tstate[t]);
        track_file_name_bed <- get_track_file_name(base_track_file_name, i, "bed");
        track_file_name_bigbed <- get_track_file_name(base_track_file_name, i, "bigbed");
        x = cbind(np[, 1:3], state_name[as.integer(np[,4])+1], 1000, ".", np[,2:3], state_color[as.numeric(np[,4])+1]);
        write.table(as.matrix(x), track_file_name_bed, quote=F, row.names=F, col.names=F);
        system(paste("bedToBigBed ", track_file_name_bed, chrom_len_file, " ", track_file_name_bigbed));
        system(paste("rm ", track_file_name_bed));
    }
    return(cells);
}

create_track_db = function(galaxy_url, encoded_dataset_id, input_dir_para, input_dir_state, build, chrom_len_file, tracks_dir, hub_name, short_label, long_label) {
    base_track_file_name <- paste(tracks_dir, hub_name, sep="");
    cells = create_track(input_dir_para, input_dir_state, chrom_len_file, base_track_file_name);
    cell_info = cbind(cells, cells, cells, "#000000");
    cell_info = array(cell_info, dim=c(length(cells), 4));
    cell_info = as.matrix(cell_info);
    track_db = NULL;
    for (i in 1:length(cells)) {
        big_data_url <- get_big_data_url(galaxy_url, encoded_dataset_id, tracks_dir, i, build);
        ii = which(cells[i] == cell_info[,1]);
        if (length(ii) == 0) {
            next;
        }
        ii = ii[1];
        # trackDb.txt track entry.
        track_db = c(track_db, paste("hub ", hub_name, "_track_", i, sep=""));
        track_db = c(track_db, "type bigBed");
        track_db = c(track_db, paste("bigDataUrl", big_data_url, sep=" "));
        track_db = c(track_db, paste("shortLabel", short_label, sep=" "));
        track_db = c(track_db, paste("longLabel", long_label, sep=" "));
        track_db = c(track_db, paste("priority", ii));
        track_db = c(track_db, "itemRgb on");
        track_db = c(track_db, "maxItems 100000");
        track_db = c(track_db, paste("color", paste(c(col2rgb(cell_info[ii,4])), collapse=","), sep=" "));
        track_db = c(track_db, "visibility dense");
        track_db = c(track_db, "");
    }
    return(track_db);
}

get_big_data_url = function(galaxy_url, encoded_dataset_id, tracks_dir, index, build) {
    track_files <- list.files(path=tracks_dir, pattern="\\.bigbed");
    s <- paste(galaxy_url, 'datasets/', encoded_dataset_id, '/display/myHub/', build, '/', track_files[index], sep="");
    return(s)
}

get_rgb<-function(statemean, markcolor=NULL)
{
    if(length(markcolor) == 0) {
        markcolor = rep("",dim(statemean)[2]);
        markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1)
        markcolor = t(col2rgb(markcolor));
    }

    rg = apply(statemean, 1, range);
    mm = NULL;
    for(i in 1:dim(statemean)[1]) {
        mm = rbind(mm, (statemean[i,]-rg[1, i]+1e-10)/(rg[2, i]-rg[1, i]+1e-10));
    }
    mm = mm^5;
    if (dim(mm)[2]>1) {
        mm = mm / (apply(mm, 1, sum)+1e-10);
    }
    mycol = mm%*%markcolor;
    s = apply(statemean, 1, max);
    s = (s-min(s))/(max(s)-min(s)+1e-10);
    mycol = round(255-(255-mycol)*s/0.5);
    mycol[mycol<0] = 0;
    rt = paste(mycol[,1], mycol[,2], mycol[,3], sep=",");
    h = t(apply(mycol, 1, function(x){rgb2hsv(x[1], x[2], x[3])}));
    h = apply(h, 1, function(x){hsv(x[1], x[2], x[3])});
    rt = cbind(rt, h);
    return(rt);

    h = t(apply(mycol, 1, function(x){rgb2hsv(x[1], x[2], x[3])}));
    h[,2] = h[,2]*s;
    h = apply(h, 1, function(x){hsv(x[1], x[2], x[3])});
    rt = cbind(apply(t(col2rgb(h)), 1, function(x){paste(x, collapse=",")}), h);
    return(rt);
}

get_state_color <- function(para_file, cols=c("white", "dark blue")) {
    x = read.table(para_file, comment="!", header=T);
    k = dim(x)[2];
    l = dim(x)[1];
    p = (sqrt(9+8*(k-1))-3)/2;
    m = as.matrix(x[,1+1:p]/x[,1]);
    colnames(m) = colnames(x)[1+1:p];
    marks = colnames(m);
    rg = range(m);
    colors = 0:100/100*(rg[2]-rg[1])+rg[1];
    markcolor = t(col2rgb(terrain.colors(ceiling(p))[1:p]));
    for(i in 1:length(marks)) {
        if (regexpr("h3k4me3",tolower(marks[i]))>0) {
            markcolor[i,]=c(255,0,0);
        }
        if (regexpr("h3k4me2",tolower(marks[i]))>0) {
            markcolor[i,]=c(250,100,0);
        }
        if (regexpr("h3k4me1",tolower(marks[i]))>0) {
            markcolor[i,]=c(250,250,0);
        }
        if (regexpr("h3k36me3",tolower(marks[i]))>0) {
            markcolor[i,]=c(0,150,0);
        }
        if (regexpr("h2a",tolower(marks[i]))>0) {
            markcolor[i,]=c(0,150,150);
        }
        if (regexpr("dnase",tolower(marks[i]))>0) {
            markcolor[i,]=c(0,200,200);
        }
        if (regexpr("atac",tolower(marks[i]))>0) {
            markcolor[i,]=c(0,200,200);
        }
        if (regexpr("h3k9ac",tolower(marks[i]))>0) {
            markcolor[i,]=c(250,0,200);
        }
        if (regexpr("h3k9me3",tolower(marks[i]))>0) {
            markcolor[i,]=c(100,100,100);
        }
        if (regexpr("h3k27ac",tolower(marks[i]))>0) {
            markcolor[i,]=c(250,150,0);
        }
        if (regexpr("h3k27me3",tolower(marks[i]))>0) {
            markcolor[i,]=c(0,0,225);
        }
        if (regexpr("h3k79me2",tolower(marks[i]))>0) {
            markcolor[i,]=c(200,0,200);
        }
        if (regexpr("h4k20me1",tolower(marks[i]))>0) {
            markcolor[i,]=c(50,200,50);
        }
        if (regexpr("ctcf",tolower(marks[i]))>0) {
            markcolor[i,]=c(200,0,250);
        }
    }
    statecolor = get_rgb(m, markcolor)[,];
    return(statecolor);
}

get_track_file_name = function(base_track_file_name, index, ext) {
    track_file_name <- paste(base_track_file_name, index, ext, sep=".");
    return(track_file_name);
}

# Create the hub.txt output.
hub_name_line <- paste("hub ", opt$hub_name, sep="");
short_label_line <- paste("shortLabel ", opt$short_label, sep="");
long_label_line <- paste("longLabel ", opt$long_label, sep="");
genomes_txt_line <- paste("genomesFile genomes.txt", sep="");
email_line <- paste("email ", opt$email, sep="");
contents <- paste(hub_name_line, short_label_line, long_label_line, genomes_txt_line, email_line, sep="\n");
hub_dir <- paste(opt$output_trackhub_files_path, "/", "myHub", "/", sep="");
dir.create(hub_dir, showWarnings=FALSE);
hub_file_path <- paste(hub_dir, "hub.txt", sep="");
write.table(contents, file=hub_file_path, quote=F, row.names=F, col.names=F);

# Create the genomes.txt output.
genome_line <- paste("genome ", opt$build, sep="");
track_db_line <- paste("trackDb ", opt$build, "/", "trackDb.txt", sep="");
contents <- paste(genome_line, track_db_line, sep="\n");
genomes_file_path <- paste(hub_dir, "genomes.txt", sep="");
write.table(contents, file=genomes_file_path, quote=F, row.names=F, col.names=F);

# Create the tracks.
tracks_dir <- paste(hub_dir, opt$build, "/", sep="");
dir.create(tracks_dir, showWarnings=FALSE);
track_db <- create_track_db(opt$galaxy_url, opt$output_trackhub_id, opt$input_dir_para, opt$input_dir_state, opt$build, opt$chrom_len_file, tracks_dir, opt$hub_name, opt$short_label, opt$long_label);

# Create the trackDb.txt output.
track_db_file_path <- paste(tracks_dir, "trackDb.txt", sep="");
write.table(track_db, file=track_db_file_path, quote=F, row.names=F, col.names=F);

# Create the primary HTML dataset.
create_primary_html(opt$output_trackhub, tracks_dir, opt$build);