#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
            make_option(c("-i", "--input_dir"), action="store", dest="input_dir", help="IDEAS para files directory"),
            make_option(c("-o", "--output_dir"), action="store", dest="output_dir", help="PDF output directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

create_heatmap<-function(data_frame, output_file_name) {
    # Plot a heatmap for a .para / .state combination
    # based on the received data_frame which was created
    # by reading the .para file.
    num_columns = dim(data_frame)[2];
    num_rows = dim(data_frame)[1];
    p = (sqrt(9 + 8 * (num_columns - 1)) - 3) / 2;
    data_matrix = as.matrix(data_frame[,1+1:p] / data_frame[,1]);
    colnames(data_matrix) = colnames(data_frame)[1+1:p];
    marks = colnames(data_matrix);
    rownames(data_matrix) = paste(1:num_rows, " (", round(data_frame[,1]/sum(data_frame[,1])*10000)/100, "%)", sep="");

    # Open the output PDF file.
    pdf(file=output_file_name);
    # Set graphical parameters.
    par(mar=c(6, 1, 1, 6));
    # Create a vector containing the minimum and maximum values in data_matrix.
    min_max_vector = range(data_matrix);
    # Define colors for the palette.
    colors = 0:100 / 100 * (min_max_vector[2] - min_max_vector[1]) + min_max_vector[1];
    # Create the color palette.
    my_palette = colorRampPalette(c("white", "dark blue"))(n=100);
    defpalette = palette(my_palette);
    # Plot the heatmap for the current .para / .state combination.
    plot(NA, NA, xlim=c(0, p+0.7), ylim=c(0, num_rows), xaxt="n", yaxt="n", xlab=NA, ylab=NA, frame.plot=F);
    axis(1, at=1:p-0.5, labels=colnames(data_matrix), las=2);
    axis(4, at=1:num_rows-0.5, labels=rownames(data_matrix), las=2);
    rect(rep(1:p-1, num_rows), rep(1:num_rows, each=p), rep(1:p, num_rows), rep(1:num_rows, each=p),
            col=round((t(data_matrix)-min_max_vector[1])/(min_max_vector[2]-min_max_vector[1])*100));
    markcolor = t(col2rgb(terrain.colors(ceiling(p))[1:p]));

    for(i in 1:length(marks)) {
        if(regexpr("h3k4me3", tolower(marks[i])) > 0) {
            markcolor[i,] = c(255, 0, 0);
        }
        if(regexpr("h3k4me2", tolower(marks[i])) > 0) {
            markcolor[i,] = c(250, 100, 0);
        }
        if(regexpr("h3k4me1", tolower(marks[i])) > 0) {
            markcolor[i,] = c(250, 250, 0);
        }
        if(regexpr("h3k36me3", tolower(marks[i]))>0) {
            markcolor[i,] = c(0, 150, 0);
        }
        if(regexpr("h2a", tolower(marks[i])) > 0) {
            markcolor[i,] = c(0, 150, 150);
        }
        if(regexpr("dnase", tolower(marks[i])) > 0) {
            markcolor[i,] = c(0, 200, 200);
        }
        if(regexpr("h3k9ac", tolower(marks[i])) > 0) {
            markcolor[i,] = c(250, 0, 200);
        }
        if(regexpr("h3k9me3", tolower(marks[i])) > 0) {
            markcolor[i,] = c(100, 100, 100);
        }
        if(regexpr("h3k27ac", tolower(marks[i])) > 0) {
            markcolor[i,] = c(250, 150, 0);
        }
        if(regexpr("h3k27me3", tolower(marks[i])) > 0) {
            markcolor[i,] = c(0, 0, 200);
        }
        if(regexpr("h3k79me2", tolower(marks[i])) > 0) {
            markcolor[i,] = c(200, 0, 200);
        }
        if(regexpr("h4k20me1", tolower(marks[i])) > 0) {
            markcolor[i,] = c(50, 200, 50);
        }
        if(regexpr("ctcf", tolower(marks[i])) > 0) {
            markcolor[i,] = c(200, 0, 250);
        }
        state_color = get_state_color(data_matrix, markcolor)[,2];
    }
    rect(rep(p+0.2,num_rows), 1:num_rows-0.8, rep(p+0.8,num_rows), 1:num_rows-0.2, col=state_color);
    palette(defpalette);
    dev.off();
}

get_state_color<-function(statemean, markcolor) {
    rg=apply(statemean, 1, range);
    mm=NULL;
    for(i in 1:dim(statemean)[1]) {
        mm = rbind(mm, (statemean[i,]-rg[1,i]+1e-10)/(rg[2,i]-rg[1,i]+1e-10));
    }
    mm = mm^5;
    if(dim(mm)[2] > 1) {
        mm = mm / (apply(mm, 1, sum) + 1e-10);
    }
    mycol = mm%*%markcolor;
    s = apply(statemean, 1, max);
    s = (s-min(s)) / (max(s)-min(s) + 1e-10);
    mycol = round(255-(255-mycol)*s/0.5);
    mycol[mycol<0] = 0;
    rt = paste(mycol[,1], mycol[,2], mycol[,3], sep=",");
    h = t(apply(mycol, 1, function(x){rgb2hsv(x[1], x[2], x[3])}));
    h = apply(h, 1, function(x){hsv(x[1], x[2], x[3])});
    rt = cbind(rt, h);
    return(rt);
}

# Read the inputs.
para_files <- list.files(path=opt$input_dir, pattern="\\.para$", full.names=TRUE);
for (i in 1:length(para_files)) {
    para_file <- para_files[i];
    para_file_base_name <- strsplit(para_file, split="/")[[1]][2]
    output_file_name <- gsub(".para", ".pdf", para_file_base_name)
    output_file_path <- paste(opt$output_dir, output_file_name, sep="/");
    data_frame <- read.table(para_file, comment="!", header=T);
    create_heatmap(data_frame, output_file_path);
}
