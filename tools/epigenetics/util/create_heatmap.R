#!/usr/bin/env Rscript

build_state_color_codes_vector <- function(data_matrix, histone_mark_color, color_code_type="rgb") {
    # Return  vector of color code strings for each state
    # in the received data_matrix.  The values will be either
    # rgb strings (e.g., 255,255,0) or hex code strings (e.g.,
    # #FFFFFF) depending on the value of color_code_type,
    # which can be one of "rgb" or "hex".
    range_vector = apply(data_matrix, 1, range);
    mm = NULL;
    for(i in 1:dim(data_matrix)[1]) {
        range_val1 = range_vector[1, i] + 1e-10;
        range_val2 = range_vector[2, i];
        mm = rbind(mm, (data_matrix[i,] - range_val1) / (range_val2 - range_val1));
    }
    mm = mm^5;
    if(dim(mm)[2] > 1) {
        mm = mm / (apply(mm, 1, sum) + 1e-10);
    }
    state_color = mm%*%histone_mark_color;
    s = apply(data_matrix, 1, max);
    s = (s - min(s)) / (max(s) - min(s) + 1e-10);
    state_color = round(255 - (255 - state_color) * s/0.5);
    state_color[state_color<0] = 0;
    if (identical(color_code_type, "rgb")) {
        # Here rgb_values is something like 255,255,255 217,98,0.
        state_colors_vector = paste(state_color[,1], state_color[,2], state_color[,3], sep=",");
    } else {
        # Here hex_code_strings is something like #FFFFFF #D96200
        # which is a one-to-one map to the above rgb_values.
        hex_code_strings = t(apply(state_color, 1, function(x){rgb2hsv(x[1], x[2], x[3])}));
        state_colors_vector = apply(hex_code_strings, 1, function(x){hsv(x[1], x[2], x[3])});
    }
    return(state_colors_vector);
}

create_heatmap <- function(data_frame, output_file_name, colors=c("white", "dark blue")) {
    # Plot a heatmap for a .para / .state combination based on the
    # received data_frame which was created by reading the .para file.
    state_colors_vector = get_state_color_codes_vector(data_frame, colors=colors, color_code_type="hex");
    # Open the output PDF file.
    pdf(file=output_file_name);
    # rownames(data_matrix) are the state indexes,
    # and will look something like this:
    # 0 (5.89%) 1 (91.78%) 2 (1.48%) 3 (0.86%)
    rownames(data_matrix) = paste(1:num_rows-1, " (", round(data_frame[,1]/sum(data_frame[,1])*10000)/100, "%)", sep="");
    # Set graphical parameters.
    par(mar=c(6, 1, 1, 6));
    # Create a vector containing the minimum and maximum values in data_matrix.
    min_max_vector = range(data_matrix);
    # Create a color palette.
    my_palette = colorRampPalette(colors)(n=100);
    default_palette = palette(my_palette);
    # Plot the heatmap for the current .para / .state combination.
    plot(NA, NA, xlim=c(0, p+0.7), ylim=c(0, num_rows), xaxt="n", yaxt="n", xlab=NA, ylab=NA, frame.plot=F);
    axis(1, at=1:p-0.5, labels=colnames(data_matrix), las=2);
    axis(4, at=1:num_rows-0.5, labels=rownames(data_matrix), las=2);
    col = round((t(data_matrix) - min_max_vector[1]) / (min_max_vector[2] - min_max_vector[1]) * 100);
    rect(rep(1:p-1, num_rows), rep(1:num_rows-1, each=p), rep(1:p, num_rows), rep(1:num_rows, each=p), col=col);
    rect(rep(p+0.2, num_rows), 1:num_rows-0.8, rep(p+0.8, num_rows), 1:num_rows-0.2, col=state_colors_vector);
    palette(default_palette);
    dev.off();
}

get_state_color_codes_vector <- function(data_frame, colors=c("white", "dark blue"), color_code_type="rgb") {
    # Return a vector of color strings for each row in data_frame.
    # These string will either be rgb (e.g., 255,255,0) or hex codes
    # (e.g., #FFFFFF), depending on the value of color_code_type.
    num_columns = dim(data_frame)[2];
    num_rows = dim(data_frame)[1];
    p = (sqrt(9 + 8 * (num_columns-1)) - 3) / 2;
    data_matrix = as.matrix(data_frame[,1+1:p] / data_frame[,1]);
    # colnames(data_matrix) will look something like this:
    # H3K4me3 H3K4me1 DNase H3K79me2
    colnames(data_matrix) = colnames(data_frame)[1+1:p];
    histone_marks = colnames(data_matrix);
    histone_mark_color = t(col2rgb(terrain.colors(ceiling(p))[1:p]));
    # Specify colors for common feature names like "h3k4me3".
    # These are histone marks frequently used to identify
    # promoter activities in a cell, and are often displayed
    # in shades of red.
    for(i in 1:length(histone_marks)) {
        if(regexpr("h3k4me3", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(255, 0, 0);
        }
        if(regexpr("h3k4me2", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(250, 100, 0);
        }
        if(regexpr("h3k4me1", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(250, 250, 0);
        }
        if(regexpr("h3k36me3", tolower(histone_marks[i]))>0) {
            histone_mark_color[i,] = c(0, 150, 0);
        }
        if(regexpr("h2a", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(0, 150, 150);
        }
        if(regexpr("dnase", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(0, 200, 200);
        }
        if(regexpr("h3k9ac", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(250, 0, 200);
        }
        if(regexpr("h3k9me3", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(100, 100, 100);
        }
        if(regexpr("h3k27ac", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(250, 150, 0);
        }
        if(regexpr("h3k27me3", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(0, 0, 200);
        }
        if(regexpr("h3k79me2", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(200, 0, 200);
        }
        if(regexpr("h4k20me1", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(50, 200, 50);
        }
        if(regexpr("ctcf", tolower(histone_marks[i])) > 0) {
            histone_mark_color[i,] = c(200, 0, 250);
        }
        state_colors_vector = build_state_color_codes_vector(data_matrix, histone_mark_color, color_code_type=color_code_type);
    }
    return(state_colors_vector);
}