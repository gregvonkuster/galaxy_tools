#!/usr/bin/env Rscript

string_to_boolean <- function(s, default=FALSE) {
    # Convert s to boolean, setting NULL to the default.
    if (is.null(s)) {
        return (default);
    } else if (s=='yes') {
        return (TRUE);
    } else {
        return (FALSE);
    }
}

string_to_character_vector <- function(s) {
    # Convert comma-separated string s to a character vector.
    i <- grep(",", s);
    if (length(i) > 0) {
        s_list <- strsplit(s, ",")[[1]];
        return (c(unlist(s_list, use.names=FALSE)));
    } else {
        return (c(s));
    }
}

string_to_integer_or_vector <- function(s) {
    s_list <- strsplit(s, ",")[[1]];
    num_s <- length(s_list)[[1]];
    if (num_s == 1) {
        return (s_list[[1]]);
    } else {
        v <- integer(num_s);
        for (i in 1:num_s) {
            v[[i]] <- s_list[[i]];
        }
        return (v);
    }
}

zero_to_null <- function(n) {
    if (n==0) {
        return (NULL);
    } else {
        return (n);
    }
}

