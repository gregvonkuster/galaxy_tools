#!/usr/bin/env Rscript

get_methylated_read_statistics <- function(obj) {
    # Descriptive statistics for methylated reads.
    methc <- lapply(obj, function(x) {
        r <- apply(cbind(x$c1, x$c2), 1, max)
        return(r)
    })

    do.call(rbind, lapply(methc, function(x) {
        q10 <- quantile(x, 0.1)
        q60 <- quantile(x, 0.6)
        q9999 <- quantile(x, 0.9999)
        idx1 <- which(x >= q60)
        idx2 <- which(x <= 500)
        q95 <- quantile(x, 0.95)
        idx <- intersect(idx1, idx2)
        return(c(round(summary(x)),
                q10,
                q60,
                quantile(x, c(0.95, 0.99, 0.999, 0.9999)),
                '#sites_ge_8' = sum(x >= 8),
                'q60_to_500' = sum((x >= q60) & (x <= 500)),
                '#sites_gt_500' = sum(x > 500)))
        })
    )
}

output_statistics <- function(output_file, csc_df=NULL, mrs_df=NULL, critical_vals_df=NULL) {
    sink(output_file);
    cat("<html><head></head><body>");
    if (!is.null(csc_df)) {
        cat("<h3>Cytosine Site Coverage</h3>");
        print(xtable(csc_df), type="html");
    }
    if (!is.null(mrs_df)) {
        cat("<h3>Methylated Reads Statistics</h3>");
        print(xtable(mrs_df), type="html");
    }
    if (!is.null(critical_vals_df)) {
        cat("<h3>Critical Values from Empirical Cumulative Probability Distributions</h3>");
        print(xtable(critical_vals_df), type="html");
    }
    cat("</body></html>");
    sink();
}

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

