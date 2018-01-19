#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--burnin_num"), action="store", dest="burnin_num", type="integer", help="Number of burnin steps"),
    make_option(c("--bychr"), action="store_true", dest="bychr", default=FALSE, help="Output chromosomes in separate files"),
    make_option(c("--hp"), action="store_true", dest="hp", default=FALSE, help="Discourage state transition across chromosomes"),
    make_option(c("--initial_states"), action="store", dest="initial_states", type="integer", default=NULL, help="Initial number of states"),
    make_option(c("--log2"), action="store", dest="log2", type="double", default=NULL, help="log2 transformation"),
    make_option(c("--maxerr"), action="store", dest="maxerr", type="double", default=NULL, help="Maximum standard deviation for the emission Gaussian distribution"),
    make_option(c("--max_cell_type_clusters"), action="store", dest="max_cell_type_clusters", type="integer", default=NULL, help="Maximum number of cell type clusters allowed"),
    make_option(c("--max_position_classes"), action="store", dest="max_position_classes", type="integer", default=NULL, help="Maximum number of position classes to be inferred"),
    make_option(c("--max_states"), action="store", dest="max_states", type="double", default=NULL, help="Maximum number of states to be inferred"),
    make_option(c("--mcmc_num"), action="store", dest="mcmc_num", type="integer", help="Number of maximization steps"),
    make_option(c("--minerr"), action="store", dest="minerr", type="double", default=NULL, help="Minimum standard deviation for the emission Gaussian distribution"),
    make_option(c("--norm"), action="store_true", dest="norm", default=FALSE, help="Standardize all datasets"),
    make_option(c("--output_log"), action="store", dest="output_log", default=NULL, help="Output log file path"),
    make_option(c("--prep_output_config"), action="store", dest="prep_output_config", help="prepMat output config file"),
    make_option(c("--prior_concentration"), action="store", dest="prior_concentration", type="double", default=NULL, help="Prior concentration"),
    make_option(c("--project_name"), action="store", dest="project_name", help="Outputs will have this base name"),
    make_option(c("--rseed"), action="store", dest="rseed", type="integer", help="Seed for IDEAS model initialization"),
    make_option(c("--save_ideas_log"), action="store", dest="save_ideas_log", default=NULL, help="Flag to save IDEAS process log"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--thread"), action="store", dest="thread", type="integer", help="Process threads"),
    make_option(c("--tmp_dir"), action="store", dest="tmp_dir", help="Directory of bed files"),
    make_option(c("--training_iterations"), action="store", dest="training_iterations", type="integer", default=NULL, help="Number of training iterations"),
    make_option(c("--training_windows"), action="store", dest="training_windows", type="integer", default=NULL, help="Number of training iterations"),
    make_option(c("--windows_bed"), action="store", dest="windows_bed", default=NULL, help="Bed file containing bed windows"),
    make_option(c("--windows_config"), action="store", dest="windows_config", default=NULL, help="Windows positions by chroms config")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

add_output_redirect <- function(cmd, save_ideas_log, output_log, default_log_name) {
    if (is.null(save_ideas_log)) {
        new_cmd = c(cmd, "&>>", default_log_name);
    }else {
        new_cmd = c(cmd, "&>>", output_log);
    }
    return(paste(new_cmd, collapse=" "));
}

combine_state <- function(parafiles, method="ward.D", mycut=0.9, pcut=1.0) {
    X = NULL;
    K = NULL;
    I = NULL;
    myheader = NULL;
    p = NULL;
    for(i in 1:length(parafiles)) {
        x = fread(parafiles[i]);
        t = max(which(is.na(x[1,])==F));
        x = as.matrix(x[,1:t]);
        if(i==1) {
            myheader = colnames(x);
            p = sqrt(9/4-2*(1-length(myheader))) - 3 / 2;
        }
        m = match(myheader[1:p+1], colnames(x)[1:p+1]);
        v = NULL;
        for(ii in 1:p) {
            for(jj in 1:ii) {
                a = max(m[ii],m[jj]);
                b = min(m[ii],m[jj]);
                v = c(v, a*(a+1)/2+b-a);
            }
        }
        X = rbind(X, array(as.matrix(x[, c(1, 1+m, 1+p+v)]), dim=c(length(x) / (1+p+length(v)), 1 + p + length(v))));
        K = c(K, dim(x)[1]);
        I = c(I, rep(i, dim(x)[1]));
    }
    N = length(parafiles);
    p = (sqrt(1 + dim(X)[2] * 8) - 3) / 2;
    omycut = mycut;
    mycut = round(length(parafiles) * mycut);
    M = array(X[,1:p+1] / X[,1], dim=c(dim(X)[1], p));
    V = array(0, dim=c(dim(X)[1] * p, p));
    for(i in 1:dim(X)[1]) {
        t = (i - 1) * p;
        l = 1;
        for(j in 1:p) {
            for(k in 1:j) {
                V[t+j, k] = V[t+k, j] = X[i,1+p+l] / X[i,1] - M[i,j] * M[i,k];
                l = l + 1;
            }
        }
        V[t+1:p,] = t(solve(chol(V[t+1:p,] + diag(1e-1,p))));
    }
    D = array(0, dim=rep(dim(X)[1],2));
    for(i in 2:dim(X)[1]) {
        for(j in 1:(i-1)) {
            D[i,j] = D[j,i] = sqrt((sum((V[(i-1)*p+1:p,]%*%(M[i,]-M[j,]))^2) + sum((V[(j-1)*p+1:p,]%*%(M[i,]-M[j,]))^2)));
        }
    }
    MM = NULL;
    kk = NULL;
    for(i in 1:N) {
        t = 1:K[i];
        if(i > 1) {
            t = t + sum(K[1:(i-1)]);
        }
        t = (1:dim(D)[1])[-t];
        h = hclust(as.dist(D[t,t]), method=method);
        k = -1;
        tM = NULL;
        for(j in min(K):(min(length(t), max(K)*2))) {
            m = cutree(h,k=j);
            tt = NULL;
            for(l in 1:j) {
                tt[l] = length(unique(I[t[which(m==l)]]));
            }
            tk = length(which(tt>=mycut));
            if(tk > k) {
                k = tk;
                tM = make_parameter(1:j, I[t], m, mycut, X[t,]);
            } else if(tk < k) {
                break;
            }
        }
        kk[i] = k;
        MM = rbind(MM, cbind(i, tM));
    }
    mysel = median(kk);
    h = hclust(as.dist(D), method=method);
    rt = rep(0, max(K)*2);
    k = -1;
    for(i in min(K):min(dim(D)[1], max(K)*2)) {
        m = cutree(h,k=i);
        tt = NULL;
        for(j in 1:i) {
            tt[j] = length(unique(I[which(m==j)]));
        }
        tk = length(which(tt>=mycut));
        if(tk==mysel | tk<k) {
            break;
        }
        k = tk;
        rt[i] = length(which(tt>=mycut));
    }
    mysel = max(k,tk);
    m = cutree(h, k=mysel);
    nn = NULL;
    for(i in 1:mysel) {
        t = which(m==i);
        nn[i] = sum(X[t,1]);
    }
    oo = order(nn, decreasing=T);
    rt = make_parameter(oo, I, m, mycut, X);
    onstate = max(rt[,1]) + 1;
    ooo = NULL;
    for(i in oo) {
        t = which(m==i);
        if(length(unique(I[t])) >= mycut) {
            ooo = c(ooo, i);
        }
    }
    d = NULL;
    for(i in 1:N) {
        d = rbind(d, compare_two(rt, MM[MM[,1]==i,-1])[1:onstate]);
    }
    dd = array(cutree(hclust(dist(c(d))), k=2), dim=dim(d));
    kk = table(c(dd));
    kk = which(as.integer(kk)==max(as.integer(kk)))[1];
    pp = apply(dd, 2, function(x){length(which(x!=kk))/length(x)});
    pp0 = apply(d, 2, function(x){length(which(x>0.5))/length(x)});
    pp[pp0<pp] = pp0[pp0<pp];
    t = which(pp > pcut);
    if(length(t) > 0) {
        j = 0;
        nrt = NULL;
        for(i in (1:onstate-1)[-t]) {
            nrt = rbind(nrt, cbind(j, rt[rt[,1]==i,-1]));
            j = j + 1;
        }
        rt = nrt;
        ooo = ooo[-t];
    }
    nrt = NULL;
    for(i in 0:max(rt[,1])) {
        t = which(rt[,1]==i);
        nrt = rbind(nrt, apply(array(rt[t,], dim=c(length(t), dim(rt)[2])), 2, sum)[-1]);
    }
    rt = nrt;
    colnames(rt) = myheader;
    O = NULL;
    Ip = NULL;
    Xp = NULL;
    k = 0;
    for(i in 1:length(parafiles)) {
        str = gsub(".para", ".profile", parafiles[i]);
        p = as.matrix(read.table(str));
        u = array(0, dim=c(dim(p)[1], length(ooo)));
        for(j in 1:length(ooo)) {
            t = which(m[k+1:K[i]] == ooo[j]);
            u[,j] = apply(array(p[,1+t], dim=c(dim(p)[1], length(t))), 1, sum);
        }
        k = k + K[i];
        u = u / (apply(u, 1, sum) + 1e-10);
        Xp = rbind(Xp, cbind(p[,1], u));
        Ip = c(Ip, rep(i,dim(u)[1]));
    }
    hp = hclust(dist(((Xp[,-1]+min(1e-3, min(Xp[,-1][Xp[,-1]>0]))))), method=method);
    ocut = min(mycut/2, length(parafiles)/2);
    t = range(as.integer(table(Ip)));
    Kp = NULL;
    for(i in t[1]:(t[2]*2)) {
        m = cutree(hp, k=i);
        tt = table(Ip,m);
        ll = apply(tt, 2, function(x){length(which(x>0))});
        Kp = c(Kp, length(which(ll>=ocut)));
    }
    oN = (t[1]:(t[2]*2))[which(Kp==max(Kp))[1]];
    m = cutree(hp, k=oN);
    tt = table(Ip,m);
    ll = apply(tt, 2, function(x){length(which(x>0))});
    tt = which(ll>=ocut);
    for(i in tt) {
        t = which(m==i);
        O = rbind(O, c(sum(Xp[t, 1]), apply(array(Xp[t,-1]*Xp[t,1], dim=c(length(t), dim(Xp)[2]-1)), 2, sum)/sum(Xp[t, 1])));
    }
    nrt = NULL;
    nrt$para = rt;
    nrt$profile = O;
    return(nrt);
}

compare_two <- function(n, m) {
    NN = get_mean(n);
    MM = get_mean(m);
    p = (-3 + sqrt(9 + 8 * (dim(n)[2] - 2))) / 2;
    dd = NULL;
    for (i in 1:dim(NN)[1]) {
        dd[i] = min(apply(array(MM[,1:p], dim=c(dim(MM)[1],p)), 1, function(x){sqrt(sum((x-NN[i,1:p])^2))}));
    }
    for (i in 1:dim(MM)[1]) {
        dd[i+dim(NN)[1]] = min(apply(array(NN[,1:p], dim=c(dim(NN)[1],p)), 1, function(x){sqrt(sum((x-MM[i,1:p])^2))}));
    }
    return(dd);
}

get_base_cmd <- function(prep_output_config, windows_bed, training_iterations, bychr, hp, norm, log2,
        max_states, initial_states, max_position_classes, max_cell_type_clusters, prior_concentration,
        burnin_num, mcmc_num, minerr, maxerr, rseed, thread) {
    base_cmd = paste("ideas", prep_output_config, sep=" ");
    if (!is.null(windows_bed)) {
        base_cmd = paste(base_cmd, windows_bed, sep=" ");
    }
    if (!is.null(training_iterations)) {
        base_cmd = paste(base_cmd, "-impute none", sep=" ");
    }
    if (bychr) {
        base_cmd = paste(base_cmd, "-bychr", sep=" ");
    }
    if (hp) {
        base_cmd = paste(base_cmd, "-hp", sep=" ");
    }
    if (norm) {
        base_cmd = paste(base_cmd, "-norm", sep=" ");
    }
    if (!is.null(log2)) {
        base_cmd = paste(base_cmd, "-log2", log2, sep=" ");
    }
    if (!is.null(max_states)) {
        base_cmd = paste(base_cmd, "-G", max_states, sep=" ");
    }
    if (!is.null(initial_states)) {
        base_cmd = paste(base_cmd, "-C", initial_states, sep=" ");
    }
    if (!is.null(max_position_classes)) {
        base_cmd = paste(base_cmd, "-P", max_position_classes, sep=" ");
    }
    if (!is.null(max_cell_type_clusters)) {
        base_cmd = paste(base_cmd, "-K", max_cell_type_clusters, sep=" ");
    }
    if (!is.null(prior_concentration)) {
        base_cmd = paste(base_cmd, "-A", prior_concentration, sep=" ");
    }
    base_cmd = paste(base_cmd, "-sample", burnin_num, mcmc_num, sep=" ");
    if (!is.null(minerr)) {
        base_cmd = paste(base_cmd, "-minerr", minerr, sep=" ");
    }
    if (!is.null(maxerr)) {
        base_cmd = paste(base_cmd, "-maxerr", maxerr, sep=" ");
    }
    base_cmd = paste(base_cmd, "-rseed", rseed, sep=" ");
    base_cmd = paste(base_cmd, "-thread", thread, sep=" ");
    return(base_cmd);
}

get_mean <- function(n) {
    N = NULL;
    for(i in sort(unique(n[,1]))) {
        t = which(n[,1]==i);
        N = rbind(N, apply(array(n[t,], dim=c(length(t), dim(n)[2])), 2, sum)[-1]);
    }
    NN = N[,-1] / N[,1];
    return(array(NN, dim=c(length(NN)/(dim(n)[2]-2), dim(n)[2]-2)));
}

get_post_training_base_cmd <- function(base_cmd, para) {
    # Change base_cmd due to training mode.
    base_cmd_items = as.list(strsplit(base_cmd[1], split=" ", fixed=TRUE))[[1]];
    if (length(which(base_cmd_items == "-G")) == 0) {
        base_cmd_items = c(base_cmd_items, "-G", length(para)-1);
    } else {
        tt = which(base_cmd_items == "-G");
        base_cmd_items[tt + 1] = length(para)-1;
    }
    tt = which(base_cmd_items == '-C');
    if(length(tt) > 0) {
        base_cmd_items = base_cmd_items[-c(tt, tt+1)];
    }
    base_cmd = paste(base_cmd_items, collapse=" ");
    return(base_cmd);
}

get_windows_by_chrom <- function(windows_config) {
    if (is.null(windows_config)) {
        windows_by_chrom = NULL;
    } else {
        fh = file(windows_config, "r");
        windows_by_chrom = readLines(fh);
        close(fh);
    }
    return(windows_by_chrom);
}

make_parameter <- function(myorder, id, mem, mycut, para) {
    rt = NULL;
    j = 0;
    for(i in myorder) {
        t = which(mem==i);
        if (length(unique(id[t])) >= mycut) {
            rt = rbind(rt, cbind(j, array(para[t,], dim=c(length(t), dim(para)[2]))));
            j = j + 1;
        }
    }
    return(rt);
}

remove_files <- function(path, pattern) {
    files = list.files(path=path, pattern=pattern);
    for (f in files) {
        unlink(f);
    }
}

run_cmd <- function(cmd, save_ideas_log, output_log, default_log_name) {
    rc = system(cmd);
    if (rc != 0) {
        if (is.null(save_ideas_log)) {
            file.rename(default_log_name, output_log);
        }
        quit(save="no", status=rc);
    }
}

default_log_name = "ideas_log.txt";
windows_by_chrom = get_windows_by_chrom(opt$windows_config);
base_cmd = get_base_cmd(opt$prep_output_config, opt$windows_bed, opt$training_iterations, opt$bychr,
        opt$hp, opt$norm, opt$log2, opt$max_states, opt$initial_states, opt$max_position_classes,
        opt$max_cell_type_clusters, opt$prior_concentration, opt$burnin_num, opt$mcmc_num, opt$minerr,
        opt$maxerr, opt$rseed, opt$thread);
output_base_name = opt$project_name;

if (is.null(opt$training_iterations)) {
    # Not performing training.
    if (is.null(windows_by_chrom)) {
        # Not performing windows by chromosome.
        output_name = output_base_name;
        cmd = paste(base_cmd, "-o", output_name, sep=" ");
        cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
    } else {
        # Performing windows by chromosome.
        for (i in 1:length(windows_by_chrom)) {
            line = windows_by_chrom[i];
            items = strsplit(line, " ")[[1]];
            chrom = items[1];
            window_start = items[2];
            window_end = items[3];
            output_name = paste(output_base_name, chrom, sep=".");
            cmd = paste(base_cmd, "-inv", window_start, window_end, sep=" ");
            cmd = paste(cmd, "-o", output_name, sep=" ");
            cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
            run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        }
    }
} else {
    # performing training.
    output_para0 = paste(output_base_name, "para0", sep=".");
    output_profile0 = paste(output_base_name, "profile0", sep=".");
    for (i in 1:opt$training_iterations) {
        cmd = paste(base_cmd, "-o", paste(output_base_name, ".tmp.", i, sep=""), sep=" ");
        cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
    }
    tpara = combine_state(paste(output_base_name, "tmp", (1:opt$training_iterations), "para", sep="."), mycut=0.5);
    write.table(tpara$profile, output_profile0, quote=F, row.names=F, col.names=F);
    para = tpara$para;
    para = apply(para, 1, function(x){paste(x, collapse=" ")});
    para = c(readLines(paste(output_base_name, "tmp", "1", "para", sep="."), n=1), para);
    writeLines(para, output_para0);
    # Now run IDEAS based on the files produced during training.
    base_cmd = get_post_training_base_cmd(base_cmd, para);
    base_cmd = paste(base_cmd, "-otherpara", output_para0[[1]], output_profile0[[1]], sep=" ");
    if (is.null(windows_by_chrom)) {
        cmd = c(base_cmd, "-o", output_base_name);
        cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
    } else {
        # Performing windows by chromosome.
        if (length(windows_by_chrom) == 1) {
            output_name = paste(output_base_name, i, sep=".");
            cmd = c(base_cmd, "-o", output_name);
            cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
            run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        } else {
            for (i in 1:length(windows_by_chrom)) {
                line = windows_by_chrom[i];
                items = strsplit(line, " ")[[1]];
                chrom = items[[1]];
                window_start = items[[2]];
                window_end = items[[3]];
                cmd = paste(base_cmd, "-inv", window_start, window_end, sep=" ");
                output_name = paste(output_base_name, chrom, sep=".");
                cmd = paste(cmd, "-o", output_name, sep=" ");
                cmd = add_output_redirect(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
                run_cmd(cmd, opt$save_ideas_log, opt$output_log, default_log_name);
            }
        }
    }
    # Remove temporary outputs.
    remove_files(path=".", pattern="tmp");
}
