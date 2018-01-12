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
    make_option(c("--window_end"), action="store", dest="window_end", type="integer", default=NULL, help="Windows positions by chromosome end value"),
    make_option(c("--window_start"), action="store", dest="window_start", type="integer", default=NULL, help="Windows positions by chromosome start value")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

add_output_redirect <- function(cmd, save_ideas_log, output_log, default_log_name) {
    if (is.null(save_ideas_log)) {
        cmd = paste(cmd, "&>>", default_log_name, sep=" ");
    }else {
        cmd = paste(cmd, "&>>", output_log, sep=" ");
    }
    return(cmd);
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

get_mean <- function(n) {
    N = NULL;
    for(i in sort(unique(n[,1]))) {
        t = which(n[,1]==i);
        N = rbind(N, apply(array(n[t,], dim=c(length(t), dim(n)[2])), 2, sum)[-1]);
    }
    NN = N[,-1] / N[,1];
    return(array(NN, dim=c(length(NN)/(dim(n)[2]-2), dim(n)[2]-2)));
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

run_cmd <- function(cmd, save_ideas_log, output_log, default_log_name) {
    cat("\n\n >>>>> cmd:\n", cmd, "\n\n");
    rc = system(cmd);
    if (rc != 0) {
        if (is.null(save_ideas_log)) {
            file.rename(default_log_name, output_log);
        }
        quit(rc);
    }
}

default_log_name = "ideas_log.txt";
output_base_name = opt$project_name;
cmd = paste("ideas", opt$prep_output_config, sep=" ");
if (!is.null(opt$windows_bed)) {
    cmd = paste(cmd, opt$windows_bed, sep=" ");
}
if (!is.null(opt$training_iterations)) {
    cmd = paste(cmd, "-impute none", sep=" ");
}
if (opt$bychr) {
    cmd = paste(cmd, "-bychr", sep=" ");
}
if (opt$hp) {
    cmd = paste(cmd, "-hp", sep=" ");
}
if (opt$norm) {
    cmd = paste(cmd, "-norm", sep=" ");
}
if (!is.null(opt$window_start) && !is.null(opt$window_end)) {
    cmd = paste(cmd, "-inv", opt$window_start, opt$window_end, sep=" ");
}
if (!is.null(opt$log2)) {
    cmd = paste(cmd, "-log2", opt$log2, sep=" ");
}
if (!is.null(opt$max_states)) {
    cmd = paste(cmd, "-G", opt$max_states, sep=" ");
}
if (!is.null(opt$initial_states)) {
    cmd = paste(cmd, "-C", opt$initial_states, sep=" ");
}
if (!is.null(opt$max_position_classes)) {
    cmd = paste(cmd, "-P", opt$max_position_classes, sep=" ");
}
if (!is.null(opt$max_cell_type_clusters)) {
    cmd = paste(cmd, "-K", opt$max_cell_type_clusters, sep=" ");
}
if (!is.null(opt$prior_concentration)) {
    cmd = paste(cmd, "-A", opt$prior_concentration, sep=" ");
}
cmd = paste(cmd, "-sample", opt$burnin_num, opt$mcmc_num, sep=" ");
if (!is.null(opt$minerr)) {
    cmd = paste(cmd, "-minerr", opt$minerr, sep=" ");
}
if (!is.null(opt$maxerr)) {
    cmd = paste(cmd, "-maxerr", opt$maxerr, sep=" ");
}
cmd = paste(cmd, "-rseed", opt$rseed, sep=" ");
cmd = paste(cmd, "-thread", opt$thread, sep=" ");

if (is.null(opt$training_iterations)) {
    final_cmd = paste(cmd, "-o", output_base_name, sep=" ");
    final_cmd = add_output_redirect(final_cmd, opt$save_ideas_log, opt$output_log, default_log_name);
    run_cmd(final_cmd, opt$save_ideas_log, opt$output_log, default_log_name);
} else {
    output_para0 = paste(output_base_name, ".para0", sep="");
    output_profile0 = paste(output_base_name, ".profile0", sep="");
    for (i in 1:opt$training_iterations) {
        final_cmd = paste(cmd, "-o", paste(output_base_name, ".tmp.", i, sep=""), sep=" ");
        final_cmd = add_output_redirect(final_cmd, opt$save_ideas_log, opt$output_log, default_log_name);
        run_cmd(final_cmd, opt$save_ideas_log, opt$output_log, default_log_name);
    }
    tpara = combine_state(paste(output_base_name, ".tmp.", (1:opt$training_iterations), ".para", sep=""), mycut=0.5);
    para = tpara$para;
    write.table(tpara$profile, output_profile0, quote=F, row.names=F, col.names=F);
    para = apply(para, 1, function(x){paste(x, collapse=" ")});
    para = c(readLines(paste(output_base_name, ".tmp.1.para", sep=""), n=1), para);
    writeLines(para, output_para0);
    cmd = c(cmd, "-otherpara", output_para0, output_profile0);
    if (length(which(cmd == "-G")) == 0) {
        cmd = c(cmd, "-G", length(para)-1);
    } else {
        tt = which(cmd == "-G");
        cmd[tt + 1] = length(para)-1;
    }
    tt = which(cmd == '-C');
    if(length(tt)>0) {
        cmd = cmd[-c(tt, tt+1)];
    }
}
