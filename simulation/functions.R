seqgen_path <- "~/apps/Seq-Gen/source/seq-gen"
mrb_path <- "./mb"

sim_alignment <- function(tree, a, b, n_site, al, dg=0, readerr=0) {
  

    seq_names <- tree$tip.label
    # init and normalize 1 matrix
    qmat <- matrix(c(-a,a,b,-b),nrow=2,byrow=TRUE)
    freq <- c(b,a) / (a+b)
    qmat <- qmat/-sum(freq*diag(qmat))
    
    # pre-mult by pi:
    qmat <- qmat%*%diag(1/freq)
    
    # set up the freqs and rates to pass to simSeq
    freq_sim <- c(freq,0,0)
    rate_sim <- c(qmat[1,2],100,100, 
    
                  100,100,
                                100 )
   
    # sim the alleles separately and combine
    if (dg == 0) {
        r <- rgamma(n_site,al,al)
        if (al == 0) r <- rep(1,n_site)

        d1 <- as.matrix(sapply(1:n_site, function(i){
            simSeq(tree, l=1, type="DNA", bf=freq_sim, Q=rate_sim, rate = r[i])
        }))
        d2 <- sapply(1:n_site, function(i){
            simSeq(tree, l=1, type="DNA", bf=freq_sim, Q=rate_sim, rate = r[i])
        })

        d1 <- matrix(unlist(d1), ncol=n_site, byrow=FALSE)
        d2 <- matrix(unlist(d2), ncol=n_site, byrow=FALSE)
        d1 <- d1 - 1
        d2 <- d2 - 1
     } else if (dg > 0) {
        r <- phangorn::discrete.gamma(al, dg)        
        if (al == 0) r <- rep(1, dg) 
        s1 <- lapply(r, function(ri) simSeq(tree, l=n_site/dg, type="DNA", bf=freq_sim, Q=rate_sim, rate = ri))
        s2 <- lapply(r, function(ri) simSeq(tree, l=n_site/dg, type="DNA", bf=freq_sim, Q=rate_sim, rate = ri))
        d1 <- lapply(s1, as.character)
        d2 <- lapply(s2, as.character)
        d1 <- matrix(unlist(d1), byrow=FALSE, ncol=n_site, nrow=length(tree$tip.label))
        d2 <- matrix(unlist(d2), byrow=FALSE, ncol=n_site, nrow=length(tree$tip.label))
        d1 <- ifelse(d1 == "a", 0, 1)
        d2 <- ifelse(d2 == "a", 0, 1)    
    }
    sim_aln <- d1 + d2
    # table(sim_aln)/sum(table(sim_aln))
    # c(b^2, 2*a*b, a^2)
   
    # add readerrs 
    if (readerr > 0) {
        nerrs <- 0
        # make the readerrs
        for (i in 1:nrow(sim_aln)) for (j in 1:ncol(sim_aln)) {
                if (runif(1) > 1-2*readerr) {
                        nerrs <- nerrs + 1
                        sim_aln[i,j] <- sample((0:2)[-(sim_aln[i,j]+1)],1)
                }
        }
        print(sprintf("number of readerrs: %s",nerrs))
    }

    seq_strings <-  apply(sim_aln, 1, function(x) paste0(c("E", "M", "D")[x+1], collapse=""))
    seq_lines <- paste(seq_names, seq_strings, sep="    ")
    seq_lines
}


sim_dna <- function(tree,n_site,rate,
                     tree_file="sim.tree",aln_file="aln.txt", 
                     al=1, freq=c(1,1,1,1)/4, 
                     missing_rate = 0) {

    tree_string <- TreeTools::NewickTree(tree)   

    # store tree w/ branch lengths to file
    writeLines(tree_string, con=tree_file)

    # simulate some data using seqgen  
    sg_call <- ifelse(al > 0, 
                 paste(seqgen_path, " -mGTR ", 
                    " -f ", paste0(freq,collapse=" "),  
                    " -r ", paste0(rate,collapse=" "), 
                    " -l ", n_site, 
                    " -a ", al,
                    " -z ", floor(runif(1)*3928109+2817615), 
                    " < ", tree_file,
                    " > ", aln_file, sep=""),
                 paste(seqgen_path, " -mGTR ", 
                    " -f ", paste0(freq,collapse=" "),  
                    " -r ", paste0(rate,collapse=" "), 
                    " -l ", n_site, 
                    " -z ", floor(runif(1)*3928109+2817615), 
                    " < ", tree_file,
                    " > ", aln_file, sep=""))
    system(sg_call)
    
    alignment_char <- phybase::read.dna.seq(aln_file,"phylip")$seq
    alignment <- apply(alignment_char, 2, function(x) 
                       match(x, c("A","C","G","T")))

    if (missing_rate > 0) {
        tot_chars <- nrow(alignment) * ncol(alignment)
        del_sites <- sample(1:tot_chars, round(tot_chars * missing_rate,0))
        alignment[del_sites] <- "-"
        alignment_char[del_sites] <- "-"
    }

    alignment_str <- apply(alignment_char, 1, 
                       function(x) paste0(x,collapse=""))

    return(list(aln=alignment, aln_str=alignment_str))

}
tip_labels <- c(letters,paste0(letters,letters))

# TODO:
#    1. make this function write rates, freqs, tree, to file(s)
#    
#    2. New function that will:
#        2a. add mb block w fixed rate prior to nexus file with alignment characters    
#        2b. fixed rates estimated by quasi (within R I guess)

# simulate an alignment, given alpha, n_taxa, and n_sites
# assumes all rate params = 1, all freqs = 0.25, and all branch lengths = 0.1.
gen_tree <- function(n_taxa, bl = 0.1) {

        tree <- ape::rtree(n_taxa)
        tree$edge.length <- rexp(length(tree$edge.length), 1/bl)
        tree$edge.length <- ifelse(tree$edge.length > 0.85, 0.85, tree$edge.length) 
        return(tree)
}

setup_mrb <- function(aln_str,
                      n_taxa, n_site,
                      model_lines,
                      mcmc_lines,
                      out_dir="sim",
                      fname="",
                      datatype="dimethyl",
                      interleave="no"
                      ) {

    ## now write the alignment:
    writeLines(c("#NEXUS\n",

    # write the data block
                "begin data;", 
                sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
                sprintf("format datatype=%s interleave=%s gap=- missing=?;",datatype, interleave),
                "matrix",
                aln_str,
                "\t;\nend;\n",

    # write the mrbayes block
                 "begin mrbayes;",
                 model_lines,
                 mcmc_lines,
                 sprintf("sump outputname=%s/sump.%s burninfrac=%s;",out_dir,fname,burninfrac),
                 sprintf("sumt outputname=%s/sumt.%s burninfrac=%s conformat=simple;",out_dir,fname,burninfrac),

                 # sprintf("comparetree outputname=treecomp"),
                 #         "filename1=tree.summ.con.tre filename2=sim.tree.nex;",
                 "end;"

                 ),
                sprintf("%s/%s.nex",out_dir,fname))
}

run_mrb <- function(mrb_path, nex_file, out_file) {
        system(sprintf("%s %s &> %s", mrb_path, nex_file,  out_file), intern=FALSE)
}

parse_mb_output <- function(ntax, nsite, lkhd, sim, out_dir) {

    # first, parse params
    parms_fname <- sprintf("%s/parms.summ.%s.%s.%s.pstat", out_dir, nsite, lkhd,sim);

    # second, summarize trees    
    contree_fname <- sprintf("%s/tree.summ.%s.%s.%s", out_dir,nsite,lkhd,sim)

    # output to get the running time
    out_fname <- sprintf()
}

compare_contree <- function(tr_true, tr_est, bl_est) {

    if (ape::is.rooted(tr_true)) tr_true <- unroot(tr_true)
    comp <- comparePhylo(tr_est, tr_true, use.edge.length=TRUE)
    splits_in_common <- as.numeric(stringr::str_extract(comp$message[7], "[0-9]+"))
    total_splits <- length(tr_est$edge.length) - length(tr_est$tip.label)

    all_equal <- all.equal.phylo(unroot(tr_true), tr_est, use.edge.length=FALSE, use.tip.label=TRUE)

    t1_bls_tr  <- reorder(tr_true, "postorder")$edge.length 
    t1_bls_tr <- tr_true$edge.length
    n_bls <- length(t1_bls_tr)
    t1_bls_est <- bl_est$Mean[1:n_bls]
    names(t1_bls_tr) <- paste0("bl_", 1:n_bls)
    names(t1_bls_est) <- paste0("bl_", 1:n_bls, "_est")

    # con_bls <- tr_est$edge.length
    t1_bls_lower <- bl_est$CredInt_Lower[1:n_bls]
    names(t1_bls_lower) <- paste0("bl_", 1:n_bls, "_est_lb") 

    t1_bls_upper <- bl_est$CredInt_Upper[1:n_bls]
    names(t1_bls_upper) <- paste0("bl_", 1:n_bls, "_est_ub")

    res <- data.frame(splits = total_splits, splits_in_common = splits_in_common)
    cbind(res, rbind(t1_bls_tr), rbind(t1_bls_est), rbind(t1_bls_lower), rbind(t1_bls_upper))
}

get_parm_ests <- function(sump_fname) {
    parm_ests <- read.csv(sump_fname, skip=1, head=TRUE, sep="\t")
    return(parm_ests)
}

# compute time
get_compute_time <- function(out_fname) {
    # get compute time 
    l <- readLines(out_fname)
    time_line <- l[grepl("Analysis used", l)]
    if (length(time_line) == 0) {
            print(out_fname)
            return(NA)
    }
    #min_sec <- try(as.numeric(stringr::str_extract_all(time_line, "[1-9\\.]+")[[1]]))
    min_sec <- try(regmatches(time_line, gregexpr("[[:digit:]]+(\\.[[:digit:]]+)?", time_line)))
    if (inherits(min_sec, "try_error")) return(NA) 
    else min_sec <- unlist(min_sec)
    mins <- ifelse(length(min_sec) > 1, min_sec[1], 0)
    secs <- ifelse(length(min_sec) > 1, min_sec[2], min_sec)
    sec <- try(as.numeric(mins)*60 + as.numeric(secs))
    return(sec)
}


get_bl_ests <- function(bl_est_fname) {
    bl_ests <- read.csv(bl_est_fname, skip=1, head=TRUE, sep="\t")
    return(bl_ests)
}



# gather_results <- function(out_dir, n_taxa, nsites, perc_missing, bl_avg, al, nsim=10) {
# 
# 
#     # set up output dataframe:
#     nrows <- nsim * length(perc_missing) * length(bl_avg) * length(nsites)
#     res_names <- c("tax",  "site", "mrate", "lkhd", "bl", "i", "splits",          
#     "splits_in_common", "bl_1"         ,   "bl_2"         , "bl_3"         ,  "bl_4"        ,  "bl_5"           , 
#     "bl_1_est"        , "bl_2_est"     ,   "bl_3_est"     , "bl_4_est"     ,  "bl_5_est"    ,  "bl_1_est_lb"    , 
#     "bl_2_est_lb"     , "bl_3_est_lb"  ,   "bl_4_est_lb"  , "bl_5_est_lb"  ,  "bl_1_est_ub" ,  "bl_2_est_ub"    , 
#     "bl_3_est_ub"     , "bl_4_est_ub"  ,   "bl_5_est_ub"  ,
#     "rate_1"          , "rate_2"       ,   "rate_3"       , "rate_4"       ,  "rate_5"       ,  "rate_6"         , 
#     "rate_est_1"      , "rate_est_2"   ,   "rate_est_3"   , "rate_est_4"   ,  "rate_est_5"   , 
#     "rate_est_6"      , "rate_est_1_lb",   "rate_est_2_lb", "rate_est_3_lb",  "rate_est_4_lb",  "rate_est_5_lb"  , 
#     "rate_est_6_lb"   , "rate_est_1_ub",   "rate_est_2_ub", "rate_est_3_ub",  "rate_est_4_ub",  "rate_est_5_ub"  , 
#     "rate_est_6_ub"   , "freq_1"       ,   "freq_2"       , "freq_3"       ,  "freq_4"       ,  "freq_est_1"     , 
#     "freq_est_2"      , "freq_est_3"   ,   "freq_est_4"   , "freq_est_1_lb",  "freq_est_2_lb",  "freq_est_3_lb"  , 
#     "freq_est_4_lb"   , "freq_est_1_ub",   "freq_est_2_ub", "freq_est_3_ub",  "freq_est_4_ub",  
#     # "alpha"           ,  "alpha_est"       , "alpha_est_lb" ,   "alpha_est_ub" , 
#     "tl"              ,  "tl_est"      ,  "tl_est_lb"     , "tl_est_ub"    , 
#     "cmp_seconds")  
#     
#     res <- data.frame(matrix(nrow=nrows, ncol=length(res_names)))
#     names(res) <- res_names
#     
#     k <- 1 
#     for (tax in n_taxa) {
#         for (bl in bl_avg) {
#     
#             tr_fname <- sprintf("%s/%s.%s.tree.nex", out_dir, tax, bl)
#             tr_true <- ape::read.nexus(tr_fname)   
#     
#             for (site in nsites) {
#                 for (mrate in perc_missing) {
#                     for (lkhd in c("full", "pw", "pwh")) {
#                         for (i in 1:nsim) {
# 
#                             sim_fname <- sprintf("%s.%s.%s.%s.%s.%s",site,tax,lkhd,mrate,bl,i)
# 
#                             sump_fname <- sprintf("%s/parms.summ.%s.%s.%s.%s.%s.%s.pstat",out_dir,site,tax,lkhd,mrate,bl,i)
#                             contree_fname <- sprintf("%s/tree.summ.%s.%s.%s.%s.%s.%s.con.tre",out_dir,site,tax,lkhd,mrate,bl,i)
#                             res_fname <- sprintf("%s/res.%s.%s.%s.%s.%s.%s",out_dir,site,tax,lkhd,mrate,bl,i)
#                             bl_est_fname <- sprintf("%s/tree.summ.%s.%s.%s.%s.%s.%s.vstat",out_dir,site,tax,lkhd,mrate,bl,i)
#     
#                             freq_fname <- sprintf("%s/freq.%s.%s.%s.%s.%s.%s",out_dir,site,tax,lkhd,mrate,bl,i)
#                             rate_fname <- sprintf("%s/rate.%s.%s.%s.%s.%s.%s",out_dir,site,tax,lkhd,mrate,bl,i)
#     
#                             tr_freqs <- as.numeric(readLines(freq_fname))
#                             tr_rates <- as.numeric(readLines(rate_fname))
#     
#                             # compute time 
#                             cmp_seconds <- get_compute_time(res_fname)
#     
#                             # parameter estimates
#                             parm_ests <- get_parm_ests(sump_fname)
#     
#                             # consensus tree estimate
#                             tr_est <- ape::read.nexus(contree_fname)[[2]]
#                             bl_est <- get_bl_ests(bl_est_fname)
#                             tr_comp <- compare_contree(tr_true, tr_est, bl_est)
#                             if (ncol(tr_comp) !=  22) cat(sprintf("k = %s, tree_comp wrong num cols... \n", k))
#     
#                             tl <- sum(tr_true$edge.length)
#                             parm_comp <- make_parm_est(tr_rates, tr_freqs, parm_ests, tl)
#                             if (ncol(parm_comp) !=  44) cat(sprintf("k = %s, parm_comp wrong num cols... \n", k))
#     
#                             sim_res <- cbind(tr_comp, parm_comp, cmp_seconds)
#                             res_row <- cbind(tax, site, mrate, lkhd, bl, i, sim_res)
#     
#                             if (ncol(res_row) > ncol(res)) {cat("breaking... \n ");}
#                             res[k,] <- res_row
#                             k <- k+1
#                         }
#                     }
#                 }
#             }
#         }
#     }
# 
#     return(res)
# }
# 
# read_lkhds <- function(out_dir,site,tax,lkhd,mrate,bl,sim, numchains=2) {
# 
#     sim_fname <- sprintf("%s.%s.%s.%s.%s.%s",site,tax,lkhd,mrate,bl,sim)
#     lkhd_list <- lapply(1:numchains, function(i) {
#         fname <- sprintf("%s/%s.nex.run%s.p",out_dir,sim_fname, i)
#         res <- read.csv(fname, skip=1, row.names=NULL, sep="\t") 
#         data.frame(gen = res$Gen, lnLike=res$lnLike, chain=i)
#                    
#     })
#     out <- dplyr::bind_rows(lkhd_list)
#     out$site <- site
#     out$tax <- tax
#     out$mrate <- mrate
#     out$bl <- bl
#     out$sim <- sim
#     out$lkhd <- lkhd
#     return(out)
# }

# gather_lkhds <- function(out_dir, n_taxa, nsites, perc_missing, bl_avg, al, nsim=10, lkhds=c("pwh", "full")) {
# 
#     nouts <- nsim * length(perc_missing) * length(bl_avg) * length(nsites) * length(lkhds) * length(n_taxa)
#     lkhd_list_out <- vector("list", nouts)
# 
#     j <- 1
#     for (tax in n_taxa) {
#         for (bl in bl_avg) {
#             for (site in nsites) {
#                 for (mrate in perc_missing) {
#                     for (lkhd in lkhds) {
#                          for (i in 1:nsim) {
#                             res <- read_lkhds(out_dir,site,tax,lkhd,mrate,bl,i,2)
#                             lkhd_list_out[[j]] <- res
#                             j <- j + 1
#                          }
#                     }
#                 }
#             }
#         }
#     }
#     dplyr::bind_rows(lkhd_list_out)
# }



#h n_taxa <- 7
# n_site <- 10000
# 
# tree <- ape::rtree(n_taxa)
# tree$edge.length <- rep(0.1, length(tree$edge.length))
# sum(tree$edge.length)
# 
# plot(tree)
# edgelabels(round(tree$edge.length,2))
# 
# al <- 0.3
# be <- 0.7
# 
# Q <- matrix(c(-2*al, 2*al, 0, 
#               be, -(al+be), al, 
#               0, 2*be, -2*be), byrow=TRUE, nrow=3)
# # scale Q:
# qscaler <- -sum(diag(Q)*freq)
# Q <- Q/qscaler
# 
# a_scl <- Q[2,1]
# b_scl <- Q[2,3]
# freq <- c(b_scl^2, 2*a_scl*b_scl, a_scl^2)/((a_scl+b_scl)^2)
# 
# # simSeq func multiplies columns by frequencies -- pre-multiply to cancel out and use the correct Q matrix  
# Qsim <- Q %*% diag(1/freq)
# 
# 
# tree$edge.length
# 
# aln <- phangorn::simSeq(tree, n_site, Q=Qsim, bf=freq, type="USER", levels=c("E", "M", "D"))
# aln_char <- as.character(aln)
# 
# table(aln_char)/(nrow(aln_char)*ncol(aln_char))
# freq
# 
# seq_strings <- apply(aln_char, 1, function(x) paste0(x, collapse=""))
# seq_lines <- paste(names(seq_strings), seq_strings, sep="    ")
# 
# fname <- "phangorn"
# writeLines(c("#NEXUS\n",
# 
# # write the data block
#            "begin data;", 
#            sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
#            "format datatype=dimethyl interleave=no missing=?;",
#            "matrix",
# 
#            seq_lines,
#            "\t;\n end;\n",
# 
# # write the mrbayes block
#             "begin mrbayes;",
#             "mcmc ngen=30000 samplefreq=600 burnin=6000;",
#             "sump;",
#             "sumt;",
#             "end;"),
#          sprintf("%s.nex", fname))
# 
# system(sprintf("./mb %s.nex > %s.out",fname,fname))
# 


# #######
# ## now sim using iqtree
# #######
# 
# n_taxa <- 6
# n_site <- 30000
# 
# a <- 0.3
# b <- 0.7
# 
# # a <- a/b; b <- b/b
# fr <- c(b, a)/(a+b)
# 
# qmat <- matrix(c(-a, a, b, -b), nrow=2, byrow=TRUE)
# scaler <- -sum(diag(qmat)*fr)
# 
# qmat <- qmat/scaler
# qmat <- qmat %*% diag(1/fr)
# 
# tree_str <- write.tree(tree, "sim.tree")
# 
# # start with the same ancestral 
# iqtree_sim_str <- sprintf("./iqtree2 --alisim %s --length %s -m GTR2\\{%s\\}+F\\{%s,%s\\} -t sim.tree --num-alignments 2 > alnsim.out", "sim_aln_iq", n_site, a/b,fr[1],fr[2])
# system(iqtree_sim_str)
# 
# 
# # combine files:
# allele1 <- matrix(as.numeric(ape::read.dna("sim_aln_iq_1.phy", as.character=TRUE)),
#                     nrow=6, ncol=n_site)
# allele2 <- matrix(as.numeric(ape::read.dna("sim_aln_iq_2.phy", as.character=TRUE)),
#                     nrow=6,ncol=n_site)
# sim_aln <- allele1 + allele2
# seq_strings <-  apply(sim_aln, 1, function(x) paste0(c("E", "M", "D")[x+1], collapse=""))
# 
# seq_names <- paste("t", 1:6, sep="")
# seq_lines <- paste(seq_names, seq_strings, sep="    ")
# 
# fname <- "iq"
# # write mrbayes file
# writeLines(c("#NEXUS\n",
# 
# # write the data block
#            "begin data;", 
#            sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
#            "format datatype=dimethyl interleave=no missing=?;",
#            "matrix",
# 
#            seq_lines,
#            "\t;\n end;\n",
# 
# # write the mrbayes block
#             "begin mrbayes;",
#             "mcmc ngen=30000 samplefreq=600 burnin=6000;",
#             "sump;",
#             "sumt;",
#             "end;"),
#          sprintf("%s.nex", fname))
# 
# system(sprintf("./mb %s.nex > %s.out",fname,fname))




########
### try with seq-gen  <-  set  G/T freqs to 0, then rate a -> c is the methylization rate (a) 
########
#
#seqgen_path <- "~/Lib/SeqGen/source/seq-gen"
#
#
#tree_string <- TreeTools::NewickTree(tree)
#tree_file <- "sim.tree"
#
## simulate some data using seqgen  
#sg_call <- paste(seqgen_path, " -mGTR ", 
#                " -f ", paste0(freq,collapse=" "),  
#                " -r ", paste0(rate,collapse=" "), 
#                " -l ", n_site, " -a ", alpha,
#                " -z ", floor(runif(1)*3928109+2817615), 
#                ifelse(dg>0, sprintf(" -g %s ",dg), ""),
#                " < sim.tree ",
#                " > alignment.txt",sep="")
#system(sg_call)
#
#alignment_char <- phybase::read.dna.seq("alignment.txt","phylip")$seq
#alignment <- matrix(as.numeric(factor(alignment_char,levels=c("A","C","G","T"))),
#                    nrow=n_taxa,ncol=n_site)
#
#
## 
# setup_mrb_part <- function(part_aln_str,
#                            part_sites,
#                       n_taxa, n_site, n_parts,
#                       out_dir="sim",
#                       fname="",
#                       nopart=FALSE
#                       ) {
# 
#     site_ends <- cumsum(part_sites)
#     site_starts <- c(1, site_ends[-n_parts]+1)
#     ## now write the alignment:
#     writeLines(c(
#                  "#NEXUS\n",
# 
#         # write the data block
#         "begin data;", 
#         sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
#         "format datatype=dimethyl interleave=yes gap=- missing=?;",
#         "matrix",
#         unlist(part_aln_str),
#         "\t;\nend;\n",
# 
#         # write the partition block
#         "begin mrbayes;",
#         unlist(lapply(1:n_parts, function(i) sprintf("charset part%s = %s-%s;", i, site_starts[i], site_ends[i]))),
#         sprintf("partition allparts = %s : %s;", n_parts, paste0("part",1:n_parts, collapse=", ")),
# 
#         # write the model block
#         "begin mrbayes;",
#         "set partition=allparts;",
#         "lset rates=gamma;",
#         ifelse(nopart, "", "unlink shape=(all) ;"),
#         ifelse(nopart, "", "unlink dimethylrate=(all);"),
#         "mcmc ngen=50000 samplefreq=400 burnin=10000;",
#         sprintf("sump outputname=%s/sump.%s;",out_dir,fname),
#         sprintf("sumt outputname=%s/sumt.%s conformat=simple;",out_dir,fname),
# 
#         # sprintf("comparetree outputname=treecomp"),
#         #         "filename1=tree.summ.con.tre filename2=sim.tree.nex;",
#         "end;"
#               ),
#         sprintf("%s/%s.nex",out_dir,fname))
# }
# 
# 
# sim_seqgen <- function(tree, a, b, n_site,
#                      al=1,  
#                      dg=0,
#                      tree_file="sim.tree", aln_file="aln.txt" 
#                      ) {
#     seq_names <- tree$tip.label
#     seqgen_path <- "~/Lib/SeqGen/source/seq-gen"
# 
#     qmat <- matrix(c(-a,a,b,-b),nrow=2,byrow=TRUE)
#     freq <- c(b,a) / (a+b)
#     qmat <- qmat/-sum(freq*diag(qmat))
#     
#     # pre-mult by pi:
#     qmat <- qmat%*%diag(1/freq)
#  
#     rate <- c(qmat[1,2],qmat[2,1],1,1,1,1)
#     freq <- c(b,a,0,0)
# 
#     tree_string <- TreeTools::NewickTree(tree)   
#     # store tree w/ branch lengths to file
#     writeLines(tree_string, con=tree_file)
# 
#     al_line <- ifelse(dg == 0, sprintf(" -a %s",al), sprintf(" -a %s -g %s", al, dg))
#     # simulate some data using seqgen  
#     sg_call <- ifelse(al > 0, 
#                  paste(seqgen_path, " -mGTR ", 
#                     " -f", paste0(freq,collapse=" "),  
#                     " -r", paste0(rate,collapse=" "), 
#                     " -l", n_site, 
#                     al_line,
#                     " -on", 
#                     " -n2",
#                     " -z1000", 
#                     " < ", tree_file,
#                     sep=""),
#                  paste(seqgen_path, " -mGTR ", 
#                     " -f", paste0(freq,collapse=" "),  
#                     " -r", paste0(rate,collapse=" "), 
#                     " -l", n_site, 
#                     " -on",
#                     " -n2",
#                     " -z1000", 
#                     " < ", tree_file,
#                     sep=""))
# 
#     system(sprintf("%s > aln.txt", sg_call))
# 
#     # delete comment lines:
#     l1 <- readLines("aln.txt")
#     l1 <- l1[c(1,22:(22+3+length(tree$tip.label)+2))]
#     writeLines(l1, "aln1.txt")
#     l2 <- readLines("aln.txt")
#     l2 <- l2[c(1,((22+3+length(tree$tip.label)+2+1):(22+3+length(tree$tip.label)+2+1+3+length(tree$tip.label)+2+1)))]
#     writeLines(l2, "aln2.txt")
# 
#     aln1_char <- ape::read.nexus.data("aln1.txt") 
#     aln2_char <- ape::read.nexus.data("aln2.txt") 
# 
#     aln1 <- apply(matrix(unlist(aln1_char), ncol=n_site), 2, function(x) 
#                        match(x, c("a","c")))-1
#     aln2 <- apply(matrix(unlist(aln2_char), ncol=n_site), 2, function(x) 
#                        match(x, c("a","c")))-1
# 
#     aln <- aln1+aln2
#     seq_strings <-  apply(aln, 1, function(x) paste0(c("E", "M", "D")[x+1], collapse=""))
#     seq_lines <- paste(seq_names, seq_strings, sep="    ")
# 
# 
# #     if (missing_rate > 0) {
# #         tot_chars <- nrow(alignment) * ncol(alignment)
# #         del_sites <- sample(1:tot_chars, round(tot_chars * missing_rate,0))
# #         alignment[del_sites] <- "-"
# #         alignment_char[del_sites] <- "-"
# #     }
# # 
# #     alignment_str <- apply(alignment_char, 1, 
# #                        function(x) paste0(x,collapse=""))
#  
#     return(seq_lines)
# 
# }


# write_mrb <- function(aln, n_taxa, n_site, fname) {
# 
#     # write mrbayes file
#     writeLines(c("#NEXUS\n",
#     
#     # write the data block
#                "begin data;", 
#                sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
#                sprintf("format datatype=dimethyl interleave=%s missing=?;",interleave),
#                "matrix",
#     
#                seq_lines,
#                "\t;\n end;\n",
#     
#     # write the mrbayes block
#                 "begin mrbayes;",
#                 "mcmc ngen= samplefreq= burnin=;",
#                 sprintf("sump outputname=%s;", fname),
#                 sprintf("sumt outputname=%s conformat=simple;", fname),
#                 "end;"),
#              sprintf("%s.nex", fname))
#    
# }
# run_simulation <- function(out_dir, n_taxa, nsites, perc_missing, bl_avg, al, nsim=10, 
#                            recook=TRUE, mrb_path="./mb") {
# 
#     for (tax in n_taxa) {
#         for (bl in bl_avg) {
#     
#             tr <- gen_tree(tax, bl)
#             tr_str <- TreeTools::NewickTree(tr)
#             tr_fname <- sprintf("%s/%s.%s.tree.nex", out_dir, tax, bl)
#             phybase::write.tree.string(tr_str, format="Nexus", file=tr_fname) 
#     
#             for (site in nsites) {
#                 for (mrate in perc_missing) {
#                     for (i in 1:nsim) {
#                         rates <- c(1,2,3,4,7,8) # rlnorm(6, 0.1, 0.2) 
#                         rates <- rates/sum(rates)
#     
#                         freqs <- c(.3,.2,.1,.4) # rlnorm(4, 0.1, 0.2)
#                         freqs <- freqs/sum(freqs)
#     
#                         aln <- sim_aln(tr,site,al=al,rate=rates,freq=freqs,missing_rate=mrate)
#     
#                         for (lkhd in c("full", "pw", "pwh")) {
#                     
#                             sim_fname <- sprintf("%s.%s.%s.%s.%s.%s",site,tax,lkhd,mrate,bl,i)
#                             sump_fname <- sprintf("parms.summ.%s",sim_fname)
#                             sumt_fname <- sprintf("tree.summ.%s",sim_fname)
#                             out_fname <- sprintf("%s/res.%s",out_dir,sim_fname)
#                              
#                             freq_fname <- sprintf("%s/freq.%s.%s.%s.%s.%s.%s",out_dir,site,tax,lkhd,mrate,bl,i)
#                             rate_fname <- sprintf("%s/rate.%s.%s.%s.%s.%s.%s",out_dir,site,tax,lkhd,mrate,bl,i)
#                             
#                             writeLines(paste0(freqs), freq_fname)
#                             writeLines(paste0(rates), rate_fname)
#     
#                             setup_mrbayes_file(aln$aln_str,tr,nrow(aln$aln),ncol(aln$aln),
#                                    lkhd,gamma_rates=(al>0.0),out_dir,sim_fname,sump_fname,sumt_fname)
#     
#                             if (!recook) if (file.exists(out_fname)) next 
#     
#                             run_mrbayes_file(mrb_path, sprintf("%s/%s.nex",out_dir,sim_fname), out_fname)
#                         }
#                     }
#                 }
#             }
#         }
#     }
# }



