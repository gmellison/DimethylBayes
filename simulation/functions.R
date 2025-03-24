
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

    seq_strings <-  apply(sim_aln, 1, function(x) paste0(c("0", "1", "2")[x+1], collapse=""))
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


