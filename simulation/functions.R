
sim_alignment <- function(tree, a, b, n_site, al, dg=0, readerr=0) {
    #  tree: phylo object 
    # a, b : methylation add and gain rates
    # n_site: length to generate
    # al: alpha for site rate heterogeneity
    # dg: number of discrete rate catefories
    # readerr: read error rate

    seq_names <- tree$tip.label

    # init and normalize q matrix
    qmat <- matrix(c(-a,a,b,-b),nrow=2,byrow=TRUE)
    freq <- c(b,a) / (a+b)
    qmat <- qmat/-sum(freq*diag(qmat))
    
    # pre-mult by pi (internally simseq will premult by diag(freq) to get back to the correct q matrix)
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

setup_mrb <- function(aln_str,
                      n_taxa, n_site,
                      model_lines,
                      mcmc_lines,
                      out_dir="sim",
                      fname="",
                      datatype="dimethyl",
                      interleave="no"
                      ) {

    # aln_str: alignment as a character vector
    # n_taxa, n_site: number of taxa, number of sites
    # model_lines, mcmc_lines: lines to write to nexus file specifying model and mcmc parameters
    # out_dir: directory to write file to
    # fname: name of nexus file

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
                 sprintf("sump outputname=%s/%s burninfrac=0.5;",out_dir,fname),
                 sprintf("sumt outputname=%s/%s burninfrac=0.5 conformat=figtree;",out_dir,fname),

                 # sprintf("comparetree outputname=treecomp"),
                 #         "filename1=tree.summ.con.tre filename2=sim.tree.nex;",
                 "end;"

                 ),
                sprintf("%s/%s.nex",out_dir,fname))
}

# wrapper for system call to run nexus file
run_mrb <- function(mrb_path, nex_file, out_file) {
        system(sprintf("%s %s &> %s", mrb_path, nex_file,  out_file), intern=FALSE)
}

