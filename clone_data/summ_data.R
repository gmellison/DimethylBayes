library(dplyr)
library(phangorn)
library(ape)
library(phybase)
library(TreeTools)
library(ggplot2)
library(stringr)

#source("functions.R")

nout <- 12 

age_root_est <- matrix(NA, nrow=nout, ncol=3)
clock_rate_est <- matrix(NA, nrow=nout, ncol=3)

nparts <- rep(NA, nout)
clock_mod <- rep(NA, nout) 
clone_mod <- rep(NA, nout) 
clock_rate_out <- rep(NA, nout)
like_out <- rep(NA, nout)
tree_out <- rep(NA, nout)

k <- 1
for (k in 1:nout) {
    data_dir <- sprintf("server_res/%s",k)

    log_fname <- list.files(data_dir, "*.log")
    fname     <- str_split_fixed(log_fname, "\\.", 2)[,1]
    run_split <- str_split_fixed(fname, "_", 3) 

    which_clone  <- run_split[,1]
    which_part   <- run_split[,2]
    which_clock  <- run_split[,3]

    nex_file <- sprintf("%s.nex",fname)

    #tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
    #tree_tr <- read.tree(tree_f)

    #tree_lines <- readLines(sprintf("%s/sumt.%s.trprobs",out_dir, fname))
    #tree_lines <- tree_lines[stringr::str_detect(tree_lines, "&W")]
    #tree_probs <- as.numeric(stringr::str_split_fixed(tree_lines, "[\\]W]", 4)[,3])

    #trees_est <- ape::read.nexus(sprintf("%s/sumt.%s.trprobs",out_dir,fname), force.multi=TRUE)
    #is_true_tree <- unlist(lapply(trees_est, function(sampled_tree) {
    #                                      all.equal.phylo(tree_tr, sampled_tree, 
    #                                                      use.edge.length=FALSE)}))
    #true_tree_prob[k] <- ifelse(any(is_true_tree), tree_probs[is_true_tree], 0)

    ## get age of root
    if (which_clock != "") {
        #vstat <- read.csv(sprintf("%s/%s.vstat", data_dir, nex_file), sep="\t", skip=1)
        #age_root <- as.matrix(vstat[vstat$Parameter == "age[0]", c(2,4,5)])
        #age_root_est[k,] <- age_root
    
        ## get clock rate
        pstat <- read.csv(sprintf("%s/%s.pstat", data_dir, nex_file), sep="\t", skip=1)
        clock_rate <- as.matrix(pstat[pstat$Parameter == "clockrate", c(2,4,5)])
        clock_rate_est[k,] <- clock_rate
    } else {
        clock_rate_est <- c("","","")
    }
    

    # get the marginal likelihoods of clock and nonclock
    #     models
    like <-    read.csv(sprintf("%s/%s.lstat", data_dir, nex_file), sep="\t", skip=1)
    like_out[k] <- like$harm[like$run=="all"]
    clock_mod[k] <- which_clock
    nparts[k] <- which_part
    clone_mod[k] <- which_clone
}

res <- data.frame(clone=clone_mod, part=nparts, clock=clock_mod, like=like_out)

