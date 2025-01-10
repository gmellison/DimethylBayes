library(phangorn)
source("functions.R")
library(dplyr)
library(ggplot2)

#
#      CLOCK       ####
#
#

nsim <- 500 

if (!file.exists("res_cr.csv")) {

tree_dir <- "trees"
out_dir <- "clock"

# set sim params
a <- b <- 0.5
sites <- c(100, 1000, 10000)
taxa <- c(4,20)
bl <- 0.1
ages <- c(20, 100)
crs <- c(0.0001, 0.0007)

# setup output structures
# site <- 1000
# j <- 1
# tax <- 4
# bl <- 0.1
# cr <- 0.1
# cr_ratio <- 0.1

nout <- nsim * length(sites) * length(crs) * length(taxa) * length(ages)
true_tree_prob <- rep(NA, nout)

age_root_est <- matrix(NA, nrow=nout, ncol=3)
clock_rate_est <- matrix(NA, nrow=nout, ncol=3)

clock_rate_out <- rep(NA, nout)
tree_age_out <- rep(NA, nout)
rep_out <- rep(NA, nout)
sites_out <- rep(NA, nout)
age_out <- rep(NA, nout)
tax_out <- rep(NA, nout)
like_out <- rep(NA, nout)
like_nc_out <- rep(NA, nout)
time_out <- rep(NA, nout)

k <- 1
for (site in sites) {
    for (tax in taxa) {
        for (age in ages) {
            for (cr in crs) {
                for (j in 1:nsim) {
                    fname <- sprintf("%s.%s.%.1e.%s.%s.m",tax,age,cr,site,j)

                    tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
                    tree_tr <- read.tree(tree_f)

                    tree_lines <- readLines(sprintf("%s/sumt.%s.trprobs",out_dir, fname))
                    tree_lines <- tree_lines[stringr::str_detect(tree_lines, "&W")]
                    tree_probs <- as.numeric(stringr::str_split_fixed(tree_lines, "[\\]W]", 4)[,3])

                    trees_est <- ape::read.nexus(sprintf("%s/sumt.%s.trprobs",out_dir,fname), force.multi=TRUE)
                    is_true_tree <- unlist(lapply(trees_est, function(sampled_tree) {
                                                          all.equal.phylo(tree_tr, sampled_tree, 
                                                                          use.edge.length=FALSE)}))
                    true_tree_prob[k] <- ifelse(any(is_true_tree), tree_probs[is_true_tree], 0)

                    ## get age of root
                    vstat <- read.csv(sprintf("%s/sumt.%s.vstat", out_dir,fname), sep="\t", skip=1)
                    age_root <- as.matrix(vstat[vstat$Parameter == "age[0]", c(2,4,5)])
                    age_root_est[k,] <- age_root
    
                    ## get clock rate
                    pstat <- read.csv(sprintf("%s/sump.%s.pstat", out_dir,fname), sep="\t", skip=1)
                    clock_rate <- as.matrix(pstat[pstat$Parameter == "clockrate", c(2,4,5)])
                    clock_rate_est[k,] <- clock_rate


                    # get the marginal likelihoods of clock and nonclock
                    #     models
                    like <-    read.csv(sprintf("%s/sump.%s.lstat", out_dir,fname), sep="\t", skip=1)
                    like_nc <- read.csv(sprintf("%s/sump.%s.nc.lstat", out_dir,fname), sep="\t", skip=1)
                    like_out[k] <- like$harm#[like$run == "all"]
                    like_nc_out[k] <- like_nc$harm#[like_nc$run == "all"]

                    time_out[k] <- get_compute_time(sprintf("%s/%s.log", out_dir, fname))

                    # fill in the true vals
                    clock_rate_out[k] <- cr 
                    sites_out[k] <- site
                    rep_out[k] <- j 
                    tax_out[k] <- tax
                    age_out[k] <- age 

                    k <- k + 1
                }
            }
        }
    }
}

# make data frame
res_cr <- data.frame(sites = sites_out,
                  tax = tax_out,
           rep = rep_out,
           cr = clock_rate_out,
           age = age_out,
           true_tree_prob  = true_tree_prob,
           age_est   =  age_root_est[,1],
           age_est_lb  =  age_root_est[,2],
           age_est_ub  =  age_root_est[,3],
           cr_est   = clock_rate_est[,1],
           cr_est_lb   = clock_rate_est[,2],
           cr_est_ub   = clock_rate_est[,3],
           like        = like_out,
           like_nc     = like_nc_out,
           time = time_out
)

write.csv(res_cr, "./res_cr.csv")

} else {
    res_cr <- read.csv("res_cr.csv")
}

tree_sum_cr <- res_cr %>% group_by(tax, cr, age, sites) %>% summarize(tr_tree=mean(true_tree_prob)) 
cr_sum_cr <- res_cr %>% group_by(tax, cr, age, sites) %>% 
        summarize(cr_in_ci = mean(cr_est_lb < cr & cr_est_ub > cr),
                  cr_rmse = sqrt(mean((cr_est - cr)^2)),
                  avg_int_len = mean(cr_est_ub - cr_est_lb))

bf_sum_cr <- res_cr %>% group_by(tax, cr, age, sites) %>% 
          summarize(scm = sum(like - like_nc >= 2 & like - like_nc < 6),
                    scs  = sum(like - like_nc > 6),
                    sncm = sum(like_nc - like > 2 & like_nc - like < 6),
                    sncs = sum(like_nc - like > 6),
                    cm = mean(like - like_nc >= 2 & like - like_nc < 6),
                    cs  = mean(like - like_nc > 6),
                    ncm = mean(like_nc - like > 2 & like_nc - like < 6),
                    ncs = mean(like_nc - like > 6))


time_cr <- res_cr %>% group_by(tax, cr, age, sites) %>% 
        summarise(compute_time = mean(time, na.rm=TRUE))

#
#
#      AGE            ####
#
#

if (!file.exists("res_age.csv")) {
out_dir <- "age"
tree_dir <- "trees"

crs <- c(0.0001, 0.0007)

a <- b <- 0.5
al <- 0.5
sites <- c(100, 1000, 10000)
taxa <- c(4,20)
bl <- 0.1
nsim <- 500
ages <- c(20,100)

nout <- nsim * length(sites) * length(crs) * length(ages) * length(taxa)
true_tree_prob <- rep(NA, nout)

age_root_est <- matrix(NA, nrow=nout, ncol=3)
clock_rate_est <- matrix(NA, nrow=nout, ncol=3)

clock_rate_out <- rep(NA, nout)
tree_age_out <- rep(NA, nout)
rep_out <- rep(NA, nout)
sites_out <- rep(NA, nout)
like_out <- rep(NA, nout)
like_nc_out <- rep(NA, nout)
tax_out <- rep(NA, nout)
age_out <-  rep(NA, nout)
time_out <- rep(NA, nout)


k <- 1
for (site in sites) {
    for (tax in taxa) {
        for (j in 1:nsim) {
            for (cr in crs) {
                for (age in ages) {
                    fname <- sprintf("%s.%s.%.1e.%s.%s.m",tax,age,cr,site,j)

                    tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
                    tree_tr <- read.tree(tree_f)

                    tree_lines <- readLines(sprintf("%s/sumt.%s.trprobs",out_dir,fname))
                    tree_lines <- tree_lines[stringr::str_detect(tree_lines, "&W")]
                    tree_probs <- as.numeric(stringr::str_split_fixed(tree_lines, "[\\]W]", 4)[,3])

                    trees_est <- ape::read.nexus(sprintf("%s/sumt.%s.trprobs",out_dir,fname), force.multi=TRUE)
                    is_true_tree <- unlist(lapply(trees_est, function(sampled_tree) {
                                                          all.equal.phylo(tree_tr, sampled_tree, 
                                                                          use.edge.length=FALSE)}))
                    true_tree_prob[k] <- ifelse(any(is_true_tree), tree_probs[is_true_tree], 0)

                    ## get age of root
                    vstat <- read.csv(sprintf("%s/sumt.%s.vstat", out_dir, fname), sep="\t", skip=1)
                    age_root <- as.matrix(vstat[vstat$Parameter == "age[0]", c(2,4,5)])
                    age_root_est[k,] <- age_root
    
                    ## get clock rate
                    pstat <- read.csv(sprintf("%s/sump.%s.pstat", out_dir, fname), sep="\t", skip=1)
                    clock_rate <- as.matrix(pstat[pstat$Parameter == "clockrate", c(2,4,5)])
                    clock_rate_est[k,] <- clock_rate
    
                    # get the marginal likelihoods of clock and nonclock
                    #     models
                    like <-    read.csv(sprintf("%s/sump.%s.lstat", out_dir, fname), sep="\t", skip=1)
                    like_nc <- read.csv(sprintf("%s/sump.%s.nc.lstat", out_dir, fname), sep="\t", skip=1)
                    like_out[k] <- like$harm
                    like_nc_out[k] <- like_nc$harm
                    time_out[k] <- get_compute_time(sprintf("%s/%s.log", out_dir, fname))

   
                    # fill in the true vals
                    #mod_out[k] <- mod
                    clock_rate_out[k] <- cr 
                    sites_out[k] <- site
                    rep_out[k] <- j 
                    tax_out[k] <- tax
                    age_out[k] <- age
                    k <- k + 1
                }
            }
        }
    }
}

# make data frame
res_age <- data.frame(sites = sites_out,
                  tax = tax_out,
           rep = rep_out,
           cr = clock_rate_out,
           age = age_out,
           true_tree_prob  = true_tree_prob,
           age_est   =  age_root_est[,1],
           age_est_lb  =  age_root_est[,2],
           age_est_ub  =  age_root_est[,3],
           cr_est   = clock_rate_est[,1],
           cr_est_lb   = clock_rate_est[,2],
           cr_est_ub   = clock_rate_est[,3],
           like = like_out,
           like_nc = like_nc_out,
           time        = time_out
)

write.csv(res_age, "./res_age.csv")

} else  {
        res_age <- read.csv("res_age.csv")
}
res_age <- res_age %>% filter(sqrt((age_est - age)^2) < 100)

tree_sum_age <- res_age %>% group_by(tax, cr, age, sites) %>% summarize(tr_tree=mean(true_tree_prob))
ta_sum_age <- res_age %>% group_by(tax, cr, age, sites) %>% 
          summarize(ta_in_ci = mean(age_est_lb < age & age_est_ub > age),
                    ta_rmse = mean(sqrt((age_est - age)^2)),
                    avg_int_len = mean(age_est_ub - age_est_lb))

bf_sum_age <- res_age %>% group_by(tax, cr, age, sites) %>% 
          summarize(cm = mean(like - like_nc > 2),
                    cs  = mean(like - like_nc > 6),
                    ncm = mean(like_nc - like > 2),
                    ncs = mean(like_nc - like > 6),
                    n=nsim)

time_age <- res_age %>% group_by(tax, cr, age, sites) %>% 
        summarise(compute_time = mean(time, na.rm=TRUE))


#
#
#      SNP COMP          ####
#
#

# setup output structures
# site <- 1000
# j <- 1
# tax <- 4
# bl <- 0.1
# cr <- 0.1
# cr_ratio <- 0.1
# methyl & dna clock rates

if (!file.exists("res_snp.csv")) {

out_dir <- "snp_comp"
tree_dir <- "trees"

dm_cr <- c(0.0007)
snp_cr <- c(4*10^(-9)) * c(100, 1000, 10000)

clocks <- list(list(mod="m", cr=dm_cr[1]),
               #list(mod="m", cr=dm_cr[3]),
               list(mod="dna", cr = snp_cr[1]),
               list(mod="dna", cr=snp_cr[2])
            )

a <- b <- 0.5
sites <- c(1000, 10000, 100000)
j <- 1
taxa <- c(4, 20)
bl <- 0.1
nsim <- 500 
ages <- c(20)

numruns <- nsim * length(sites) * length(taxa) * length(clocks) * length(ages)

nout <- nsim * length(sites) * length(clocks) * length(taxa) * length(ages)
true_tree_prob <- rep(NA, nout)

age_root_est <- matrix(NA, nrow=nout, ncol=3)
clock_rate_est <- matrix(NA, nrow=nout, ncol=3)

clock_rate_out <- rep(NA, nout)
tree_age_out <- rep(NA, nout)
rep_out <- rep(NA, nout)
sites_out <- rep(NA, nout)
mod_out <- rep(NA, nout)
tax_out <- rep(NA, nout)
age_out <- rep(NA, nout)
like_out <- rep(NA, nout)
time_out <- rep(NA, nout)


k <- 1
for (site in sites) {
    for (age in ages) {
        for (tax in taxa) {
            for (clock in clocks) {
                for (j in 1:nsim) {

                    mod <- clock$mod
                    cr  <- clock$cr 

                    if  (mod == "m") {
                        fname <- sprintf("%s.%s.%.1e.%s.%s.m",tax,age,cr,site,j)                             
                        tree_f <- sprintf("%s/%s.%s.%s.tree", tree_dir, tax, cr, age)
                        tree_tr <- read.tree(tree_f)

                        tree_lines <- readLines(sprintf("%s/sumt.%s.trprobs",out_dir,fname))
                        tree_lines <- tree_lines[stringr::str_detect(tree_lines, "&W")]
                        tree_probs <- as.numeric(stringr::str_split_fixed(tree_lines, "[\\]W]", 4)[,3])

                        trees_est <- ape::read.nexus(sprintf("%s/sumt.%s.trprobs",out_dir,fname), force.multi=TRUE)
                        is_true_tree <- unlist(lapply(trees_est, function(sampled_tree) {
                                                              all.equal.phylo(tree_tr, sampled_tree, 
                                                                              use.edge.length=FALSE)}))
                        true_tree_prob[k] <- ifelse(any(is_true_tree), tree_probs[is_true_tree], 0)

                        ## get age of root
                        vstat <- read.csv(sprintf("%s/sumt.%s.vstat", out_dir, fname), sep="\t", skip=1)
                        age_root <- as.matrix(vstat[vstat$Parameter == "age[0]", c(2,4,5)])
                        age_root_est[k,] <- age_root
        
                        ## get clock rate
                        pstat <- read.csv(sprintf("%s/sump.%s.pstat", out_dir,fname), sep="\t", skip=1)
                        clock_rate <- as.matrix(pstat[pstat$Parameter == "clockrate", c(2,4,5)])
                        clock_rate_est[k,] <- clock_rate
   
                        # get the marginal likelihoods of clock and nonclock
                        #     models
                        like <-    read.csv(sprintf("%s/sump.%s.lstat", out_dir,fname), sep="\t", skip=1)
                        like_out[k] <- like$harm[like$run == "all"]
                        #like_nc_out[k] <- like_nc$harm[like_nc$run == "all"]

                        time_out[k] <- get_compute_time(sprintf("%s/%s.log", out_dir, fname))

                        # fill in the true vals
                        mod_out[k] <- mod
                        clock_rate_out[k] <- cr 
                        sites_out[k] <- site
                        age_out[k] <- age
                        rep_out[k] <- j 
                        tax_out[k] <- tax
                        k <- k + 1
                    }

                    if  (mod == "dna") {

                        tree_f <- sprintf("%s/%s.%s.%s.tree", tree_dir, tax, cr, age)
                        tree_tr <- read.tree(tree_f)
                        fname     <- sprintf("%s.%s.%.1e.%s.%s.dna",tax,age,cr,site,j)
                        tree_lines <- readLines(sprintf("%s/sumt.%s.trprobs",out_dir,fname))
                        tree_lines <- tree_lines[stringr::str_detect(tree_lines, "&W")]
                        tree_probs <- as.numeric(stringr::str_split_fixed(tree_lines, "[\\]W]", 4)[,3])

                        trees_est <- ape::read.nexus(sprintf("%s/sumt.%s.trprobs",out_dir,fname), force.multi=TRUE)
                        is_true_tree <- unlist(lapply(trees_est, function(sampled_tree) {
                                                              all.equal.phylo(tree_tr, sampled_tree, 
                                                                              use.edge.length=FALSE)}))
                        true_tree_prob[k] <- ifelse(any(is_true_tree), tree_probs[is_true_tree], 0)

                        ## get age of root
                        vstat <- read.csv(sprintf("%s/sumt.%s.vstat", out_dir,fname), sep="\t", skip=1)
                        age_root <- as.matrix(vstat[vstat$Parameter == "age[0]", c(2,4,5)])
                        age_root_est[k,] <- age_root
        
                        ## get clock rate
                        pstat <- read.csv(sprintf("%s/sump.%s.pstat", out_dir,fname), sep="\t", skip=1)
                        clock_rate <- as.matrix(pstat[pstat$Parameter == "clockrate", c(2,4,5)])
                        clock_rate_est[k,] <- clock_rate
        
                        time_out[k] <- get_compute_time(sprintf("%s/%s.log", out_dir, fname))

                        # fill in the true vals
                        mod_out[k] <- mod
                        age_out[k] <- age
                        clock_rate_out[k] <- cr
                        sites_out[k] <- site
                        rep_out[k] <- j 
                        tax_out[k] <- tax
                        k <- k + 1
                    }
                }
            }
        }
    }
}

# make data frame
res_snp <- data.frame(sites = sites_out,
                  tax = tax_out,
           rep = rep_out,
           mod = mod_out,
           cr = clock_rate_out,
           age = age_out,
           true_tree_prob  = true_tree_prob,
           age_est   =  age_root_est[,1],
           age_est_lb  =  age_root_est[,2],
           age_est_ub  =  age_root_est[,3],
           cr_est   = clock_rate_est[,1],
           cr_est_lb   = clock_rate_est[,2],
           cr_est_ub   = clock_rate_est[,3],
           like        = like_out,
           cmp_time        = time_out
)
write.csv(res_snp, "./res_snp.csv")
} else {
    res_snp <- read.csv("res_snp.csv")
}


tree_sum_snp <- res_snp %>% 
        group_by(mod, tax, cr, age, sites) %>% 
        summarize(tr_tree=mean(true_tree_prob)) 

cr_sum_snp <- res_snp %>% 
        group_by(mod, tax, cr, age, sites) %>% 
        summarize(cr_in_ci = mean(cr_est_lb < cr & cr_est_ub > cr),
                  cr_rmse  = sqrt(mean((cr_est - cr)^2)),
                  avg_int_len = mean(cr_est_ub - cr_est_lb),
                  cv = sqrt(mean((cr_est - cr)^2))/mean(cr_est),
                  cci = mean(cr_est_ub - cr_est_lb)/mean(cr_est) )

time_snp <- res_snp %>% 
        group_by(mod, tax, cr, age, sites) %>% 
        summarise(comp_time = mean(time,na.rm=TRUE))


#
#
#      PLOTS            ####
#
#
## CR
## plot 1: mse
plt_cr_mse <- ggplot(cr_sum_cr) +
       geom_point(aes(x=factor(sites),
                      y=cr_rmse,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=cr_rmse,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       labs(x = "Sites", y = "CR rMSE") + 
       scale_color_grey(name="Taxa")
ggsave("plots/sim_cr_rmse.pdf", plt_cr_mse, w=3,h=3,units="in")

# plot 2: ci coverage rate
plt_cr_cov <- ggplot(cr_sum_cr) +
       geom_point(aes(x=factor(sites),
                      y=cr_in_ci,
                      col=factor(tax),
                      shape=interaction(cr,age, sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=cr_in_ci,
                     col=factor(tax),
                     linetype=interaction(cr,age, sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "CR CI Coverage") + 
       scale_color_grey(name="Taxa") +
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       coord_cartesian(ylim=c(0.5,1)) 
plt_cr_cov
ggsave("plots/sim_cr_cov.pdf",plt_cr_cov,w=3,h=3,units="in")

# plot 3: avg ci width
plt_cr_ci_wid <- ggplot(cr_sum_cr) +
       geom_point(aes(x=factor(sites),
                      y=avg_int_len,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=avg_int_len,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "CI Width") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_cr_ci_width.pdf",plt_cr_ci_wid,w=3,h=3,units="in")

# plot 4: true tree posterior %
plt_tr_topo_cr <-  ggplot(tree_sum_cr) +
       geom_point(aes(x=factor(sites),
                      y=tr_tree,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=tr_tree,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "True Topo. Cov.") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa") 
ggsave("plots/sim_cr_tr_topo_prob.pdf",plt_tr_topo_cr,w=3,h=3,units="in")

# plot 5 
plt_cr_bf <- ggplot(bf_sum_cr) +
       geom_point(aes(x=factor(sites),
                      y=(cm + cs) - (ncm + ncs),
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                      y=(cm + cs) - (ncm + ncs),
                      col=factor(tax),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "Clock Preference") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_cr_bf.pdf",plt_cr_bf,w=3,h=3,units="in")


## AGE
# plot 1: mse
plt_age_rmse <- ggplot(filter(ta_sum_age))+
       geom_point(aes(x=factor(sites),
                      y=ta_rmse,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=ta_rmse,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "Tree Age rMSE") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa") 
ggsave("plots/sim_age_rmse.pdf",plt_age_rmse,w=3,h=3,units="in")

# plot 2: ci coverage rate
plt_age_cov <- ggplot(ta_sum_age) +
       geom_point(aes(x=factor(sites),
                      y=ta_in_ci,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=ta_in_ci,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       coord_cartesian(ylim=c(0,1)) + 
       labs(x = "Sites", y = "Tree Age CI Cov.") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_age_cov.pdf",plt_age_cov,w=3,h=3,units="in")

# plot 3: avg ci width
plt_age_ci_wid <- ggplot(ta_sum_age) +
       geom_point(aes(x=factor(sites),
                      y=avg_int_len,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                     y=avg_int_len,
                     col=factor(tax),
                     linetype=interaction(cr,age,sep="/"),
                     group=interaction(age,cr,tax))) + 
       theme_minimal() + 
       #coord_cartesian(ylim=c(0,1)) + 
       labs(x = "Sites", y = "CI Width") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_age_ci_width.pdf",plt_age_ci_wid,w=3,h=3,units="in")

# plot 4: true tree posterior %
plt_age_tr_topo <- ggplot(tree_sum_age) +
       geom_point(aes(x=factor(sites),
                      y=tr_tree,
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                      y=tr_tree,
                      col=factor(tax),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "True Topo. Cov.") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_age_tr_topo_prob.pdf",plt_age_tr_topo,w=3,h=3,units="in")

# plot 5 
plt_age_bf <- ggplot(bf_sum_age) +
       geom_point(aes(x=factor(sites),
                      y=cm+cs - (ncm + ncs),
                      col=factor(tax),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                      y=cm+cs - (ncm + ncs),
                      col=factor(tax),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) + 
       theme_minimal() + 
       labs(x = "Sites", y = "Clock Preference") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa")
ggsave("plots/sim_age_bf.pdf",plt_age_bf,w=3,h=3,units="in")


## SNP COMP
# plot 3: avg ci width
# plot 3: avg ci width

plt_snp_rmse <- ggplot(cr_sum_snp) +
       geom_point(aes(x=factor(sites),
                      y=cr_rmse,
                      col=factor(mod),
                      shape=interaction(cr,age,sep="/"))
                      ) +
       geom_line(aes(x=factor(sites),
                      y=cr_rmse,
                      col=factor(mod),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) +
       theme_minimal() + 
       labs(x = "Sites", y = "CR rMSE") + 
       scale_color_grey(name="Taxa") +
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       facet_wrap(~paste("Taxa:", tax)) +
       scale_y_log10()

ggsave("plots/sim_snp_rmse.pdf", plt_snp_rmse,w=3,h=3,units="in")

plt_snp_cr_ci_cov <- ggplot(cr_sum_snp) +
       geom_point(aes(x=factor(sites),
                      y=cr_in_ci,
                      col=factor(mod),
                      shape=interaction(cr,age,sep="/"))
                      ) +
       geom_line(aes(x=factor(sites),
                      y=cr_in_ci,
                      col=factor(mod),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) +
       theme_minimal() + 
       labs(x = "Sites", y = "CR CI Cov.") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Taxa") +
       facet_wrap(~paste("Taxa:", tax))
ggsave("plots/sim_snp_cr_ci_cov.pdf", plt_snp_cr_ci_cov,w=3,h=3,units="in")


plt_snp_cv <- ggplot(cr_sum_snp) +
       geom_point(aes(x=factor(sites),
                      y=cv,
                      col=factor(mod),
                      shape=interaction(cr,age,sep="/"))
                      ) +
       geom_line(aes(x=factor(sites),
                      y=cv,
                      col=factor(mod),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) +
       theme_minimal() + 
       labs(x = "Sites", y = "Clock Rate CV") + 
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       scale_color_grey(name="Data") +
       facet_wrap(~paste("Taxa:", tax))
ggsave("plots/sim_snp_cv.pdf", plt_snp_cv,w=3,h=3,units="in")

#plt_snp_ci_wid <- ggplot(cr_sum_snp) +
#       geom_point(aes(x=factor(sites),
#                      y=cv,
#                      col=factor(mod),
#                      shape=interaction(cr,age))
#                      ) +
#       geom_line(aes(x=factor(sites),
#                      y=cv,
#                      col=factor(mod),
#                      linetype=interaction(cr,age),
#                      group=interaction(tax,cr,age))) +
#       theme_minimal() + 
#       scale_color_grey() +
#       labs(x = "Sites", y = "Clock Rate/CI Width ") + 
#       facet_wrap(~tax)
#ggsave("plots/sim_snp_cci.pdf",plt_snp_ci_wid,w=3,h=3,units="in")

# plot 4: true tree posterior %
plt_snp_tr_topo <- ggplot(tree_sum_snp) +
       geom_point(aes(x=factor(sites),
                      y=tr_tree,
                      col=factor(mod),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                      y=tr_tree,
                      col=factor(mod),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) + 
       theme_minimal() + 
       scale_color_grey(name="Data") +
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       labs(x = "Sites", y = "True Tree Probability") + 
       facet_wrap(~paste("Taxa:", tax))
ggsave("plots/sim_snp_tr_topo_prob.pdf",plt_snp_tr_topo,w=3,h=3,units="in")

# plot 4: true tree posterior %
plt_snp_time <- ggplot(time_snp) +
       geom_point(aes(x=factor(sites),
                      y=comp_time,
                      col=factor(mod),
                      shape=interaction(cr,age,sep="/"))) +
       geom_line(aes(x=factor(sites),
                      y=comp_time,
                      col=factor(mod),
                      linetype=interaction(cr,age,sep="/"),
                      group=interaction(tax,cr,age))) + 
       theme_minimal() + 
       scale_color_grey(name="Data") +
       scale_shape_discrete(name = "CR/Age") +
       scale_linetype_discrete(name="CR/Age") +
       labs(x = "Sites", y = "Computation Time (s)") + 
       scale_y_log10() +
       facet_wrap(~paste("Taxa:", tax), scales="free_y")
ggsave("plots/sim_snp_time.png",plt_snp_time,w=3,h=3,units="in")

library(ggpubr)

print("arranging plots")
top <- ggarrange(
       plt_cr_mse,
       plt_cr_cov,  
       plt_cr_ci_wid, 
       plt_tr_topo_cr,
       plt_cr_bf, 
       ncol=5, common.legend=TRUE)
bottom <- ggarrange(
       plt_age_rmse,
       plt_age_cov,
       plt_age_ci_wid,
       plt_age_tr_topo,
       plt_age_bf, ncol=5, common.legend=TRUE) 

print("plot1")
plt_sim1 <- ggarrange(top, bottom,
   nrow=2, common.legend=TRUE)

print("plot2")
plt_sim2 <- ggarrange(
   plt_snp_rmse,
   plt_snp_cr_ci_cov,
   plt_snp_cv,
   plt_snp_tr_topo, 
   plt_snp_time,
   ncol=5, common.legend=TRUE)

ggsave("plots/plt_sim1_results.pdf", plt_sim1, w=12,h=4, units="in")
ggsave("plots/plt_sim2_results.pdf", plt_sim2, w=12,h=4, units="in")

