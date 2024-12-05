
######## 
# 
########
library(dplyr)
library(phangorn)
library(ape)
library(phybase)
library(TreeTools)
# library(ggplot2)

print("running snp comp")
source("functions.R")

out_dir <- "snp_comp"
tree_dir <- "trees"

if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(tree_dir)) dir.create(tree_dir)

tax_names <- c(
"t01",
"t02",
"t03",
"t04",
"t05",
"t06",
"t07",
"t08",
"t09",
"t10",
"t11",
"t12",
"t13",
"t14",
"t15",
"t16",
"t17",
"t18",
"t19",
"t20")

make_rc_tree <- function(tree, cr, age, clock="relaxed", model="iln", cv=0.1) {
    tree <- chronos(tree, model="relaxed",
                    calibration=makeChronosCalib(tree, age.max=age, age.min=age))
    rc_tree <- simclock::relaxed.tree(tree,model,r=cr,s2=cv)
    return(rc_tree)
}

make_sc_tree <- function(tree, cr, age) {
    tree <- chronos(tree,
                    calibration=makeChronosCalib(tree, age.max=age, age.min=age))
    tree$edge.length <- tree$edge.length * cr
    return(tree)
}

# methyl & dna clock rates
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
nsim <- 400 
ages <- c(20)

numruns <- nsim * length(sites) * length(taxa) * length(clocks) * length(ages)

s <- 1
recook <- FALSE 

# make the trees, if !exist yet 
par(mfrow=c(1,length(clocks)))
for (tax in taxa) {

    tree_file <- sprintf("%s/%s.%s.tree", tree_dir, tax, bl)
    
    # if tree exits, get it, else generate it
    if (file.exists(tree_file)) {
        tree <- read.tree(tree_file)
    } else {
        tree <- rtopology(tax, bl)    
        write.tree(tree, tree_file)
    }

    for (age in ages) {
        for (clock in clocks) {
            mod <- clock$mod
            cr <- clock$cr

            cr_tree_f <- sprintf("%s/%s.%s.%s.tree", tree_dir, tax, cr, age)
            print(cr_tree_f)
            # if tree exits, get it, else generate it
            if (file.exists(cr_tree_f)) {
                tree_dm <- read.tree(cr_tree_f)
            } else {
                tree_dm <- make_sc_tree(tree, cr, age)    
                tree_dm$tip.label <- tax_names[1:tax]
                write.tree(tree_dm, cr_tree_f)
            }
        }
    }
}

for (j in 1:nsim) {
    for (clock in clocks) {
        cr <- clock$cr
        mod <- clock$mod
        for (age in ages) {
    
            for (site in sites) {
                for (tax in taxa) {
                    print(sprintf("Running simulation %s / %s", s, numruns))

                    cr_line <- sprintf("prset clockratepr=normal(%s, %s);",0.1,0.1)
                    ta_line <- sprintf("prset treeagepr=fixed(%s);", age)
                    brlen <-   "prset brlenspr=clock:uniform;"
                    rates   <- "" # "lset rates=gamma;" 
                    model_lines <- c(cr_line, ta_line, brlen, rates)
                    mcmc <- c("mcmc ngen=50000 samplefreq=500 burnin=20000;")

                    if (mod == "m") {
                        fname_methyl <- sprintf("%s.%s.%.1e.%s.%s.m",tax,age,cr,site,j)
                        cr_tree_f <- sprintf("%s/%s.%s.%s.tree", tree_dir, tax, cr, age)
                        tree_dm <- read.tree(cr_tree_f)
                        if (!file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl)) | recook) { 
                            if (site > 50000) {
                                intlv <- "yes"
                                site_per_chunk <- 10000
                                num_chunks <- site/site_per_chunk
                                if (site %% site_per_chunk != 0 ) stop("site > 50000 should be in a multiple of 10000")
                                sims <- lapply(1:num_chunks, function(chunk) sim_alignment(tree_dm,a,b,site_per_chunk,al=0,dg=0))
                                aln_methyl <- unlist(sims)
                            } else {
                                intlv <- "no"
                                aln_methyl <- sim_alignment(tree_dm,a,b,site,al=0,dg=0)
                            }

                            setup_mrb(aln_methyl, tax, sprintf("%d",site), model_lines, mcmc, fname=fname_methyl, 
                                    datatype="dimethyl", out_dir=out_dir, interleave=intlv)
                            run_mrb("./mb", sprintf("%s/%s.nex", out_dir, fname_methyl), sprintf("%s/%s.log", out_dir, fname_methyl))
                        }  
                    }

                    if (mod == "dna") {
                        fname_dna <- sprintf("%s.%s.%.1e.%s.%s.dna",tax,age,cr,site,j)
                        cr_tree_f <- sprintf("%s/%s.%s.%s.tree", tree_dir, tax, cr, age)
                        tree_dna  <- read.tree(cr_tree_f)

                        if (!file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_dna)) | recook) { 
                            if (site > 50000) {
                                intlv <- "yes"
                                site_per_chunk <- 10000
                                num_chunks <- site/site_per_chunk
                                if (site %% site_per_chunk != 0 ) stop("site > 50000 should be in a multiple of 10000")
                                sims <- lapply(1:num_chunks, function(chunk) {
                                    aln_dna <- sim_dna(tree_dna,site_per_chunk,rep(1,6),al=0)
                                    aln_dna_str <- paste(tax_names[1:tax], aln_dna$aln_str, sep="   ")
                                    })
                                aln_dna_str <- unlist(sims)
                            } else {
                                intlv <- "no"
                                aln_dna <- sim_dna(tree_dna,site,rep(1,6),al=0)
                                aln_dna_str <- paste(tax_names[1:tax], aln_dna$aln_str, sep="   ")
                            }

                            # put a prior on the clockrate, assume tree age is known
                            setup_mrb(aln_dna_str, tax, sprintf("%d",site), model_lines, mcmc, 
                                      fname=fname_dna, datatype="dna", out_dir=out_dir,
                                      interleave=intlv) 
                            run_mrb("./mb", sprintf("%s/%s.nex", out_dir, fname_dna), sprintf("%s/%s.log", out_dir, fname_dna))
                        } 
                    }
                    s <- s + 1
                }
            }
        }
    }
}

