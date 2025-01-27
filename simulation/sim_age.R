######## 
# Clock Trees?
########
library(dplyr)
library(phangorn)
library(ape)
library(phybase)
library(TreeTools)
library(ggplot2)


source("functions.R")

out_dir <- "age"
tree_dir <- "trees"

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
"t20",
"t21",
"t22",
"t23",
"t24")

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
crs <- c(0.0001, 0.0007)
# snp_cr <- c(4*10^(-9))

a <- b <- 0.5
al <- 0.0
sites <- c(100, 1000, 10000)#, 50000)
taxa <- c(4,20)
bl <- 0.1
nsim <- 400 
ages <- c(20,100)

numruns <- nsim * length(crs) * length(sites) * length(taxa) * length(ages)

s <- 1
recook <- FALSE  

# set up the trees 
for (tax in taxa) {

    tree_file <- sprintf("%s/%s.%s.tree", tree_dir, tax, bl)

    # if tree exits, get it, else generate it
    if (file.exists(tree_file)) {
        tree <- read.tree(tree_file)
    } else {
        tree <- rtopology(tax, bl)    
        write.tree(tree, tree_file)
    }

    for (age in ages)  {
        for (cr in crs) {
            cr_tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
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
    for (cr in crs) {
        for (age in ages) {
            for (tax in taxa) {
                for (site in sites) {
                    print(sprintf("Running simulation %s / %s", s, numruns))

                    cr_line <- sprintf("prset clockratepr=normal(%s, %s);", cr, cr*0.1)
                    ta_line <- sprintf("prset treeagepr=lognormal(%s, %s);", age,age/4)
                    brlen <-   "prset brlenspr=clock:uniform;"

                    rates   <- "" # "lset rates=gamma;" 
                    model_lines <- c(cr_line, ta_line, brlen, rates)
                    mcmc <- c("mcmc nruns=1 ngen=100000 samplefreq=1000 burnin=20000;")

                    fname_methyl <- sprintf("%s.%s.%.1e.%s.%s.m",tax,age,cr,site,j)
                    cr_tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
                    tree_dm <- read.tree(cr_tree_f)

                    # check for necessary files in output directory
                    vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl))
                    tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl))
                    pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl))
                    lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl))
                    trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl))
                    contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl))
                    log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl))
                    nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl))

                    if ((!vstat | !tstat | !pstat | !lstat | !trprobs | !contre | !log | !nex) | 
                        recook) { 
                       print(sprintf("recooking %s", fname_methyl))
                       aln_methyl <- sim_alignment(tree_dm,a,b,site,al=0,dg=0)

                       setup_mrb(aln_methyl, tax, site, model_lines, mcmc, fname=fname_methyl, 
                                datatype="dimethyl", out_dir=out_dir)
                       run_mrb("./mb", sprintf("%s/%s.nex", out_dir, fname_methyl), sprintf("%s/%s.log", out_dir, fname_methyl))
                    }

                    # fit the non-clock model:
                    # if not exists, copy exising nex file, delete clock prior lines and run
                    fname_methyl_nc <- sprintf("%s.%s.%.1e.%s.%s.m.nc",tax,age,cr,site,j)

                    # check for proper outputs
                    vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl_nc))
                    tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl_nc))
                    pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl_nc))
                    lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl_nc))
                    trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl_nc))
                    contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl_nc))
                    log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl_nc))
                    nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl_nc))

                    if ((!vstat | !tstat | !pstat | !lstat | !trprobs | !contre | !log | !nex) | 
                        recook) { 
                        print(sprintf("recooking nc %s", fname_methyl_nc))
                        nex_c <- sprintf("%s/%s.nex", out_dir,fname_methyl)
                        nex_nc <- sprintf("%s/%s.nex", out_dir,fname_methyl_nc) 
                        system(sprintf("cp %s %s", nex_c, nex_nc))
                        fc <- file(nex_nc)
                        l <- readLines(fc)
                        l[stringr::str_detect(l, "prset")] <- ""
                        l <- stringr::str_replace_all(l, "\\.m", "\\.m\\.nc")
                        writeLines(l, fc)
                        close(fc)
                        run_mrb("./mb", nex_nc, sprintf("%s/%s.log", out_dir, fname_methyl_nc))
                    }
                    s <- s + 1
                }
            }
        }
    }
}

