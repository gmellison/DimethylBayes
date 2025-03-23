######## 
# Clock Trees?
########
library(dplyr)
library(phangorn)
library(ape)
library(phybase)
library(TreeTools)
library(ggplot2)

j <- 1
arg <- commandArgs(trailingOnly=TRUE)
j <- arg[1]

source("functions.R")

out_dir <- sprintf("sim_out/%s", j)
print(sprintf("output going into %s", out_dir))
if (!dir.exists(out_dir)) dir.create(out_dir)
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


make_sc_tree <- function(tree, cr, age) {
    tree <- chronos(tree,
                    calibration=makeChronosCalib(tree, age.max=age, age.min=age))
    tree$edge.length <- tree$edge.length * cr
    return(tree)
}

# methyl & dna clock rates
crs <- c(0.0001, 0.0007)
# snp_cr <- c(4*10^(-9))

a <- 0.6
b <- 1-a
al <- 0.0
#sites <- c(100, 1000, 10000, 50000)
sites <- c(50000)
taxa <- c(4,20)
bl <- 0.1
nsim <- 1 
ages <- c(20)

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
                #tree_dm <- read.tree(cr_tree_f)
                #pdf(sprintf("plots/tree%s.pdf", i))
                #par(xpd=NA)
                #plot(tree_dm)
                #edgelabels(round(tree_dm$edge.length/cr,3), frame="none", adj=c(0,-1))
                #dev.off()
                #i <- i+1
            } else {
                tree_dm <- make_sc_tree(tree, cr, age)    
                tree_dm$tip.label <- tax_names[1:tax]
                write.tree(tree_dm, cr_tree_f)
            }
        }
    }
}
#for (j in 1:nsim) {
for (cr in crs) {
    for (age in ages) {
        for (tax in taxa) {
            for (site in sites) {
            
                sim_settings <- sprintf("%s.%s.%.1e.%s.%s",tax,age,cr,site,j)
                print(sprintf("Running simulation %s / %s", s, numruns))
                s <- s + 1

                fname_methyl_age <- sprintf("%s.%s.%.1e.%s.%s.m.age",tax,age,cr,site,j)
                fname_methyl_cr <- sprintf("%s.%s.%.1e.%s.%s.m.cr",tax,age,cr,site,j)

                fname_methyl_nc_age <- sprintf("%s.%s.%.1e.%s.%s.m.age.nc",tax,age,cr,site,j)
                fname_methyl_nc_cr <- sprintf("%s.%s.%.1e.%s.%s.m.cr.nc",tax,age,cr,site,j)

                run1 <- run2 <- run3 <- run4 <- TRUE

                # global
                cr_tree_f <- sprintf("%s/%s.%s.%.1e.%s.tree", tree_dir, tax, bl, cr, age)
                tree_dm <- read.tree(cr_tree_f)
                # run for necessary files in output directory
                vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl_age))
                tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl_age))
                pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl_age))
                lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl_age))
                trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl_age))
                contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl_age))
                log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl_age))
                nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl_age))

                if (vstat & tstat & pstat & lstat & trprobs & contre & log & nex) run1 <- FALSE 

                # run for proper outputs
                vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl_nc_age))
                tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl_nc_age))
                pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl_nc_age))
                lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl_nc_age))
                trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl_nc_age))
                contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl_nc_age))
                log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl_nc_age))
                nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl_nc_age))

                if (vstat & tstat & pstat & lstat & trprobs & contre & log & nex) run2 <- FALSE 

                vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl_cr))
                tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl_cr))
                pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl_cr))
                lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl_cr))
                trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl_cr))
                contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl_cr))
                log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl_cr))
                nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl_cr))

                if (vstat & tstat & pstat & lstat & trprobs & contre & log & nex) run3 <- FALSE 

                vstat   <- file.exists(sprintf("%s/sumt.%s.vstat", out_dir, fname_methyl_nc_cr))
                tstat   <- file.exists(sprintf("%s/sumt.%s.tstat", out_dir, fname_methyl_nc_cr))
                pstat   <- file.exists(sprintf("%s/sump.%s.pstat", out_dir, fname_methyl_nc_cr))
                lstat   <- file.exists(sprintf("%s/sump.%s.lstat", out_dir, fname_methyl_nc_cr))
                trprobs <- file.exists(sprintf("%s/sumt.%s.trprobs", out_dir, fname_methyl_nc_cr))
                contre  <- file.exists(sprintf("%s/sumt.%s.con.tre", out_dir, fname_methyl_nc_cr))
                log     <- file.exists(sprintf("%s/%s.log", out_dir, fname_methyl_nc_cr))
                nex     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_methyl_nc_cr))

                if (vstat & tstat & pstat & lstat & trprobs & contre & log & nex) run4 <- FALSE 

                # if we need to run anything, get the alignment
                if (any(c(run1,run2,run3,run4))) {
                    aln_file <- sprintf("sim_out/alns/aln%s.txt",j)
                    if (!file.exists(aln_file)) {
                        aln_methyl <- sim_alignment(tree_dm,a,b,site,al=0,dg=0)
                        writeLines(aln_methyl, aln_file)
                    } else {
                        aln_methyl <- readLines(aln_file)
                    }
                }

                if (run1) {
                    print (sprintf("running %s", fname_methyl_age))
                    ## first fix clock rate ##
                    cr_line <- sprintf("prset clockratepr=fixed(%s);",cr)
                    ta_line <- sprintf("prset treeagepr=lognormal(%s, %s);", age,age/2)
                    brlen <-   "prset brlenspr=clock:uniform;"

                    rates   <- "" # "lset rates=gamma;" 
                    model_lines <- c(cr_line, ta_line, brlen, rates)
                    mcmc <- c("mcmc nruns=1 ngen=400000 samplefreq=1000 burnin=400000;")
                    # run fixed clock rate sim
                    setup_mrb(aln_methyl, tax, site, model_lines, mcmc, fname=fname_methyl_age, 
                                datatype="dimethyl", out_dir=out_dir, burninfrac=0.5)
                    run_mrb("./mb", sprintf("%s/%s.nex", out_dir, fname_methyl_age), sprintf("%s/%s.log", out_dir, fname_methyl_age))

                }
                # fit the non-clock model:
                # if not exists, copy exising nex file, delete clock prior lines and run

                if (run2) {
                    print(sprintf("running %s", fname_methyl_nc_age))
                    nex_c <- sprintf("%s/%s.nex", out_dir,fname_methyl_age)
                    nex_nc <- sprintf("%s/%s.nex", out_dir,fname_methyl_nc_age) 
                    system(sprintf("cp %s %s", nex_c, nex_nc))
                    fc <- file(nex_nc)
                    l <- readLines(fc)
                    l[stringr::str_detect(l, "prset")] <- ""
                    l <- stringr::str_replace_all(l, "\\.m\\.age", "\\.m\\.age\\.nc")
                    writeLines(l, fc)
                    close(fc)
                    run_mrb("./mb", nex_nc, sprintf("%s/%s.log", out_dir, fname_methyl_nc_age))
                }

                if (run3) {
                ## now fix tree age ##
                    cr_line <- sprintf("prset clockratepr=normal(%s,%s);", cr, cr*2)
                    ta_line <- "prset treeagepr=fixed(20);"
                    brlen <-   "prset brlenspr=clock:uniform;"

                    rates   <- "" # "lset rates=gamma;" 
                    model_lines <- c(cr_line, ta_line, brlen, rates)
                    mcmc <- c("mcmc nruns=1 ngen=400000 samplefreq=1000;")

                    # check for necessary files in output directory
                    setup_mrb(aln_methyl, tax, site, model_lines, mcmc, fname=fname_methyl_cr, 
                                datatype="dimethyl", out_dir=out_dir, burninfrac=0.5)
                    run_mrb("./mb", sprintf("%s/%s.nex", out_dir, fname_methyl_cr), sprintf("%s/%s.log", out_dir, fname_methyl_cr))
                }

                # fit the non-clock model:
                # if not exists, copy exising nex file, delete clock prior lines and run

                if (run4) {
                # check for proper outputs
                    print(sprintf("%s/%s.nex", out_dir,fname_methyl_cr))
                    nex_c <- sprintf("%s/%s.nex", out_dir,fname_methyl_cr)
                    nex_nc <- sprintf("%s/%s.nex", out_dir,fname_methyl_nc_cr) 
                    system(sprintf("cp %s %s", nex_c, nex_nc))
                    fc <- file(nex_nc)
                    l <- readLines(fc)
                    l[stringr::str_detect(l, "prset")] <- ""
                    l <- stringr::str_replace_all(l, "\\.m\\.cr", "\\.m\\.cr\\.nc")
                    writeLines(l, fc)
                    close(fc)
                    run_mrb("./mb", nex_nc, sprintf("%s/%s.log", out_dir, fname_methyl_nc_cr))
                }

            }
        }
    }
}
#}

