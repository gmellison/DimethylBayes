# usage: 
# Rscript sim.R i
# where i is an integer indicating the simulation replicate index.
# the script will create the trees if they don't exist (in the trees/ directory and write the 
# simulation files into the sim/i directory. 

library(phangorn)
library(ape)

source("functions.R")

# get simulation replicate number from the command line
j <- commandArgs(trailingOnly=TRUE)[1]

sim_dir <- "sim"
tree_dir <- "trees"
out_dir <- sprintf("%s/%s", sim_dir, j)

if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive=TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
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

# simulation params
cr <- 0.001
a <- 0.6
b <- 1-a
al <- 0.0
sites <- c(100, 1000, 10000)
taxa <- c(5,10,15,20)
age <- c(20)

recook <- FALSE  
sigma <- 0.2375
ilnvar <- sigma^2
cr <- 0.001


print("making trees")
# set up the trees 
for (tax in taxa) {

    sc_tree_file  <- sprintf("%s/%s.sc.tree", tree_dir, tax)
    iln_tree_file <- sprintf("%s/%s.iln.tree", tree_dir, tax)

    # if tree exits, get it, else generate it
    if (file.exists(iln_tree_file)) {
        tree_sc <- read.tree(sc_tree_file)
        tree_iln <- read.tree(iln_tree_file)
    } else {
        tree_dm <- rtree(tax)    

        # age of tree is 20:
        cal <- makeChronosCalib(tree_dm, node="root", age.min=20, age.max=20)
        tree_dm <- chronos(tree_dm, calibration = cal)
        tree_dm$tip.label <- tax_names[1:tax]
        tree_time <- tree_dm

        # "strict clock tree"
        tree_sc <- tree_time
        tree_sc$edge.length <- tree_time$edge.length * cr

        # use iln mutation rates to modify branch lengths
        tree_iln <- tree_sc
        bl_mrates <- rlnorm(length(tree_sc$edge.length), -log(ilnvar + 1)/2, sqrt(log(ilnvar + 1)))
        tree_iln$edge.length <- tree_iln$edge.length * bl_mrates
        write.tree(tree_sc, sc_tree_file)
        write.tree(tree_iln, iln_tree_file)
    }
}

clocks <- c("iln", "none")

for (site in sites) {
    for (tax in taxa) {
        iln_tree_file <- sprintf("%s/%s.iln.tree", tree_dir, tax)
        tree_iln  <- read.tree(iln_tree_file)
        aln_file <- sprintf("%s/%s.%s.aln.txt",out_dir,tax,site)

        # simulate alignment with some readerr
        if (!file.exists(aln_file)) {
            aln_methyl2 <- sim_alignment(tree_iln,a,b,site,al=0,dg=0,readerr=0.005)
            save(aln_methyl2,file=aln_file)
        } else {
            load(aln_file)
        }

        for (clock in clocks) {

                # filenames for param est and ss runs
                fname <- sprintf("%s.%s.%s",tax,clock,site)
                fname_ss <- sprintf("%s.%s.%s.ss",tax,clock,site)

                print(sprintf("Running %s/%s ", out_dir, fname))

                # clock model priors
                brlen_pr <- "prset brlenspr=clock:uniform;"
                cr_pr    <- sprintf("prset clockratepr=normal(%s, %s);", cr, cr*10)
                ta_pr    <- sprintf("prset treeagepr=fixed(%s);", age)
                clvar_pr <- sprintf("prset clockvarpr=%s;", clock)

                # unconstrained bl priors
                if (clock == "none") {
                    cr_pr     <- ""
                    ta_pr     <- ""
                    brlen_pr  <- ""
                    clvar_pr  <- ""
                }

                rates   <- "" # "lset rates=gamma;" # not using gamma rates for sim
                model_lines <- c(brlen_pr, ta_pr, cr_pr, clvar_pr, "prset readerrpr = uniform(0,0.05); ")
                #model_lines <- c(brlen_pr, ta_pr, cr_pr, clvar_pr, "prset readerrpr = fixed(0.005); ")
                mcmc <- c("mcmc nruns=2 ngen=500000 samplefreq=200 burnin=300000;")
                ss <- c("ss nruns=1 ngen=2000000;")

                # check for necessary files in output directory
                vstat   <- file.exists(sprintf("%s/%s.vstat", out_dir, fname))
                pstat   <- file.exists(sprintf("%s/%s.pstat", out_dir, fname))
                lstat   <- file.exists(sprintf("%s/%s.lstat", out_dir, fname))
                log     <- file.exists(sprintf("%s/%s.log", out_dir, fname))

                # run if the output files don't exist. 
                # only run the iln model for parameter estimation
                if (((!vstat | !pstat | !lstat | !log) | recook) & clock == "iln") {  

                   print(sprintf("recooking %s", fname))
                   # write nexus file and run binary
                   setup_mrb(aln_methyl2, tax, site, model_lines, mcmc, fname=fname, 
                           datatype="dimethyl", out_dir=out_dir)
                   print(sprintf("./dmb %s > %s", sprintf("%s/%s.nex", out_dir, fname), sprintf("%s/%s.log", out_dir, fname)))
                   run_mrb("./dmb", sprintf("%s/%s.nex", out_dir, fname), sprintf("%s/%s.log", out_dir, fname))
                }
                
                # run the ss files 
                # check for existing output
                lstat_ss   <- file.exists(sprintf("%s/%s.lstat", out_dir, fname_ss))
                nex_ss     <- file.exists(sprintf("%s/%s.nex", out_dir, fname_ss))
                log_ss     <- file.exists(sprintf("%s/%s.log", out_dir, fname_ss))

                if ((!lstat_ss | !log_ss | !nex_ss) | recook) { 
                    print(sprintf("recooking %s", fname_ss))
                 
                    # write nexus file and run binary for ss run
                    setup_mrb(aln_methyl2, tax, site, model_lines, ss, fname=fname_ss, 
                            datatype="dimethyl", out_dir=out_dir)
                    run_mrb("../src/dmb", sprintf("%s/%s.nex", out_dir, fname_ss), sprintf("%s/%s.log", out_dir, fname_ss))
                }
        }
    }
}

