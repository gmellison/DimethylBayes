library(stringr)
library(ape)
library(tidysq)

options(scipen=999)
data_dir <- "data"
out_dir <- "sapelo"

nreps <- 1000000
burnin <- 50000

if (!dir.exists(out_dir)) dir.create(out_dir)

clones <- c("G", "R")
for (clone in clones) {
    parts <- list.files(sprintf("%s/%s",data_dir,clone), pattern="*.fa")
    parsed <- lapply (parts, function(f) {
        gene_name <- str_split_fixed(f,"\\.",2)[1]
        filepath <- sprintf("%s/%s/%s",data_dir,clone,f)
        fast <- tidysq::read_fasta(filepath,alphabet=c("0","1","2","?"))
    
        xx <- str_split_fixed(f, "[ag]", 3)[2]
        ntax <- length(fast$name) 
        
        chars <- fast$sq %>%
            tidysq::substitute_letters(c("0"="E","1"="M","2"="D")) %>%
            as.character
    
        nsite <- length(strsplit(chars[1],"")[[1]])
        charmat <- stringr::str_split_fixed(chars, "", n=nsite)
        nmissing <- sum(unlist(strsplit(chars, ""))=="?")
    
        # get mrate
        mmat <- match(charmat, c("?", "E", "M", "D"))
        mmat <- ifelse(mmat <= 2, 0, mmat)
        mmat <- ifelse(mmat == 0, 0, mmat-2)
        mrate <- sum(mmat) / (2 * nsite * ntax)
    
        list(ntax=ntax, nsite=nsite, chars=chars, gene_name=gene_name, xx=xx, nmissing=nmissing, name=fast$name, charmat=charmat, mrate=mrate)
    
    })
    
    charmat <- parsed[[1]]$mrate
    
    charset_idx <- 1
    for (i in 1:length(parsed)) {
    
        start_char <- charset_idx
        end_char <- start_char + parsed[[i]]$nsite - 1
        charset_idx <- end_char + 1
    
        parsed[[i]]$charset <- sprintf("charset %s = %s-%s;", parsed[[i]]$gene_name, start_char, end_char) 
        parsed[[i]]$matrix <- paste(parsed[[i]]$name, parsed[[i]]$chars, sep="      ")
    
    }
    
    # write a summary of the dataset:
    summ <- data.frame(name=sapply(parsed, function(x) x$gene_name),
               nsite=sapply(parsed, function(x) x$nsite),
               ntaxa=sapply(parsed, function(x) x$ntax),
               nmissing=sapply(parsed, function(x) x$nmissing),
               xx=sapply(parsed, function(x) x$xx),
               mrate=sapply(parsed, function(x) x$mrate))
    summ$missingrate <-with(summ, nmissing/ntaxa/nsite)
    write.csv(summ, sprintf("%s/%s/summary.csv", data_dir,clone))

    for (clock in c("","_rc","_sc")) {
        for (partition in c(1, 5, 10)) {

            nex_file <- sprintf("%s/%s_%s_%s.nex", out_dir, clone, partition, clock)
            nex_file <- str_replace(nex_file, "__", "_") 
            nex_file <- str_replace(nex_file, "_\\.", "\\.") 

            fc <- file(nex_file, "w")
            log_file <- sprintf("%s_%s_%s.log", clone, partition, clock)
            log_file <- str_replace(log_file, "__", "_") 
            log_file <- str_replace(log_file, "_\\.", "\\.") 

            # make partitions: 
            if (partition > 1) {
                qs <- quantile(summ$mrate, 1:partition/partition)
                cuts <- cut(summ$mrate, c(-0.1,qs))
                part <- match(cuts, levels(cuts))
                summ$mrate_part <- part
        
                # setup list to hold partitioned charset information
                partition_names <- paste0("p",1:partition) 
        
                # charset_idx <- 1
                part_charset <- lapply(partition_names, function(part) list(partition = part, charset=""))
        
                # loop through parsed list, add char ranges to charset list
                for (i in 1:length(parsed)) {
                    which_part <- summ$mrate_part[i]
                    gene_charset <- parsed[[i]]$charset
                    new_charset <- str_replace_all(str_split_fixed(gene_charset,"=", 2)[2],";","")
                    part_charset[[which_part]]$charset <- paste0(part_charset[[which_part]]$charset, new_charset, sep=" ")
                }
        
                for(i in 1:length(part_charset)) {
                    part_charset[[i]]$charset <- sprintf("charset %s = %s;", partition_names[i], part_charset[[i]]$charset)
                    part_charset[[i]]$charset_name <- partition_names[i]
                }
            }

            ## write matrix info
            writeLines(con=fc, 
                c(        "begin data;",
            	  sprintf("dimensions ntax=%s nchar=%s;", parsed[[1]]$ntax, charset_idx-1),
            	          "format datatype=dimethyl interleave=yes gap=- missing=?;", "", ""))
            
            # write the partitioned charsets
            writeLines(con=fc, " matrix")
            for (i in 1:length(parsed)) {
                writeLines(con=fc, 
                           parsed[[i]]$matrix)
            }
            writeLines(con=fc, "; \n")
            writeLines(con=fc, c("end;","",""))
           
            # write the charsets if partitioned:
            if (partition > 1) {
                writeLines(con=fc, "begin mrbayes;")
                for (i in 1:length(part_charset)) {
                    writeLines(con=fc, 
                               part_charset[[i]]$charset)
                }
                writeLines(con=fc, 
                           c("", sprintf("partition mrate_part = %s: %s;", 
                                             length(part_charset),
                                             paste0(sapply(part_charset, function(x) x$charset_name), collapse=","))))
                
                writeLines(con=fc, "end;")
            }
        
            # clone specific lines - only tree age is different 
            if (clone == "Finnish") {
                    tree_age_lines <- "prset treeagepr=lognormal(1000,500);"
            } else {
                    tree_age_lines <- "prset treeagepr=fixed(17);"
            }
   
            # clock specific lines 
            if (clock == "") {
                clock_rate_lines <- ""
            } else  {
                clock_rate_lines <- c("prset brlenspr=clock:uniform;")
                if ( clone == "Finnish" ) {
                    clock_rate_lines <- c(clock_rate_lines, "prset clockratepr=fixed(0.007);")
                } else {
                    clock_rate_lines <- c(clock_rate_lines, "prset clockratepr=normal(0.1,0.1);")
                }
            }
            if (clock == "_rc") clock_rate_lines <- c(clock_rate_lines, "prset clockvarpr=wn;")

            # if partitioned, set partition
            if (partition > 1) {
                partition_set <- "set partition=mrate_part;"
                unlink_lines <- c("unlink dimethylrate=(all);", 
                                  "prset applyto=(all) ratepr=variable;")

            } else {
                partition_set <- ""
                unlink_lines <- ""
            }

            # write the mrbayes model block
            mod_mcmc_lines <- c("begin mrbayes;",
               sprintf("log start filename=%s;",log_file),
               clock_rate_lines,
               tree_age_lines,
               partition_set,
               "lset rates=gamma;",
               "prset readerrpr=fixed(0.0009);",
               unlink_lines, 
               "report siterates=yes;",
               "mcmcp temp=0.3;",
               "mcmcp checkpoint=yes checkfreq=20000 append=no ;",
               sprintf("mcmc ngen=%s samplefreq=1000 burnin=%s;", nreps, burnin),
               sprintf("log stop filename=%s;",log_file),
               "end;")
            writeLines(con=fc, mod_mcmc_lines[mod_mcmc_lines!=""])
            close(fc)


            #### now write summary nex file
            sum_file <- sprintf("sum_%s_%s_%s.nex", clone, partition, clock)
            sum_file <- str_replace(sum_file, "__", "_") 
            sum_file <- str_replace(sum_file, "_\\.", "\\.") 
            fc <- file(sum_file,"w")

            log_file <- sprintf("sum_%s_%s_%s.log", clone, partition, clock)
            log_file <- str_replace(log_file, "__", "_") 
            log_file <- str_replace(log_file, "_\\.", "\\.") 


            writeLines(con=fc, 
                c(        "begin data;",
            	  sprintf("dimensions ntax=%s nchar=%s;", parsed[[1]]$ntax, parsed[[1]]$nsite),
            	          "format datatype=dimethyl interleave=yes gap=- missing=?;", "", ""))
            
            ## write matrix for tree reading 
            writeLines(con=fc, " matrix")
            writeLines(con=fc, parsed[[1]]$matrix)
            writeLines(con=fc, ";")
            writeLines(con=fc, c("end;",""))

            # write the model summary block
            writeLines(con=fc, 
            c("begin mrbayes;",
               sprintf("log start filename=%s;",log_file),
               sprintf("sump burninfrac=%s filename=%s;", nreps/burnin, nex_file),
               sprintf("sumt burninfrac=%s filename=%s;", nreps/burnin, nex_file),
               sprintf("log stop filename=%s;",log_file),
               "end;"))
            close(fc)

        }
    }
}

# list filenames 
fc <- file(sprintf("%s/files.lst", out_dir),"w")
for (clone in clones) {
    for (clock in c("","_rc","_sc")) {
        for (partition in c(1,5,10)) {
            out_file <- sprintf("%s_%s_%s.nex", clone, partition, clock)
            out_file <- str_replace(out_file, "__", "_") 
            out_file <- str_replace(out_file, "_\\.", "\\.") 
            writeLines(out_file, fc)
        }
    }
}
close(fc)

# list summary filenames 
fc <- file(sprintf("%s/files_sum.lst", out_dir),"w")
for (clone in clones) {
    for (clock in c("","_rc","_sc")) {
        for (partition in c(1,5,10)) {
            out_file <- sprintf("sum_%s_%s_%s.nex", clone, partition, clock)
            out_file <- str_replace(out_file, "__", "_") 
            out_file <- str_replace(out_file, "_\\.", "\\.") 
            writeLines(out_file, fc)
        }
    }
}
close(fc)

cluster_dir <- "dim"
jobs <- list( list(job_name="dim", 
                   num_files=18,
                   files_list="files.lst",
                   sub_file="sub_all.sh",
                   run_dir="${SLURM_JOB_NAME}"), 

              list(job_name="dim_sum",
                   num_files=18,
                   files_list="files_sum.lst",
                   sub_file="sum_all.sh",
                   run_dir="dim") )

## write submission scripts
for (job in jobs) {
    job_name <- job$job_name
    num_files <- job$num_files
    files_list <- job$files_list
    sub_file <- job$sub_file
    run_dir <- job$run_dir

    fc <- file(sprintf("%s/%s", out_dir, sub_file), "w")
    writeLines(
      c("#!/bin/bash",
        sprintf("#SBATCH --job-name=%s ", job_name),
        "#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)",
        "#SBATCH --ntasks=1		    	# 1 task (process) for below commands",
        "#SBATCH --cpus-per-task=1	 	# CPU core count per task, by default 1 CPU core per task",
        "#SBATCH --nodes=1",
        "#SBATCH --mem=16G		    	# Memory per node (4GB); by default using M as unit",
        "#SBATCH --time=5-00:00:00         # Time limit hrs:min:sec or days-hours:minutes:seconds",
        "#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out",
        "#SBATCH --mail-user=gme88782@uga.edu    # Where to send mail",
        "#SBATCH --mail-type=END          	# Mail events (BEGIN, END, FAIL, ALL)",
        "#SBATCH --constraint=EDR",
        sprintf("#SBATCH --array=1-%s", num_files),
        "",
        "# cd $SLURM_SUBMIT_DIR",
        "",
        sprintf("file=$(awk \"NR==${SLURM_ARRAY_TASK_ID}\" %s)", files_list),
        "",
        sprintf("# mkdir -p /scratch/${USER}/%s",run_dir),
        sprintf("cp $file /scratch/${USER}/%s",run_dir),
        sprintf("cd /scratch/${USER}/%s", run_dir),
        "",
        "./mb $file"),
      con=fc
    )
    close(fc)
}

## now copy files to sapelo and submit job

# copy files
# system(sprintf("cp ./mb %s", out_dir))
system(sprintf("scp -r %s/* gme88782@sapelo2.gacrc.uga.edu:/home/gme88782/dim/", out_dir))

# submit
cmd <- sprintf("ssh gme88782@sapelo2.gacrc.uga.edu \'cd /home/gme88782/%s/ ; sbatch %s \'", cluster_dir, "sub_all.sh")
system(cmd)





