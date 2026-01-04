#! /bin/bash

# this will run 3 replicates of the simulation in the background
# using nohup. 

for ii in 1 2 3
  do
  nohup Rscript sim.R $ii & # run in background
  done

