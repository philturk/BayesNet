#!/bin/bash
##
## simulation.slurm.sh: submits an array of jobs, each job is one simulation
## create a Jobs directory before running
##
#SBATCH --job-name=simulation                      #job name
#SBATCH --output=Jobs/simulation-%j.out            #output file (with job ID %j)
#SBATCH --error=Jobs/simulation-%j.err             #error file (with job ID %j)
#SBATCH --partition=defq                           #partition of Hyalite to run the job on (generally set to defq, has 24 hour time limit, 62 nodes, 32+ cores, 64+ GB of memory)
#SBATCH --nodes=1                                  #number of nodes to allocate (keep to 1 for higher priority in Hyalite queue?)
#SBATCH --ntasks-per-node=1                        #number of descrete tasks - keep at one except for MPI (??)
#SBATCH --cpus-per-task=1                          #number of cores to allocate per job (generally set to 8 to 32, which is 0.3 to 1 node)
#SBATCH --mem-per-cpu=2000                         #system memory to use (per core) (generally set at 2GB/core = 2000MB)
#SBATCH --time=23:59:00                            #maximum job run time (dd-hh:mm:ss, generally set at 20% longer than expect, cannot exceed 24 hours on defq)
#SBATCH --array=1-50                               #indices for number of jobs in array (SLURM_ARRAY_TASK_ID)
#SBATCH --mail-user=sallyslipher@montana.edu       #send emails to this address
#SBATCH --mail-type=FAIL                           #email on when a job fails (options are ALL, BEGIN, END, FAIL, REQUEUE)

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID ##pastes SLURM_ARRAY_TASK_ID: followed by array index number

module load R/3.5.3 ##load R into the environment (3.5.3 is latest version on Hyalite at the time)

Rscript --vanilla BayesInf_v4_sim_SKS.R $SLURM_ARRAY_TASK_ID ##run R script, first input is task ID from environment
