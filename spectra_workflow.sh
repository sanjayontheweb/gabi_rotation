#!/bin/bash                                                                                               
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --mem=1G                                                                                         
#SBATCH --time=1-00:00:00                                                                                 
#SBATCH --exclude=c4-n20                                                                                  
#SBATCH --output=%x-%j.out                                                                                
#SBATCH --error=%x-%j.err                                                                                 
#SBATCH --gres=scratch:20G 


job_output=$(sbatch --partition=freecycle,krummellab,common,dscolab -x c4-n20 run_spectra.sh)
job_id=$(echo $job_output | awk '{print $4}')

log_stats_job_output=$(sbatch --dependency=afterok:$job_id --partition=freecycle,krummellab,common,dscolab -x c4-n20 log_job_stats.sh "$job_id")
