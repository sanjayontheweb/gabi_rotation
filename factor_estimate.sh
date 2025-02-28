#!/bin/bash                                                                                               
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --mem=80G                                                                                         
#SBATCH --time=5-00:00:00                                                                                 
#SBATCH --exclude=c4-n20,c4-n26                                                                                  
#SBATCH --output=%x-%j.out                                                                                
#SBATCH --error=%x-%j.err                                                                                 
#SBATCH --gres=scratch:20G                                                                                

source /software/c4/cbi/software/miniforge3-24.11.0-0/etc/profile.d/conda.sh
conda activate spectra_kernel
python factor_estimate.py