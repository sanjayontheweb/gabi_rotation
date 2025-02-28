#!/bin/bash                                                                                               
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=4                                                                                      
#SBATCH --mem-per-cpu=33G                                                                                         
#SBATCH --time=5-00:00:00                                                                                 
#SBATCH --exclude=c4-n20,c4-n26                                                                                  
#SBATCH --output=%x-%j.out                                                                                
#SBATCH --error=%x-%j.err                                                                                 
#SBATCH --gres=scratch:20G                                                                                

source /software/c4/cbi/software/miniforge3-24.11.0-0/etc/profile.d/conda.sh
conda activate spectra_kernel
python spectra_skin_autoimmune.py --n_epochs 10000  --bulk_factors True --n_genes 9000
