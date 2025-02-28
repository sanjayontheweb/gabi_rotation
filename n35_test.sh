#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --nodelist=c4-n35
#SBATCH --output=%x.out

sleep 15
echo "this job is ${SLURM_JOB_ID} and running on ${SLURMD_NODENAME}"
