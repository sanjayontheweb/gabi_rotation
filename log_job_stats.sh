#!/bin/bash                                                                                               
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --mem=2G                                                                                         
#SBATCH --time=1-00:00:00                                                                                 
#SBATCH --exclude=c4-n20                                                                                  
#SBATCH --output=%x-%j.out                                                                                
#SBATCH --error=%x-%j.err                                                                                 
#SBATCH --gres=scratch:20G 

# Check if job ID is provided as an argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide the previous job ID"
    exit 1
fi

job_id=$1

RESULTS_FILE="spectra_job_stats.csv"

if [ ! -f "$RESULTS_FILE" ]; then
    echo "jobid,n_samples,highly_variable,n_epochs,eigen_factor_used,wall_clock_time,memory_used,node,cpu_efficiency" > "$RESULTS_FILE"
fi

# Get job ID, wall clock time, and memory used from seff
seff_output=$(/software/c4/utilities/seff $job_id)
wall_clock_time=$(echo "$seff_output" | grep "Job Wall-clock time:" | awk '{print $4}')
memory_used=$(echo "$seff_output" | grep "Memory Utilized:" | awk '{print $3}')
cpu_efficiency=$(echo "$seff_output" | grep "CPU Efficiency:" | awk '{print $3}')
node=$(sacct -j $job_id --format=NodeList --noheader | awk 'NR==1 {print $1}')

# Get metadata from Python JSON output

# Wait for the Python script to finish
if [ ! -f "job_metadata.json" ]; then
    echo "Error: Results JSON file not found!"
    exit 1
fi

# Read JSON file
n_sample_cells=$(jq -r '.n_sample_cells' job_metadata.json)
highly_variable=$(jq -r '.highly_variable' job_metadata.json)
n_epochs=$(jq -r '.n_epochs' job_metadata.json)
bulk_eigen=$(jq -r '.bulk_factors' job_metadata.json)

echo "$job_id,$n_sample_cells,$highly_variable,$n_epochs,$bulk_eigen,$wall_clock_time,$memory_used,$node,$cpu_efficiency" >> "$RESULTS_FILE"

# Create a directory to store all .err and .out files, named logs_{job_id}
log_dir="logs_${job_id}"
mkdir -p "$log_dir"

# Move all corresponding .err and .out files into the folder
mv *"$((job_id - 1))".err *"$((job_id - 1))".out \
   *"${job_id}".err *"${job_id}".out \
   *"$((job_id + 1))".err *"$((job_id + 1))".out \
   "$log_dir" 2>/dev/null