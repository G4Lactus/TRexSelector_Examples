#!/bin/bash
#SBATCH --job-name=trex_gvs_mc
#SBATCH -A p0020087
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12              # 12 parallel threads per file
#SBATCH --mem-per-cpu=31000M            # 31,000 MB (~31 GB) per core * 12 = 372 GB total RAM
#SBATCH -C avx512                       # Force Phase 2 nodes
#SBATCH --array=1-12                    # Array job for files 01 to 12
#SBATCH --time=00-23:00:00              # 23 hours time limit
#SBATCH -D /home/uo31suqa/trex_software_tests  # Working directory

# Save error and output files to our structured logs folder
#SBATCH -o logs/out_Job_%A_%02a.txt
#SBATCH -e logs/err_Job_%A_%02a.txt

# Enable e-mail notification
#SBATCH --mail-user=fabian.scheidt@tu-darmstadt.de
#SBATCH --mail-type=begin,end

# Load required modules for the job
module purge
module load gcc/11.5.0-z7mc
module load openmpi/4.1.8-6xzv
module load r/4.5.2-uxdk

# Set locale
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Working directory
cd ~/trex_software_tests

# ---------------------------------------------------------
# Dynamic File Mapping using SLURM Array ID
# ---------------------------------------------------------
# Format the array task ID to have a leading zero (1 -> 01, 12 -> 12)
FILE_IDX=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

# Define the script to run based on the formatted ID
SCRIPT2Run="src/simulation_trex_gvs_${FILE_IDX}.R"

# Submit the script
srun Rscript "${SCRIPT2Run}" "${SLURM_CPUS_PER_TASK}"
