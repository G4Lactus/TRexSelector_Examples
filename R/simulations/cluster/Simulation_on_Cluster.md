# Cluster Usage

## 1. Directory Architecture

Log into the Lichtenberg II cluster and create a dedicated, structured environment in your home 
directory.
We separate the R scripts, SLURM submission scripts, and log files.

Command on the login node to generate the structure:

```bash
mkdir -p ~/trex_software_tests/{src,slurm,logs,data}
```

- `src/`: Holds all your `.R` script files.
- `slurm/`: Holds the corresponding `.slurm` batch scripts.
- `logs/`: SLURM will automatically pipe your console outputs here.
- `data/`: Used for any small dataset dependencies (if using massive datasets, create a symlink here pointing to `$HPC_SCRATCH`).


## 2. Syncing Your Local Files

From your **local** terminal (not the cluster), use `rsync` to push the R scripts you have already created into the cluster's `src` directory. Replace `<TU-ID>` with your actual ID:

```bash
rsync -avz /home/*my_id_doesn_matter*/trex_software_tests/ <TU-ID>@lcluster1.hrz.tu-darmstadt.de:~/trex_software_tests/src/
```

```bash
rsync -avz /home/*my_id_doesn_matter*/trex_software_tests/xyz_{02..12}.R <TU-ID>@lcluster1.hrz.tu-darmstadt.de:~/trex_software_tests/src/
```

## 3. The SLURM R Template

To capture your console outputs exactly as you requested, we will use SLURM's native standard output redirection (`--output`). When you run `Rscript`, all `print()` and `cat()` statements from the `TRexSelector` package will be written directly to a `.txt` or `.out` file in your `logs/` directory.

Here is the baseline template we will adapt for each of your files:

```bash
#!/bin/bash
#SBATCH --job-name=trex_gvs_base
#SBATCH --output=logs/%x_%j.txt       # Saves console output as JobName_JobID.txt
#SBATCH --time=02:00:00               # Adjust based on expected runtime
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4             # Adjust if your R script uses parallel processing
#SBATCH --mem=16G                     # Adjust memory based on dataset size

# 1. Clean environment and load R
module purge
module load R                         # You may need to append a version, e.g., R/4.3.1

# 2. Execute the script
# NOTE: Always submit this script from the ~/trex_software_tests/ directory
Rscript src/YOUR_SCRIPT_NAME.R
```
