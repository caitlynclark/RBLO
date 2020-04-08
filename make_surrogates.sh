#!/bin/bash -l
#SBATCH --ntasks=73                             # Request number of CPU cores
#SBATCH --time=4:00:00                          # Job should run for time
#SBATCH --account=windse                        # Accounting
#SBATCH --job-name=make_surrogates              # Job name
#SBATCH --mail-user caitlyn.clark@nrel.gov      # user email for notifcations
#SBATCH --mail-type ALL                         # ALL will notify for BEIGN,END,FAIL
#SBATCH --output=make_surrogates.%j.out         # %j will be replaced with job ID

srun --unbuffered -n 73 python3 -m mpi4py.futures example_make_surrogates.py

