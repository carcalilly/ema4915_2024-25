#!/bin/bash
#SBATCH --partition=normal                   #Partition/Queue name
#SBATCH --ntasks-per-node=48                 #number of CPU cores.
#SBATCH --job-name=vasp_std                  #Job Name
#SBATCH --nodes=1                            #Number of nodes to use
#SBATCH --time=48:00:00                      #Maximum amount of time to run
#SBATCH --output=slurm.%j.out                #File to send standard output
#SBATCH --error=slurm.%j.err                 #File to send standard error

ulimit -s unlimited
module load vasp/vasp-6.4.2-WITH_VTST-oneapi-2023.1.0
srun vasp_std