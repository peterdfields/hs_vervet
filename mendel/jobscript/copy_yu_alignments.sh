#!/bin/bash
#PBS -N stagein
#PBS -P vervetmonkey
#PBS -q staging
#PBS -l select=2:ncpus=2:mpiprocs=4 -l place=scatter
#PBS -l walltime=12:00:00
module load mutil
# $PBS_O_WORKDIR contains the submission location (where qsub was called from)
cd /home/GMI/hannes.svardal/vervet_project/data/db/individual_alignment/
mpirun mcp -up --direct-read --direct-write --threads=1 --double-buffer --mpi /mnt/tmp/Hannes/individual_alignment/* .



