#!/bin/bash

date=$(date "+%Y_%m_%d")
echo $date

module load R
module load singularity

mkdir -p slurm_out


snakemake --jobs 200 --use-conda \
--use-singularity \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -p {cluster.p} -o {cluster.o} --mem {cluster.mem} --time {cluster.time} --job-name {cluster.name}"
