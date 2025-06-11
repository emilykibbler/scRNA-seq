# This script will be following this tutorial:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in

# It was run on Ohio Supercomputer Center using a combination of interactive mode and SLURM jobs where notated

# Create a conda environment
module load miniconda3

conda create -n scRNA
conda activate scRNA

mkdir scRNA-seq
cd scRNA-seq

# Module 0: install

mkdir apps
cd apps

wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1748422152&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=oXxScLUdV1EhBcO4hO5aBjLm87znVv1NtpX3LgGaytKjZqJvAvFHkt6GcePDkN02aAzCzZAlR47Wf2laD23vmySlJIlo1-mqP-Ap-aauv4~AePr-Ev~mSx1i~GChsIjqtOc8O1d~9gdGt1tjA-hrh55HAT9aYp5xfBxzGHGgQb8jAMm~pCSBgKgqW2RKQfgsyXzTfM1HRv253UZ-XQTjnRmDJ4Z7Jf5hHlqDeZLbC2UIyjU2vKtmBTMpUNRdEMoO-XM0mXK8a28zUCedHPytBsZ241JXoNgxvj0eesHvO1I1DAmkJZLt-FlACy3bbVm5BmJD1yRbGJ5Q6qzahRAKgw__"

# Check it's there
ls -1
# Unzip
tar -zxvf cellranger-9.0.1.tar.gz
# Navigate into it
cd cellranger-9.0.1

# Save the path to cellranger in the environment for easy access in the future
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
nano $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
export PATH=/users/PUOM0012/emilykibbler/scRNA-seq/apps/cellranger-9.0.1:$PATH
# then hit control+0 and enter, then control+x (control on a Mac not command)

# Confirm it's ready to go: this will print a version if successful up to here
which cellranger

cellranger
# Collect system specs
cellranger sitecheck > sitecheck.txt

# Print all your system specs
# Check 10X website to see if compatible (if using OSC, more than enough)
less sitecheck.txt
q # to quit file reader mode

cellranger testrun --id=check_install

# Module 1: mkfastq (demultiplex to create fastq files)

# Need an Illumina demultiplexer called bcl2fastq

# This file will need to be requested from Illumina by creating an account on their website
unzip bcl2fastq2-v2-20-0-linux-x86-64.zip

env -i rpm2cpio ./bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm | cpio -idmv

mkdir run_cellranger_mkfastq
cd run_cellranger_mkfastq

wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -zxvf cellranger-tiny-bcl-1.2.0.tar.gz

cat cellranger-tiny-bcl-simple-1.2.0.csv

tree -L 2 cellranger-tiny-bcl-1.2.0/

# Display help page; confirms successful installation
cellranger mkfastq --help

# This makes bcl2fastq accessible
export PATH=/users/PUOM0012/emilykibbler/scRNA-seq/usr/local/bin:$PATH

cellranger mkfastq --id=tutorial_walk_through \
  --run=./cellranger-tiny-bcl-1.2.0 \
  --csv=./cellranger-tiny-bcl-simple-1.2.0.csv

# Navigate into the folder created by previous step
cd ./tutorial_walk_through/outs/fastq_path
# Checks to make sure it worked
ls -1
ls -1 H35KCBCXY/test_sample

# Module 2: Cell Ranger Count
cd ~
cd scRNA-seq

mkdir run_cellranger_count
cd run_cellranger_count

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar

wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
tar -zxvf refdata-gex-GRCm39-2024-A

# Check if command is available by showing help page
cellranger count --help

# Cellranger count will likely take more than the time limit on interactive mode

#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/scRNA-seq
module load miniconda3
conda activate scRNA

cellranger count --id=run_count_1kpbmcs \
   --fastqs=./run_cellranger_count/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=./run_cellranger_count/refdata-gex-GRCh38-2024-A \
   --localcores=48 \
   --create-bam false

