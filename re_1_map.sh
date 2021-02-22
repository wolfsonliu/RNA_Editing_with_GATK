#! /bin/bash
#SBATCH --job-name=re_1_map
#SBATCH --output=job.%x.%j_%A_%a.%N.array.out
#SBATCH --error=job.%x.%j_%A_%a.%N.array.err
#SBATCH --partition=C032M0256G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 32
#SBATCH --cpu-freq=high
#SBATCH --time=120:00:00

set -e -u -o pipefail

thread=32

source settings.sh

lab=$(awk -F "\t" 'FNR == '$SLURM_ARRAY_TASK_ID'{print $1;}' ${FQs})

####################
# 1. mapping with STAR
####################

echo "[$(date)] STAR..."
STAR --genomeDir ${STAR_REF_GENOME} \
    --runThreadN ${thread} --readFilesCommand cat \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${DIR_0_QC}/clean-${lab}.1.fq ${DIR_0_QC}/clean-${lab}.2.fq \
    --outFileNamePrefix ${DIR_1_MAP}/${lab}.
echo "[$(date)] STAR finished."

gzip ${DIR_0_QC}/clean-${lab}.1.fq
gzip ${DIR_0_QC}/clean-${lab}.2.fq
####################
