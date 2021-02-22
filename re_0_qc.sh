#! /bin/bash
#SBATCH --job-name=re_0_qc
#SBATCH --output=job.%x.%j_%A_%a.%N.array.out
#SBATCH --error=job.%x.%j_%A_%a.%N.array.err
#SBATCH --partition=C032M0256G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 1
#SBATCH --cpu-freq=high
#SBATCH --time=120:00:00

set -e -u -o pipefail

thread=1

source settings.sh

lab=$(awk -F "\t" 'FNR == '$SLURM_ARRAY_TASK_ID'{print $1;}' ${FQs})
rfqz1=$(awk -F "\t" 'FNR == '$SLURM_ARRAY_TASK_ID'{print $2;}' ${FQs})
rfqz2=$(awk -F "\t" 'FNR == '$SLURM_ARRAY_TASK_ID'{print $3;}' ${FQs})

####################
# 0. qc
####################

zcat ${rfqz1} > ${DIR_0_QC}/${lab}.1.fq
zcat ${rfqz2} > ${DIR_0_QC}/${lab}.2.fq

echo "[$(date)] 1st FastQC..."
fastqc -t ${thread} -f fastq -o ${DIR_0_QC} \
    ${DIR_0_QC}/${lab}.1.fq ${DIR_0_QC}/${lab}.2.fq
echo "[$(date)] 1st FastQC finished."

forward=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
backward=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
echo "[$(date)] cutadapt..."
cutadapt -a ${forward} -A ${backward} \
    -u 6 \
    -q 10,20 --minimum-length 90 \
    -o ${DIR_0_QC}/clean-${lab}.1.fq \
    -p ${DIR_0_QC}/clean-${lab}.2.fq \
    ${DIR_0_QC}/${lab}.1.fq ${DIR_0_QC}/${lab}.2.fq
echo "[$(date)] cutadapt finished."

echo "[$(date)] 2nd FastQC..."
fastqc -t ${thread} -f fastq -o ${DIR_0_QC} \
    ${DIR_0_QC}/clean-${lab}.1.fq ${DIR_0_QC}/clean-${lab}.2.fq
echo "[$(date)] 2nd FastQC finished."

rm -f ${DIR_0_QC}/${lab}.1.fq
rm -f ${DIR_0_QC}/${lab}.2.fq
####################
