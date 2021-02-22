#! /bin/bash
#SBATCH --job-name=re_2_bam
#SBATCH --output=job.%x.%j_%A_%a.%N.array.out
#SBATCH --error=job.%x.%j_%A_%a.%N.array.err
#SBATCH --partition=C032M0256G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -N 1
#SBATCH --cpu-freq=high
#SBATCH --time=120:00:00

set -e -u -o pipefail

thread=1

source settings.sh

lab=$(awk -F "\t" 'FNR == '$SLURM_ARRAY_TASK_ID'{print $1;}' ${FQs})

####################
# 2. bam modification
####################

# Picard AddOrReplaceReadGroups
echo "[$(date)] Picard AddOrReplaceReadGroups..."
java -jar $HOME/.local/bin/picard.jar \
    AddOrReplaceReadGroups \
    I=${DIR_1_MAP}/${lab}.Aligned.sortedByCoord.out.bam \
    O=${DIR_2_BAM}/${lab}.rgadd.bam \
    SO=coordinate \
    RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
echo "[$(date)] Picard AddOrReplaceReadGroups finished."

# Picard MarkDuplicates
echo "[$(date)] Picard MarkDuplicates..."
java -jar $HOME/.local/bin/picard.jar MarkDuplicates \
    I=${DIR_2_BAM}/${lab}.rgadd.bam \
    O=${DIR_2_BAM}/${lab}.rgadd.dedup.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=${DIR_2_BAM}/${lab}.rgadd.dedup.metrics
echo "[$(date)] Picard MarkDuplicates finished."

# GATK SplitNCigarReads
echo "[$(date)] GATK SplitNCigarReads..."
gatk SplitNCigarReads \
    -R ${GENOME_FA} \
    -I ${DIR_2_BAM}/${lab}.rgadd.dedup.bam \
    -O ${DIR_2_BAM}/${lab}.rgadd.dedup.split.bam
echo "[$(date)] GATK SplitNCigarReads finished."

# GATK BaseRecalibrator
echo "[$(date)] GATK BaseRecalibrator..."
gatk  BaseRecalibrator \
    -R ${GENOME_FA} \
    --known-sites ${dbSNP} \
    -I ${DIR_2_BAM}/${lab}.rgadd.dedup.split.bam \
    -O ${DIR_2_BAM}/${lab}.rgadd.dedup.split.recal_gatk4.grv
echo "[$(date)] GATK BaseRecalibrator finished."

echo "[$(date)] GATK ApplyBQSR..."
gatk ApplyBQSR \
    -R ${GENOME_FA} \
    --bqsr-recal-file ${DIR_2_BAM}/${lab}.rgadd.dedup.split.recal_gatk4.grv \
    -I ${DIR_2_BAM}/${lab}.rgadd.dedup.split.bam \
    -O ${DIR_2_BAM}/${lab}.rgadd.dedup.split.recal.bam
echo "[$(date)] GATK ApplyBQSR finished."
####################
