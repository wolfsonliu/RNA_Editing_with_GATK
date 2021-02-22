#! /bin/bash
#SBATCH --job-name=re_3_va
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
# 3. variant
####################

# GATK HaplotypeCaller
echo "[$(date)] GATK HaplotypeCaller..."
gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R ${GENOME_FA} \
    -I ${DIR_2_BAM}/${lab}.rgadd.dedup.split.recal.bam \
    -O ${DIR_3_VA}/${lab}.vcf.gz \
    -D ${dbSNP} \
    --minimum-mapping-quality 0 \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 0
echo "[$(date)] GATK HaplotypeCaller finished."

# Filter
# select SNP
echo "[$(date)] GATK SelectVariants SNP..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    --select-type-to-include SNP \
    -V ${DIR_3_VA}/${lab}.vcf.gz \
    -O ${DIR_3_VA}/${lab}.snp.vcf.gz
tabix -f ${DIR_3_VA}/${lab}.snp.vcf.gz
echo "[$(date)] GATK SelectVariants SNP finished."

echo "[$(date)] GATK VariantFiltration..."
# gatk VariantFiltration \
#     -R ${GENOME_FA} \
#     -V ${DIR_3_VA}/${lab}.vcf.gz \
#     -O ${DIR_3_VA}/${lab}.cluster.vcf.gz \
#     --cluster-window-size 35 --cluster-size 3
#      --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
#      --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
#      --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
#      --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
#      --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
#      --filter-expression "SOR > 3.0" --filter-name "SOR_gt_3"
gatk VariantFiltration \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.snp.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filter.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
    --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
    --filter-expression "MQ < 30.0" --filter-name "MQ_lt_30" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
    --filter-expression "DP < 20.0" --filter-name "DP_lt_20" \
    --filter-expression "QUAL < 20.0" --filter-name "QUAL_lt_20"
tabix -f ${DIR_3_VA}/${lab}.filter.vcf.gz
echo "[$(date)] GATK VariantFiltration finished."

# select SNP
echo "[$(date)] GATK SelectVariants rm filter..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.filter.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filtered.vcf.gz \
    --exclude-filtered true
tabix -f ${DIR_3_VA}/${lab}.filtered.vcf.gz
echo "[$(date)] GATK SelectVariants rm filter finished."

echo "[$(date)] GATK SelectVariants rm dbSNP..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.filtered.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filtered.rmdbsnp.vcf.gz \
    --discordance ${dbSNP}
tabix -f ${DIR_3_VA}/${lab}.filtered.rmdbsnp.vcf.gz
echo "[$(date)] GATK SelectVariants rm dbSNP finished."

echo "[$(date)] GATK SelectVariants rm 1000G..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.filtered.rmdbsnp.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.vcf.gz \
    --discordance ${G1000}
tabix -f ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.vcf.gz
echo "[$(date)] GATK SelectVariants rm 1000G finished."

echo "[$(date)] GATK SelectVariants rm EVS..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.rmevs.vcf.gz \
    --discordance ${EVS}
tabix -f ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.rmevs.vcf.gz
echo "[$(date)] GATK SelectVariants rm EVS finished."

echo "[$(date)] GATK SelectVariants rm Alu..."
gatk SelectVariants \
    -R ${GENOME_FA} \
    -V ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.rmevs.vcf.gz \
    -O ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.rmevs.rmalu.vcf.gz \
    -XL ${ALU}
tabix -f ${DIR_3_VA}/${lab}.filtered.rmdbsnp.rm1000g.rmevs.rmalu.vcf.gz
echo "[$(date)] GATK SelectVariants rm Alu finished."
####################
