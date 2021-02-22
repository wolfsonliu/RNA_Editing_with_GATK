# Directory
DIR_WK=$(pwd)
DIR_FQ=${DIR_WK}/fqdata

DIR_0_QC=${DIR_WK}/0_qc && mkdir -p ${DIR_0_QC}
DIR_1_MAP=${DIR_WK}/1_map && mkdir -p ${DIR_1_MAP}
DIR_2_BAM=${DIR_WK}/2_bam && mkdir -p ${DIR_2_BAM}
DIR_3_VA=${DIR_WK}/3_variant && mkdir -p ${DIR_3_VA}

# Data
# <label>\t<fq1.gz>\t<fq2.gz>
FQs=${DIR_WK}/fqs.txt

# Reference hg38
STAR_REF_GENOME=${DIR_WK}/Reference/STARIndex
GENOME_FA=${DIR_WK}/Reference/genome.fa
dbSNP=${DIR_WK}/Reference/dbSNP/dbsnp.38.unsorted.vcf
G1000=${DIR_WK}/Reference/1000genome/hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
EVS=${DIR_WK}/Reference/esp/EVS.new.vcf.gz
ALU=${DIR_WK}/Reference/ucsc_repeatmasker_alu_hg38_filtered.bed
