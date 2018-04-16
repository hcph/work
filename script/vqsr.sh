 cat /zfssz3/MGI_BIT/RUO/raojunhua/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/shell/VVNWGSpilot-YH29/variant.0.1.vqsr.sh

      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/bcftools view -O z /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/input/variantion/genotype/VVNWGSpilot-YH29/VVNWGSpilot-YH29.GaeaGenotyper.sorted.vcf.gz --types snps > /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.snp.vcf.gz
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/tabix -f -p vcf /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.snp.vcf.gz
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/java -Xms8g -Xmx8g -jar /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hg19.fasta \
      --input /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.snp.vcf.gz \
      --mode SNP \
      --recal_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.snp.recal \
      --tranches_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.snp.tranches \
      -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hapmap_3.3.vcf.gz \
      -resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/1000G_omni2.5.vcf.gz \
      -resource:1000g,VCF,known=false,training=true,truth=false,prior=10.0 /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/1000G_phase1.snps.high_confidence.vcf.gz \
      -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/dbsnp-147.vcf.gz \
      -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 -tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 \
      -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 90.0 \
      -an DP -an QD -an FS -an ReadPosRankSum -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/java -Xms8g -Xmx8g -jar /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hg19.fasta \
      --input /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.snp.vcf.gz \
      --out /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/3.snp.vcf.gz \
      --recal_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.snp.recal \
      --tranches_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.snp.tranches \
      --mode SNP --ts_filter_level 99.0 -nt 4 \
      --disable_auto_index_creation_and_locking_when_reading_rods

      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/bcftools view -O z /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/input/variantion/genotype/VVNWGSpilot-YH29/VVNWGSpilot-YH29.GaeaGenotyper.sorted.vcf.gz --types indels > /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.indel.vcf.gz
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/tabix -f -p vcf /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.indel.vcf.gz
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/java -Xms8g -Xmx8g -jar /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hg19.fasta \
      --input /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.indel.vcf.gz \
      --mode INDEL \
      --recal_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.indel.recal \
      --tranches_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.indel.tranches \
      -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/Mills_and_1000G_gold_standard.indels.vcf.gz \
      --maxGaussians 4 \
      -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 -tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 \
      -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 90.0 \
      -an DP -an QD -an FS -an ReadPosRankSum -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment
      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/java -Xms8g -Xmx8g -jar /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hg19.fasta \
      --input /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/1.indel.vcf.gz \
      --out /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/3.indel.vcf.gz \
      --recal_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.indel.recal \
      --tranches_file /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/2.indel.tranches \
      --mode INDEL --ts_filter_level 99.0 -nt 4 \
      --disable_auto_index_creation_and_locking_when_reading_rods

      /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/java -Xms8g -Xmx8g -jar /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/lib/GenomeAnalysisTK.jar \
      -T CombineVariants \
      -R /zfssz3/MGI_BIT/RUO/raojunhua/00.private/software/software_self/20170210.Module_collection/20170601.wat/database/hg19/hg19.fasta \
      -o /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/4.vqsr.vcf.gz \
      --variant:v0 /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/3.snp.vcf.gz \
      --variant:v1 /hwfssz1/BIGDATA_COMPUTING/huweipeng/deligation/raojunhua/20170711.VV2_V3_PE100_PI_1_NA12878/20170727.stat/output/variant/VVNWGSpilot-YH29/vqsr/3.indel.vcf.gz \
      --rod_priority_list v0,v1 --genotypemergeoption PRIORITIZE --suppressCommandLineHeader --setKey null -nt 4
