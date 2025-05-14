#!/bin/bash -l
# From http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
# 1000 Genomes 2500 + 700 relatives re-sequenced at NYGC to 30X. Aligned and jointly called on Hg38 reference build

cd /wrk/resources/genomes/1kGP_highcoverage_hg38/original 

for chr in chr{1..22}; do

echo "Downloading ${chr}..."
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

done

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf

# Filter for allele frequency, SNV only and drop individual genotypes
for chr in chr{1..22}; do

echo "Filtering ${chr}..."

bcftools view -G original/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
bcftools filter -i 'AF_EUR_unrel >= 0.005 && strlen(REF)=1 && strlen(ALT)=1' -O z > 1kg_highcoverage.${chr}.snv_filtered.vcf.gz

done

bcftools view -G original/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz | \
bcftools filter -i 'AF_EUR_unrel >= 0.005 && strlen(REF)=1 && strlen(ALT)=1' -O z > 1kg_highcoverage.chrX.snv_filtered.vcf.gz
