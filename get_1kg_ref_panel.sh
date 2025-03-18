

# 1000 Genomes and C4 Reference Panel
# Convert multi-allelic variants to new variants
#-----------------------------------------------
ml bcftools tabix parallel
OUT=/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg
for CHR in $(seq 1 22)
do
FILE=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools norm -m - $FILE | bcftools view -O b -o $OUT/1kg_chr${CHR}_norm.bcf &
done

ls $OUT/*.bcf | parallel bcftools index

# Write MAP file of variant locations
cd /sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/
for CHR in $(seq 1 22)
do
	echo -e "CHROM\tPOS\tID\tA1\tA2" > 1kg_chr${CHR}_norm.map
	bcftools view -s HG00096 1kg_chr${CHR}_norm.bcf | grep -v "#" | cut -f1-5 >> 1kg_chr${CHR}_norm.map
done

# Filter for MAF in European samples
ml bcftools tabix parallel
OUT=/sc/arion/projects/roussp01a/gabriel/ref_panels/1kg
for CHR in $(seq 1 22)
do
FILE=/sc/arion/projects/data-ark/Public_Unrestricted/1000G/phase3/VCF/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools view -i "EUR_AF > 0.001" $FILE | bcftools norm -m - | bcftools view -O b -o $OUT/1kg_chr${CHR}_norm_eur.bcf &
done

ls $OUT/*_eur.bcf | parallel bcftools index

# Write MAP file of variant locations
cd /sc/arion/projects/roussp01a/gabriel/ref_panels/1kg/
for CHR in $(seq 1 22)
do
	echo -e "CHROM\tPOS\tID\tA1\tA2" > 1kg_chr${CHR}_norm_eur.map
	bcftools view -s HG00096 1kg_chr${CHR}_norm_eur.bcf | grep -v "#" | cut -f1-5 >> 1kg_chr${CHR}_norm_eur.map
done
