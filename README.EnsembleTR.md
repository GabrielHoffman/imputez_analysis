



# Get EnsembleTR:
 Paper: [(https://www.nature.com/articles/s41467-023-42278-3](https://www.nature.com/articles/s41467-023-42278-3)
 
 GitHub with VCFs: [https://github.com/gymrek-lab/EnsembleTR](https://github.com/gymrek-lab/EnsembleTR)

## Process VCFs
```sh
cd /hpc/users/hoffmg01/www/imputez_analysis/EnsembleTR/data

# hg38
seq 1 22 | parallel -P4 wget https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v4/ensembletr_refpanel_v4_chr{}.vcf.gz
 
ml tabix bcftools
ls *.vcf.gz | parallel tabix -p vcf {}

# normalize multi-allelic TR's
for FILE in $(ls *.vcf.gz); do
  OUT=norm/$(basename $FILE .vcf.gz).norm.bcf
  bcftools norm -m - $FILE | bcftools view -O b -o $OUT &
done

# index
ls norm/*bcf | parallel bcftools index {}

# map of variant locations
for FILE in $(ls norm/*bcf); do
  OUT=norm/$(basename $FILE .bcf).map.gz
  echo -e "CHROM\tPOS\tID\tA1\tA2" | gzip > $OUT
  bcftools view -s HG00096 $FILE | zgrep -v "#" | cut -f1-5 | gzip >> $OUT &
done
```


# Brain eQTL data 
Summary statistics: [https://www.synapse.org/Synapse:syn25592272](https://www.synapse.org/Synapse:syn25592272)


## Process summary statistics
Write to parquet format:

```R
file = "/hpc/users/hoffmg01/www/imputez_analysis/mmQTL_brain_meta_eqtl_all.tsv"
df_qtl = read_tsv(file) %>%      
          filter(eQTL_order == 1)

file = "/hpc/users/hoffmg01/www/imputez_analysis/mmQTL_brain_meta_eqtl_all.parquet"
df = df_qtl %>%
  dplyr::rename(CHROM=Chr, ID = Variant, z = z_score_fixed) %>%
  select(CHROM, Gene, ID, z)

unique(df$ID) %>%
  snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, ., ifnotfound="drop") %>%
  data.frame %>%
  tibble %>%
  dplyr::rename(CHROM = seqnames, 
                POS = pos,
                ID = RefSNP_id) %>%
  mutate(CHROM = paste0('chr', CHROM)) %>%
  select(-strand, -alleles_as_ambig) %>%
  inner_join(df, by=c("CHROM", "ID")) %>%
  write_parquet( file )
```

## Run imputation on tandom repeats
```sh
Rscript run_EnsembleTR_imputation.R
```