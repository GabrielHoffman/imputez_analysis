
cat("Loading packages...\n")
suppressPackageStartupMessages({
library(GenomicDataStream)
library(imputez)
library(tidyverse)
library(arrow)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(qtlPlots)
library(remaCor)
library(arrow)
library(ggbio)
})

# Compute -log10 P from z-statistic
# stable for large z-statistics
z_to_score = function(z){

  # -log10 P from z-statistic
  # -log10(2*pnorm(abs(z.stat), lower.tail=FALSE))
  -(pnorm(abs(z), lower.tail=FALSE, log.p=TRUE) + log(2)) / log(10)
}

ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

setwd("/hpc/users/hoffmg01/www/imputez_analysis/EnsembleTR/data/")


cat("Reading data...\n")

# Target Genes:
genes = c("HTT", "C9orf72", "TCF4", "CNNM4", "TSPAN14", "HRCT1")

# get ENSEMBL ids
df_gene = ensembldb::select(x = ensdb, keys = genes, keytype = 'SYMBOL', column = c('GENEID', "SEQNAME"))

# European samples
ids.euro = read.table("EUR.ids")$V1

# Read eQTL results
file = "/hpc/users/hoffmg01/www/imputez_analysis/mmQTL_brain_meta_eqtl_all.parquet"
df_qtl = read_parquet(file)

# for each target gene
for( gene in genes){

  cat(gene, "\n")

  # Read reference panel
  ######################

  # get chrom number
  chrom = df_gene %>% 
    filter(SYMBOL == gene) %>%
    pull(SEQNAME)

  # ENSEMBL id
  ensID = df_gene %>% 
    filter(SYMBOL == gene) %>%
    pull(GENEID)

  # Read reference panel
  file = paste0("norm/ensembletr_refpanel_v4_chr", chrom, ".norm.bcf")
  gds.1kg = GenomicDataStream(file, "GT", initialize=TRUE, MAF=0.01, samples=ids.euro)

  # Variant locations in referene panel
  file = paste0("norm/ensembletr_refpanel_v4_chr", chrom, ".norm.map.gz")
  df_map = read_tsv(file, progress=FALSE, show_col_types=FALSE)
  colnames(df_map) = c("CHROM", "POS", "ID", "REF", "ALT")

  # Extract gene from QTL results
  df2 = df_qtl %>%
          filter(Gene == ensID) %>%
    inner_join(df_map, by=c("CHROM", "POS")) %>%
    dplyr::rename(ID = ID.y, REF_A1 = REF, REF_A2 = ALT) %>%
    mutate(GWAS_A1 = REF_A1, GWAS_A2 = REF_A2)

  # run imputation
  ################
  res = run_imputez(df2, gds.1kg, 100000, 10000, quiet=TRUE)

  # merge with variant info
  res = res %>%
    dplyr::rename(REF = A1, ALT = A2) %>%
    left_join(df_map, by=c("ID", "REF", "ALT"))

  # Multi-allelics: fixed effects meta-analysis
  ##############################################
  tab = table(res$ID)

  res_meta = lapply(names(tab)[tab>1], function(id){

    cat(id, "\n")
    df.focus = res %>%
              filter(ID == id)

    region = df_map %>%
      filter(ID == id) %>%
      with(paste0(CHROM[1], ":", min(POS), "-", max(POS)))

    gds = setRegion(gds.1kg, region)
    dat <- getNextChunk(gds)

    colnames(dat$X) = dat$info$A2
    C = cor(dat$X)
    C = C[df.focus$ALT, df.focus$ALT]

    # random effects meta-analysis
    z.value = RE2C(df.focus$z.stat, df.focus$se, C) %>% 
              mutate(z=qnorm(RE2Cp/2, lower.tail=FALSE)) %>% 
              pull(z)

    tibble(ID = id, 
      A1 = df.focus$REF[1], 
      A2 = "meta", 
      z.stat = z.value, 
      se = 1, 
      r2.pred = NA, 
      lambda = NA, 
      CHROM = gsub(":.*", "", region),
      POS = as.numeric(gsub(".*:(\\d+)-.*", "\\1", region))
      )
  })
  res_meta = do.call(rbind, res_meta)

  # Combine Observed and Imputed z-statistics
  res_combined = bind_rows(res_meta, res) %>%
                  mutate(imputed = TRUE) %>%
                  bind_rows(df2 %>% 
                    mutate(imputed=FALSE) %>%
                    dplyr::rename(z.stat = z)) %>%
                  filter(!is.na(z.stat)) %>%
                  mutate( score = z_to_score(z.stat))

  file = paste0("/hpc/users/hoffmg01/www/imputez_analysis/EnsembleTR/results/", gene, ".RDS")
  saveRDS(res_combined, file=file)

  # Convert to GRanges
  gr = res_combined %>% 
      arrange(imputed) %>%
      with(GRanges(CHROM, IRanges(POS, POS), score=score, name=ID, inCandidateSet=imputed))

  wh = GRanges(paste0("chr", chrom), 
    IRanges(start(range(gr)), end(range(gr))))
  wh = resize(wh, width=6e5, fix="center")

  fig_mht = plotMht(gr, wh, recombRate=FALSE, size=12, ptSize=2)

  fig_gene = plotEnsGenes( ensdb, wh, size=20)

  fig_track = ggbio::tracks(
       gene     = fig_mht,
      'Genes'    = fig_gene,
      padding = unit(-.65, "lines"), 
      label.bg.fill="navy", 
      label.text.color="white",
      heights=c(1, .1),
      title = gene) 
  names(fig_track@labeled)[1] = gene
  fig_track@mutable['Genes'] = FALSE
  # fig_track@fixed[] = TRUE
  xlim(fig_track) = c(start(wh), end(wh))
  
  file = paste0("/hpc/users/hoffmg01/www/imputez_analysis/EnsembleTR/figures/", gene, ".pdf")
  pdf(file)
  fig_track
  dev.off()
}

