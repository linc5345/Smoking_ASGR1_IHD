# Load necessary packages
library(data.table)
library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(MVMR)

# Define output directory
out_dir = "/path/to/results/depository"

# Load in smoking variable, BMI, and protein outcome
exp <- fread("/path/to/your/instrument/exp.csv") # Expected col: rsid, beta, se, eaf, effect_allele, other_allele, pval, N, chr, pos 
bmi <- fread("/path/to/bmi.csv")                 # Expected col: rsid, beta, se, eaf, effect_allele, other_allele, pval, N, chr, pos
prot <- fread("path/to/protein.csv")             # Expected col: rsid, beta, se, eaf, effect_allele, other_allele, pval, N, chr, pos

# Format datasets for MVMR
exp_df <- as.data.frame(exp)
exp_df$TRAIT <- "SMOKING"
exp_mvmr <- format_data(exp_df,
                        type = "exposure",
                        header = TRUE,
                        phenotype_col = "TRAIT",
                        snp_col = "rsid",
                        beta_col = "beta",
                        se_col = "se",
                        eaf_col = "eaf",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        pval_col = "pval",
                        samplesize_col = "N",
                        chr_col = "chr",
                        pos_col = "pos",
                        log_pval = FALSE) %>% as_tibble()

bmi_df <- as.data.frame(bmi)
bmi_df$TRAIT <- "BMI"
bmi_mvmr <- format_data(bmi_df,
                        type = "exposure",
                        header = TRUE,
                        phenotype_col = "TRAIT",
                        snp_col = "rsid",
                        beta_col = "beta",
                        se_col = "se",
                        eaf_col = "eaf",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        pval_col = "pval",
                        samplesize_col = "N",
                        chr_col = "chr",
                        pos_col = "pos",
                        log_pval = FALSE) %>% as_tibble()

prot_df <- as.data.frame(prot)
prot_df$TRAIT <- "PROTEIN"
prot_mvmr <- format_data(prot_df,
                         type = "outcome",
                         header = TRUE,
                         phenotype_col = "TRAIT",
                         snp_col = "rsid",
                         beta_col = "beta",
                         se_col = "se",
                         eaf_col = "eaf",
                         effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele",
                         pval_col = "pval",
                         chr_col = "chr",
                         pos_col = "pos",
                         log_pval = FALSE) %>% as_tibble()

# Clump each exposure
clump_args <- list(clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
clump_threshold <- 5e-8

clump_and_keep_gws <- function(df) {
  tmp <- df %>% filter(pval.exposure < 1e-6)
  if (nrow(tmp) == 0) return(tmp)
  clp <- clump_data(tmp,
                    clump_kb = clump_args$clump_kb,
                    clump_r2 = clump_args$clump_r2,
                    clump_p1 = clump_args$clump_p1,
                    clump_p2 = clump_args$clump_p2,
                    pop = clump_args$pop)
  clp %>% filter(pval.exposure < clump_threshold)
}

exp_clumped <- clump_and_keep_gws(exp_mvmr)
bmi_clumped <- clump_and_keep_gws(bmi_mvmr)

if (nrow(exp_clumped) == 0) stop("No genome-wide significant SNPs found for the exposure after clumping.")
if (nrow(bmi_clumped) == 0) stop("No genome-wide significant SNPs found for BMI after clumping.")

# Combined SNP list and second-round clump for independence
mvmr_snps <- unique(c(exp_clumped$SNP, bmi_clumped$SNP))

comb_exp <- bind_rows(
  exp_mvmr %>% filter(SNP %in% mvmr_snps),
  bmi_mvmr %>% filter(SNP %in% mvmr_snps)
)

# second round clump across union
# TwoSampleMR::clump_data needs pval.exposure column and SNP
clump_input <- comb_exp %>% distinct(SNP, pval.exposure)
clumped_union <- clump_data(clump_input,
                            clump_kb = clump_args$clump_kb,
                            clump_r2 = clump_args$clump_r2,
                            clump_p1 = clump_args$clump_p1,
                            clump_p2 = clump_args$clump_p2,
                            pop = clump_args$pop)

comb_exp_independent <- comb_exp %>% filter(SNP %in% clumped_union$SNP)

# Ensure SNPs exist in outcome (protein)
outcome_clump <- semi_join(prot_mvmr, comb_exp_independent, by = "SNP")
missing_snps <- setdiff(comb_exp_independent$SNP, outcome_clump$SNP)
if (length(missing_snps) > 0) {
  message(length(missing_snps), " exposure SNP(s) not present in protein outcome and will be dropped.")
}

# Harmonise
mvdat <- mv_harmonise_data(comb_exp_independent, outcome_clump)

# Run MVMR
mvmr_res <- mv_multiple(mvdat, plots = FALSE)
mvmr_res$result <- mvmr_res$result %>%
  mutate(lci = b - 1.96 * se, uci = b + 1.96 * se)

# Save MVMR results
out_res_file <- file.path(config$out_dir, "mvmr_result.csv")
write.csv(mvmr_res$result, out_res_file, row.names = FALSE)