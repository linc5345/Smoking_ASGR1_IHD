library(pacman)
# Packages needed for perform_mr() function
p_load(rio, tidyverse, data.table, MendelianRandomization,
       testthat, grid,gridExtra,ieugwasr,patchwork,TwoSampleMR,gsmr)


##########################################################################################################################################
# Ensure 'perform_mr()' function is defined - see 'https://github.com/adamkvonende/mr_pipeline' for details incl. column name requirements
##########################################################################################################################################

#####User inputs (please edit these paths as necessary)#####

# File should contain protein assay names and gene coords
# Required columns: Assay, CHR, Gene_start, Gene_end
smok_prots <- fread("path/to/protein_list_with_coords.csv")


# Disease variable outcomes
IHD  <- fread("path/to/IHD.csv")
T2D  <- fread("path/to/T2D.csv")
COPD <- fread("path/to/COPD.csv")
IS   <- fread("path/to/IS.csv")
LC   <- fread("path/to/LC.csv")

# Directory where protein summary statistics are stored
prot_dir <- "/path/to/olink/proteins/" 

# Output directory
out_dir  <- "/path/to/your/results/deposit/"

# Output filenames
results_file <- file.path(out_dir, "prot_to_disease_results.csv")
het_file     <- file.path(out_dir, "prot_to_disease_het.csv")

# outcome list used in loop (can be modified)
outcomes <- list(IHD = IHD, T2D = T2D, COPD = COPD, IS = IS, LC = LC)

# Helper function to find protein file by assay name
find_prot_file <- function(assay, dir) {
  # Try exact filename containing assay (case-insensitive), prefer .gz
  patt <- paste0("(?i)", assay) # case-insensitive
  files <- list.files(dir, full.names = TRUE)
  cand <- files[grepl(patt, files)]
  if (length(cand) == 0L) return(NA_character_)
  # prefer gz
  gz <- cand[grepl("\\.gz$", cand)]
  if (length(gz) > 0L) return(gz[1L])
  return(cand[1L])
}

# For varying gene windows, modify to use case
flank_kb = 1000

# Safe column canonicalisation (toggle to your dataset)
canonicalise_outcome_cols <- function(DT) {
  # mapping of common headers -> canonical names expected by perform_mr()
  old <- c("CHR","chr","BP38","pos","EA_FREQ","EAF","BETA","beta","SE","se","P","p","EA","effect_allele","AA","other_allele","SNP","ID","rsid")
  new <- c("chr","chr","pos","pos","eaf","eaf","beta","beta","se","se","pval","pval","effect_allele","effect_allele","other_allele","other_allele","rsid","rsid","rsid")
  # perform mapping only for present columns (match first occurrence for duplicates)
  present <- intersect(old, names(DT))
  if (length(present) == 0L) return(DT)
  # for each present, pick the corresponding new name (first match)
  for (p in present) {
    idx <- which(old == p)[1]
    newname <- new[idx]
    # avoid renaming if already the canonical name
    if (p != newname && !(newname %in% names(DT))) {
      setnames(DT, p, newname)
    }
  }
  return(DT)
}


# Append per-call using data.table::fwrite(..., append = TRUE)
write_results <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(invisible())
  if (!file.exists(results_file)) {
    fwrite(as.data.table(df), results_file)
  } else {
    fwrite(as.data.table(df), results_file, append = TRUE)
  }
}
write_het <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(invisible())
  if (!file.exists(het_file)) {
    fwrite(as.data.table(df), het_file)
  } else {
    fwrite(as.data.table(df), het_file, append = TRUE)
  }
}


# Main loop: iterate proteins
n_prots <- nrow(smok_prots)
for (i in seq_len(n_prots)) {
  row <- smok_prots[i]
  assay <- as.character(row$Assay)
  chr   <- as.character(row$CHR)
  # parse potential composite coords (semicolon-separated. e.g., AMY1A_AMY1B_AMY1C) - take min start, max end
  parse_vals <- function(x) {
    x <- as.character(x)
    parts <- unlist(strsplit(x, ";", fixed = TRUE))
    nums <- suppressWarnings(as.numeric(parts))
    nums <- nums[!is.na(nums)]
    if (length(nums) == 0L) return(NA_real_)
    return(nums)
  }
  svals <- parse_vals(row$Gene_start)
  evals <- parse_vals(row$Gene_end)
  gene_start <- if (!is.na(svals[1])) min(svals) else NA_real_
  gene_end   <- if (!is.na(evals[1])) max(evals) else NA_real_
  if (is.na(gene_start) || is.na(gene_end)) {
    warning(sprintf("Protein %s: invalid gene_start/gene_end; skipping", assay))
    next
  }
  # apply flank
  gene_start_flank <- max(1L, as.integer(gene_start - flank_kb * 1000L))
  gene_end_flank   <- as.integer(gene_end + flank_kb * 1000L)
  
  # find protein file
  prot_file <- find_prot_file(assay, prot_dir)
  if (is.na(prot_file)) {
    warning(sprintf("Protein file not found for assay '%s' in %s; skipping", assay, prot_dir))
    next
  }
  
  # read file
  ss_exposure <- tryCatch({
    fread(prot_file)
  }, error = function(e) {
    warning(sprintf("Failed to read protein file %s : %s", prot_file, e$message))
    return(NULL)
  })
  if (is.null(ss_exposure)) next
  
  # canonicalise common column names
  ss_exposure <- canonicalise_outcome_cols(ss_exposure)
  
  # perform MR for each smoking outcome
  for (out_name in names(outcomes)) {
    outcome_obj <- outcomes[[out_name]]
    
    res <- tryCatch({
      perform_mr(
        ss_exposure   = ss_exposure,
        ss_outcome    = outcome_obj,
        cis           = TRUE,
        gene_chr      = as.numeric(chr),
        gene_start    = gene_start_flank,
        gene_end      = gene_end_flank,
        flank_kb      = flank_kb,
        maf_thresh    = 0.01,
        p_thresh      = 1.7e-11,
        harmo_action  = 1,
        clump_when    = "after",
        clump_locally = FALSE,
        clump_r2      = 0.001,
        clump_pop     = "EUR",
        name_exp      = assay,
        name_outcome  = out_name,
        which_methods = "basic",
        plots         = FALSE,
        outcome_binary= TRUE
      )
    }, error = function(e) {
      warning(sprintf("perform_mr failed for %s -> %s: %s", assay, out_name, e$message))
      return(NULL)
    })
    
    if (is.null(res)) next
    
    # results and het: attach metadata then write
    if (!is.null(res$results) && nrow(as.data.table(res$results)) > 0L) {
      dt_res <- as.data.table(res$results)
      # ensure useful metadata present
      dt_res[, `:=` (outcome = out_name, exposure = assay, chr = chr,
                     gene_start = gene_start_flank, gene_end = gene_end_flank)]
      # write (append)
      write_results(dt_res)
    }
    
    if (!is.null(res$het) && nrow(as.data.table(res$het)) > 0L) {
      dt_het <- as.data.table(res$het)
      dt_het[, `:=` (outcome = out_name, exposure = assay, chr = chr,
                     gene_start = gene_start_flank, gene_end = gene_end_flank)]
      write_het(dt_het)
    }
  } # outcomes loop
  
  # cleanup
  rm(ss_exposure); gc()
} # proteins loop