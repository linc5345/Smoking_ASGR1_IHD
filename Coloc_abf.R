# Load required packages
#library(optparse)
library(data.table)
library(dplyr)
library(coloc)

# Function to run coloc
run_coloc <- function(
    protein_file,
    outcome_file,
    protein_name = "PROTEIN",
    outcome_name = "OUTCOME",
    outcome_N,
    outcome_s,
    protein_N = NA,
    window_bp = 250000,
    out_dir = ".",
    save = TRUE   # if TRUE writes summary CSV of full coloc object
) {
  # minimal deps
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table")
  if (!requireNamespace("coloc", quietly = TRUE)) stop("Install coloc")
  library(data.table)
  library(coloc)
  
  # basic file existence
  if (!file.exists(protein_file)) stop("protein_file not found: ", protein_file)
  if (!file.exists(outcome_file)) stop("outcome_file not found: ", outcome_file)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # read
  prot <- data.table::fread(protein_file)
  outcome <- data.table::fread(outcome_file)
  
  # Minimal column expectations (user-responsibility)
  req_prot_cols <- c("snp","chr","position","p","beta","SE")
  req_out_cols  <- c("snp","chr","position","beta","SE")
  
  miss_prot <- setdiff(req_prot_cols, names(prot))
  miss_out  <- setdiff(req_out_cols, names(outcome))
  
  if (length(miss_prot) > 0) {
    stop("Protein file missing required columns: ", paste(miss_prot, collapse = ", "))
  }
  if (length(miss_out) > 0) {
    stop("Outcome file missing required columns: ", paste(miss_out, collapse = ", "))
  }
  
  # compute varbeta
  prot[, varbeta := SE^2]
  outcome[, varbeta := SE^2]
  
  # find lead SNP in protein by lowest p
  lead_idx <- which.min(prot$p)
  lead_snp <- prot$snp[lead_idx]
  lead_chr <- as.character(prot$chr[lead_idx])
  lead_pos <- as.numeric(prot$position[lead_idx])
  
  if (is.na(lead_chr) || is.na(lead_pos)) stop("Lead SNP missing chr/position; cannot define window")
  
  window_min <- lead_pos - window_bp
  window_max <- lead_pos + window_bp
  
  message(sprintf("Lead SNP: %s (chr%s:%d). Window: %d-%d", lead_snp, lead_chr, lead_pos, window_min, window_max))
  
  # subset
  prot_sub <- prot[chr == lead_chr & position >= window_min & position <= window_max]
  out_sub  <- outcome[chr == lead_chr & position >= window_min & position <= window_max]
  
  if (nrow(prot_sub) == 0) stop("No protein variants in window")
  if (nrow(out_sub) == 0) stop("No outcome variants in window")
  
  # match by rsid
  common_snps <- intersect(prot_sub$snp, out_sub$snp)
  if (length(common_snps) == 0) stop("No overlapping rsIDs in window. Ensure rsIDs match between files.")
  
  # reorder to same SNP order
  prot_sub2 <- prot_sub[match(common_snps, prot_sub$snp)]
  out_sub2  <- out_sub[ match(common_snps, out_sub$snp)]
  
  # prepare datasets for coloc
  ds1 <- list(
    snp = out_sub2$snp,
    beta = as.numeric(out_sub2$beta),
    varbeta = as.numeric(out_sub2$varbeta),
    MAF = if ("MAF" %in% names(out_sub2)) as.numeric(out_sub2$MAF) else rep(NA_real_, nrow(out_sub2)),
    type = "cc",
    N = as.numeric(outcome_N),
    s = as.numeric(outcome_s)
  )
  
  ds2 <- list(
    snp = prot_sub2$snp,
    beta = as.numeric(prot_sub2$beta),
    varbeta = as.numeric(prot_sub2$varbeta),
    MAF = if ("MAF" %in% names(prot_sub2)) as.numeric(prot_sub2$MAF) else rep(NA_real_, nrow(prot_sub2)),
    type = "quant",
    N = if (!is.na(protein_N)) as.numeric(protein_N) else NA
  )
  
  message("Running coloc.abf ...")
  coloc_res <- coloc.abf(dataset1 = ds1, dataset2 = ds2)
  
  summary_df <- as.data.frame(coloc_res$summary)
  # add metadata
  summary_df$protein <- protein_name
  summary_df$outcome <- outcome_name
  summary_df$chr <- lead_chr
  summary_df$lead_snp <- lead_snp
  summary_df$lead_pos <- lead_pos
  summary_df$window_min <- window_min
  summary_df$window_max <- window_max
  summary_df$overlap_snps <- length(common_snps)
  
  if (save) {
    out_summary <- file.path(out_dir, paste0("coloc_summary_", protein_name, "_vs_", outcome_name, ".csv"))
    data.table::fwrite(summary_df, out_summary)
    message("Saved summary to: ", out_summary)
  }
  
  # return a small list
  return(list(
    summary = summary_df,
    coloc = coloc_res,
    prot_window = prot_sub2,
    outcome_window = out_sub2
  ))
}

 res <- run_coloc(
   protein_file = "path/to/protein/file",     # Columns expected: snp, chr, position, MAF, p, beta, se
   outcome_file = "path/to/outcome/file",     # Columns expected: snp, chr, position, MAF, p, beta, se
   protein_name = "PROTEIN",
   outcome_name = "OUTCOME",
   outcome_N = 1764486,                       # Total sample size, modify to use case
   outcome_s = 0.224,                         # Proportion of sample that are cases, modify to use case
   protein_N = 34557,
   out_dir = "/path/to/your/results/deposit/"
 )
