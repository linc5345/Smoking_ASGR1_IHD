library(pacman)
# Packages needed for perform_mr() function
p_load(rio, tidyverse, data.table, MendelianRandomization,
       testthat, grid,gridExtra,ieugwasr,patchwork,TwoSampleMR,gsmr)


##########################################################################################################################################
# Ensure 'perform_mr()' function is defined - see 'https://github.com/adamkvonende/mr_pipeline' for details incl. column name requirements
##########################################################################################################################################

#####User inputs (please edit these paths as necessary)#####
exp_path        <- "/path/to/your/instrument/exp.csv" # Your exposure summary statistics (i.e. SmkInit, LSI, CigDay)
folder_proteins <- "/path/to/olink/proteins/"         # Folder containing Olink protein summary statistics
out_dir         <- "/path/to/your/results/deposit/"   # Output directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
final_results_path <- file.path(out_dir, "results.csv")
final_het_path     <- file.path(out_dir, "het.csv")

# Load in your instrument
exp <- fread(exp_path) # perform_mr() expects certain column names, view function

# Load in gene coordinates for proteins. Required columns: Assay, CHR, Gene_start, Gene_end
gene_coords <- fread("path/to/protein_list_with_coords.csv")

# Extract a short protein tag from filename (for ease of reading, modify to use case)
protein_tag <- function(fpath) {
  b <- basename(fpath)
  b <- str_remove(b, "\\.gz$")
  str_replace(b, "([A-Za-z0-9_-]+).*", "\\1")
}

append_write <- function(dt, out_path) {
  if (is.null(dt)) return(invisible())
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (nrow(dt) == 0) return(invisible())
  if (file.exists(out_path)) {
    fwrite(dt, out_path, append = TRUE)
  } else {
    fwrite(dt, out_path)
  }
}

# Define list of proteins to feed into loop (modify to use case)
protein_files <- list.files(folder_proteins, pattern = "\\.gz$", full.names = TRUE)


#####Main loop#####
for (f in protein_files) {
  pname <- protein_tag(f)
  message(sprintf(">> Processing: %s", pname))
  
  # Read outcome file
  outcome <- tryCatch({
    fread(f)
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", f, e$message))
    return(NULL)
  })
  if (is.null(outcome)) next
  
  # Exclude cis protein region
  # quick column normalisation for outcome to ensure chr/pos present for filtering
  if ("CHR" %in% names(outcome)) setnames(outcome, "CHR", "chr")
  if ("BP38" %in% names(outcome)) setnames(outcome, "BP38", "pos")
  
  # If gene_coords uses uppercase CHR field, normalise to 'chr'
  if ("CHR" %in% names(gene_coords)) setnames(gene_coords, "CHR", "chr")
  
  # derive gene symbol from protein name (leading token before underscore)
  gene_sym <- if (grepl("_", pname)) sub("_.*$", "", pname) else pname
  gene_sym <- trimws(as.character(gene_sym))
  
  # attempt simple matches: case-insensitive exact first, then substring
  gc_row <- gene_coords[toupper(Assay) == toupper(gene_sym)]
  if (nrow(gc_row) == 0) {
    gc_row <- gene_coords[grepl(gene_sym, Assay, ignore.case = TRUE)]
  }
  
  if (nrow(gc_row) == 0) {
    message(sprintf("  > No gene coords found for %s (gene_sym='%s'); skipping cis exclusion.", pname, gene_sym))
    outcome_filtered <- outcome   # no change
  } else {
    # take the first matching row (simple behaviour)
    gc <- gc_row[1]
    # ensure numeric positions present
    gc$Gene_start <- as.numeric(gc$Gene_start)
    gc$Gene_end   <- as.numeric(gc$Gene_end)
    gc$chr        <- as.integer(gc$chr)
    
    if (is.na(gc$Gene_start) || is.na(gc$Gene_end) || is.na(gc$chr)) {
      message(sprintf("  > Gene coords for %s contain NA; skipping cis exclusion.", pname))
      outcome_filtered <- outcome
    } else {
      start_exclude <- gc$Gene_start - 1e6
      end_exclude   <- gc$Gene_end   + 1e6
      if (start_exclude < 0) start_exclude <- 0L
      
      # ensure outcome chr/pos are integer-type for comparison
      if (!("chr" %in% names(outcome)) || !("pos" %in% names(outcome))) {
        warning(sprintf("  > outcome for %s lacks chr/pos columns; cannot apply cis exclusion. Proceeding without exclusion.", pname))
        outcome_filtered <- outcome
      } else {
        outcome[, chr := as.integer(chr)]
        outcome[, pos := as.integer(pos)]
        before_n <- nrow(outcome)
        outcome_filtered <- outcome[!(chr == gc$chr & pos >= start_exclude & pos <= end_exclude)]
        removed <- before_n - nrow(outcome_filtered)
        message(sprintf("  > Excluded %d variants in cis region chr%s:%d-%d (+/-1Mb) for %s",
                        removed, gc$chr, gc$Gene_start, gc$Gene_end, pname))
        if (nrow(outcome_filtered) == 0) {
          warning(sprintf("  > No variants remain after cis exclusion for %s — skipping this protein.", pname))
          next
        }
      }
    }
  }
  
  # Run MR
  res <- tryCatch({
    mr <- perform_mr(
      ss_exposure   = exp,
      ss_outcome    = outcome,
      cis           = FALSE,
      maf_thresh    = 0.01,
      p_thresh      = 5e-8,
      harmo_action  = 1,
      clump_when    = "after",
      clump_locally = FALSE,
      clump_kb      = 10000,
      clump_r2      = 0.001,
      clump_pop     = "EUR",
      name_exp      = "exp",
      name_outcome  = pname,
      which_methods = "basic",
      plots         = FALSE,
      outcome_binary= FALSE
    )
    list(results = mr$results, het = mr$het)
  }, error = function(e) {
    warning(sprintf("perform_mr failed for %s: %s", pname, e$message))
    return(NULL)
  })
  
  if (is.null(res)) next
  
  # Prepare and write results
  R <- as.data.table(res$results)
  if (!"outcome" %in% names(R)) R[, outcome := pname]
  if (all(c("b", "se") %in% names(R))) {
    R[, ll := b - 1.96*se]
    R[, ul := b + 1.96*se]
  }
  append_write(R, final_results_path)
  
  H <- as.data.table(res$het)
  if (!is.null(H) && nrow(H)) {
    if (!"outcome" %in% names(H)) H[, outcome := pname]
    append_write(H, final_het_path)
  }
  
  # Memory cleanup
  rm(outcome, res, R, H)
  gc()
}