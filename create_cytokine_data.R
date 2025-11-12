# Process Cytokine ELISA/Multiplex Data for Interactive Visualization
# Calculate statistics and export to JSON

library(jsonlite)

# Read the CSV
cat("Reading cytokine data...\n")
data <- read.csv("C:/Users/racacer/OneDrive - Emory/claude_onedrive/01_Active_Projects/PDX59_RNAseq/_SECREATED_ELISA_MULTIPLEX_CCHMC_.csv",
                 skip = 2,  # Skip first 2 rows (column numbers and cell line)
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

cat("Data dimensions:", nrow(data), "x", ncol(data), "\n")

# First row is treatments
treatments_raw <- as.character(data[1, -1])  # Exclude first column (empty)
cat("Treatments:", paste(treatments_raw, collapse = ", "), "\n")

# Get cytokine data (rows 2 onward)
cytokine_data_raw <- data[-1, ]

# Extract cytokine names (first column)
cytokine_names <- cytokine_data_raw[, 1]
cat("Total cytokines:", length(cytokine_names), "\n")

# Extract values (columns 2 onward)
values_matrix <- as.matrix(cytokine_data_raw[, -1])

# Convert to numeric, replacing empty strings or non-numeric with NA
values_matrix <- apply(values_matrix, c(1, 2), function(x) {
  x <- as.character(x)
  if (x == "" || is.na(x)) return(NA)
  as.numeric(x)
})

# Group treatments
# Control: columns 1-3
# CAY10566 (SCDi): columns 4-6
# CAY10566_OA: columns 7-9
# CAY10566_PO: columns 10-12

control_idx <- 1:3
scdi_idx <- 4:6
scdi_oa_idx <- 7:9
scdi_po_idx <- 10:12

# Welch's t-test function
welch_t_test <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2 || length(y) < 2) return(1.0)
  if (sd(x) == 0 && sd(y) == 0) return(1.0)

  tryCatch({
    t_result <- t.test(x, y, var.equal = FALSE)
    return(t_result$p.value)
  }, error = function(e) {
    return(1.0)
  })
}

# Calculate SEM
calc_sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(0)
  return(sd(x) / sqrt(length(x)))
}

# Process each cytokine
cat("Processing cytokines...\n")
cytokine_list <- list()

for (i in 1:length(cytokine_names)) {
  cytokine <- cytokine_names[i]

  # Skip if cytokine name is empty or NA
  if (is.na(cytokine) || cytokine == "" || cytokine == "NA") next

  # Get values for this cytokine
  control_vals <- values_matrix[i, control_idx]
  scdi_vals <- values_matrix[i, scdi_idx]
  scdi_oa_vals <- values_matrix[i, scdi_oa_idx]
  scdi_po_vals <- values_matrix[i, scdi_po_idx]

  # Remove NAs
  control_vals_clean <- control_vals[!is.na(control_vals)]
  scdi_vals_clean <- scdi_vals[!is.na(scdi_vals)]
  scdi_oa_vals_clean <- scdi_oa_vals[!is.na(scdi_oa_vals)]
  scdi_po_vals_clean <- scdi_po_vals[!is.na(scdi_po_vals)]

  # Skip if control has no values
  if (length(control_vals_clean) == 0) next

  # Calculate p-values vs Control
  p_scdi <- welch_t_test(scdi_vals_clean, control_vals_clean)
  p_scdi_oa <- welch_t_test(scdi_oa_vals_clean, control_vals_clean)
  p_scdi_po <- welch_t_test(scdi_po_vals_clean, control_vals_clean)

  # Calculate normalized values (fold change relative to control)
  control_mean <- mean(control_vals_clean, na.rm = TRUE)

  scdi_vals_norm <- if (length(scdi_vals_clean) > 0) scdi_vals_clean / control_mean else numeric(0)
  scdi_oa_vals_norm <- if (length(scdi_oa_vals_clean) > 0) scdi_oa_vals_clean / control_mean else numeric(0)
  scdi_po_vals_norm <- if (length(scdi_po_vals_clean) > 0) scdi_po_vals_clean / control_mean else numeric(0)
  control_vals_norm <- control_vals_clean / control_mean

  cytokine_list[[cytokine]] <- list(
    Control = list(
      mean_raw = mean(control_vals_clean, na.rm = TRUE),
      sem_raw = calc_sem(control_vals_clean),
      values_raw = as.list(control_vals_clean),
      mean_norm = mean(control_vals_norm, na.rm = TRUE),
      sem_norm = calc_sem(control_vals_norm),
      values_norm = as.list(control_vals_norm),
      pValue = 1.0
    ),
    SCDi = list(
      mean_raw = if (length(scdi_vals_clean) > 0) mean(scdi_vals_clean, na.rm = TRUE) else NA,
      sem_raw = calc_sem(scdi_vals_clean),
      values_raw = as.list(scdi_vals_clean),
      mean_norm = if (length(scdi_vals_norm) > 0) mean(scdi_vals_norm, na.rm = TRUE) else NA,
      sem_norm = calc_sem(scdi_vals_norm),
      values_norm = as.list(scdi_vals_norm),
      pValue = p_scdi
    ),
    `SCDi+OA` = list(
      mean_raw = if (length(scdi_oa_vals_clean) > 0) mean(scdi_oa_vals_clean, na.rm = TRUE) else NA,
      sem_raw = calc_sem(scdi_oa_vals_clean),
      values_raw = as.list(scdi_oa_vals_clean),
      mean_norm = if (length(scdi_oa_vals_norm) > 0) mean(scdi_oa_vals_norm, na.rm = TRUE) else NA,
      sem_norm = calc_sem(scdi_oa_vals_norm),
      values_norm = as.list(scdi_oa_vals_norm),
      pValue = p_scdi_oa
    ),
    `SCDi+PO` = list(
      mean_raw = if (length(scdi_po_vals_clean) > 0) mean(scdi_po_vals_clean, na.rm = TRUE) else NA,
      sem_raw = calc_sem(scdi_po_vals_clean),
      values_raw = as.list(scdi_po_vals_clean),
      mean_norm = if (length(scdi_po_vals_norm) > 0) mean(scdi_po_vals_norm, na.rm = TRUE) else NA,
      sem_norm = calc_sem(scdi_po_vals_norm),
      values_norm = as.list(scdi_po_vals_norm),
      pValue = p_scdi_po
    )
  )
}

cat("Total cytokines in final dataset:", length(cytokine_list), "\n")

# Write JSON
output_file <- "C:/Users/racacer/OneDrive - Emory/claude_onedrive/01_Active_Projects/PDX59_RNAseq/work/for github/FINAL_SITE/cytokine_data_optimized.json"
cat("Writing JSON file to:", output_file, "\n")

json_output <- toJSON(cytokine_list, auto_unbox = TRUE, pretty = FALSE)
write(json_output, file = output_file)

cat("\n✅ Done! Created", output_file, "\n")
cat("   Cytokines included:", length(cytokine_list), "\n")

# Sample output
cat("\n=== Sample: IL-6 ===\n")
if ("IL-6 (pg/ML)" %in% names(cytokine_list)) {
  print(cytokine_list[["IL-6 (pg/ML)"]])
}
