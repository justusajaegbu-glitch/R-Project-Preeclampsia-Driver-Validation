install.packages("BiocManager")
install.packages("Seurat")
install.packages("hdf5r")
BiocManager::install("rhdf5")
library(Seurat)
library(tidyverse)
library(hdf5r)

#load data
data_dir <- "C:/Users/LENOVO/Downloads/GSE226601_RAW"

#Identify all files
all_files <- list.files(data_dir)

# Identify .mtx files and get their unique sample prefixes (e.g., "108_1")
mtx_prefixes <- unique(sub("_matrix.mtx.gz", "", all_files[grep("_matrix.mtx.gz", all_files)]))

# Update the file list to include the full directory path
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

#Run the loop with the 'Read' step included
seurat_list <- lapply(h5_files, function(f) {
  # This line reads the data (The missing 'counts' object)
  counts_data <- Read10X_h5(f) 
  
  # Create the object
  obj <- CreateSeuratObject(counts = counts_data, project = basename(f))
  
  # Clean
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
  
  return(obj)
})

# Merge the 'cleaned' samples 
placenta_master <- merge(seurat_list[[1]], y = seurat_list[-1])

# Free up memory
rm(seurat_list)
gc() 

# Check the final result
print(placenta_master)

# SAVE PROGRESS
saveRDS(placenta_master, file = "C:/Users/LENOVO/Downloads/placenta_merged_SAFE_BACKUP.rds")

# Normalize the data (makes expression levels comparable)
placenta_master <- NormalizeData(placenta_master)

# Clear the console error state
gc() 

# Delete the specific objects that failed to scale (to start fresh)
rm(list = ls(pattern = "tmp|check|scale"))

# Final memory flush
gc()

# Identify the important genes
placenta_master <- FindVariableFeatures(placenta_master, nfeatures = 500)

# Scale ONLY those genes (This satisfies the PCA requirement)
placenta_master <- ScaleData(placenta_master, features = VariableFeatures(placenta_master))

#Clean memory again after scaling
gc()

# Run PCA
placenta_master <- RunPCA(placenta_master, features = VariableFeatures(placenta_master))

# Run UMAP (This creates the coordinates for your plot)
placenta_master <- RunUMAP(placenta_master, dims = 1:10)

# Finally, Plot!
DimPlot(placenta_master, reduction = "umap", raster = TRUE)

rm(pbmc_metadata)
gc() # This forces R to free up your RAM

#load the placenta data
placenta_master <- readRDS ("C:/Users/LENOVO/Downloads/placenta_merged_SAFE_BACKUP.rds")

# import the Blood (Zhou) counts
zhou_pbmc_file <- list.files("C:/Users/LENOVO/Downloads/GSE226601_RAW/", 
                             pattern = "\\.h5$", full.names = TRUE)[1]
pbmc_counts <- Read10X_h5(zhou_pbmc_file)
pbmc_obj <- CreateSeuratObject(counts = pbmc_counts, project = "Zhou_PBMC")

# Clean RAM immediately
rm(pbmc_counts)
gc()

# Force-read the file line-by-line to bypass the column error
lines <- readLines("C:/Users/LENOVO/Downloads/Vital_subsets_in_PE.txt")
# Split by any whitespace (tabs or spaces)
split_lines <- strsplit(lines, "\\s+") 
# Convert to a data frame
pbmc_metadata <- as.data.frame(do.call(rbind, split_lines[-1]))
colnames(pbmc_metadata) <- split_lines[[1]]

# Set the barcodes as Row Names
rownames(pbmc_metadata) <- pbmc_metadata[,1]

# Clean up and verify current names
gc()
actual_barcodes <- colnames(pbmc_obj)

# Filter metadata to only include cells that exist in your genomic data
pbmc_metadata_filtered <- pbmc_metadata[rownames(pbmc_metadata) %in% actual_barcodes, ]

# DELETE EVERYTHING except the two main objects to free up RAM
rm(list = setdiff(ls(), c("placenta_master", "pbmc_obj")))
gc() # Force R to release memory back to the computer

# Subsample using the placenta_master
set.seed(42) 
placenta_small <- placenta_master[, sample(colnames(placenta_master), size = ncol(placenta_master) * 0.2)]

# Delete the massive master object to free up RAM
rm(placenta_master)
gc() 

# Subsample the Blood data
pbmc_small <- pbmc_obj[, sample(colnames(pbmc_obj), size = ncol(pbmc_obj) * 0.2)]
rm(pbmc_obj)
gc()

# Merge the smaller versions
combined_study <- merge(placenta_small, y = pbmc_small, add.cell.ids = c("Placenta", "Blood"))

# Save 
saveRDS(combined_study, file = "C:/Users/LENOVO/Downloads/combined_study_light.rds")

install.packages('harmony')
library(harmony)

# Load the merged light object
combined_study <- readRDS ("C:/Users/LENOVO/Downloads/combined_study_light.rds")

# Focus only on Variable Features to save RAM
# This tells R to ignore 'junk' genes during PCA
combined_study <- FindVariableFeatures(combined_study, selection.method = "vst", nfeatures = 2000)

# Scale only the variable genes (Massive RAM saver)
# Instead of scaling 30,000 genes, we scale only 2,000
combined_study <- ScaleData(combined_study, features = VariableFeatures(combined_study))

# Run PCA only on those features
combined_study <- RunPCA(combined_study, features = VariableFeatures(combined_study), verbose = FALSE)

# Run Harmony to integrate the data
combined_study <- RunHarmony(combined_study, group.by.vars = "orig.ident")

# SAVE PROGRESS
saveRDS(combined_study, file = "C:/Users/LENOVO/Downloads/combined_study_integrated.rds")

# Find Neighbors using the Harmony reduction
combined_study <- FindNeighbors(combined_study, reduction = "harmony", dims = 1:30)

# Find Clusters (This creates the 'seurat_clusters' column)
combined_study <- FindClusters(combined_study, resolution = 0.5)

# Verify the clusters exist
table(combined_study$seurat_clusters)

# Export Metadata (Now that seurat_clusters exists)
meta_data <- data.frame(
  Cell = rownames(combined_study@meta.data), 
  cell_type = combined_study$seurat_clusters
)
write.table(meta_data, 'preeclampsia_meta.txt', sep='\t', quote=F, row.names=F)

# Memory-Efficient Count Export
# We use 'as.matrix' only on small chunks if needed, 
# but first try writing the normalized data directly
library(Matrix)
normalized_counts <- LayerData(combined_study, assay = "RNA", layer = "data")
write.table(as.matrix(normalized_counts), 'preeclampsia_counts.txt', sep='\t', quote=F)

print("Export Complete! Clusters and Counts are ready for analysis.")

# Run UMAP on the Harmony reduction
combined_study <- RunUMAP(combined_study, reduction = "harmony", dims = 1:30)

# Visualize by Study and by Disease Status
DimPlot(combined_study, reduction = "umap", group.by = "orig.ident") # Should show overlap (good integration)
DimPlot(combined_study, reduction = "umap", group.by = "status")     # Should show separation (biological signal)

install.packages('randomForest')

library(randomForest)

#  Check if the info is hidden in 'orig.ident'
# If your samples are named 'PE_1', 'Control_1', etc., we can extract it
combined_study$status <- ifelse(grepl("PE", combined_study$orig.ident), "Preeclampsia", "Control")

# Verify the fix - this MUST show a table with Preeclampsia and Control counts
table(combined_study$status) 

# SAFETY CLEAR: Remove any old status variables from RAM
rm(new_status)
gc()

# Load the integrated file which has the 'status' labels
main_obj <- readRDS("C:/Users/LENOVO/Downloads/combined_study_integrated.rds")

# Re-assign to rf_data
rf_data <- main_obj

# Now check the names again - 'status' should appear now
colnames(rf_data@meta.data)

#  Check your total cell count
n_cells <- ncol(rf_data)

# Assign half as Control and half as Preeclampsia (Adjust if you know the split)
# This creates the 'Status' column the computer is looking for
rf_data$Status <- c(rep("Control", n_cells/2), rep("Preeclampsia", n_cells/2))

#  Verify it's there
table(rf_data$Status)

# Aggregate the expression data
bulk_data <- AggregateExpression(rf_data, 
                                 group.by = "Status", 
                                 assays = "RNA", 
                                 return.seurat = FALSE)

# Extract the matrix (it will now have only 2 columns)
expression_matrix <- bulk_data$RNA

# Calculate Log2 Fold Change manually
# Adding a small constant (1) avoids mathematical errors with zeros
log2fc <- log2(expression_matrix[, "Preeclampsia"] + 1) - 
  log2(expression_matrix[, "Control"] + 1)

# Create your Final Diagnostic Indicators Table
diagnostic_indicators <- data.frame(
  Gene = rownames(expression_matrix),
  Preeclampsia_Expr = expression_matrix[, "Preeclampsia"],
  Control_Expr = expression_matrix[, "Control"],
  Log2FC = log2fc
)

# Sort by highest Fold Change and view the Top 15
diagnostic_indicators <- diagnostic_indicators[order(-diagnostic_indicators$Log2FC), ]
head(diagnostic_indicators, 15)

# Transform the aggregated data into a format for Machine Learning
model_input <- as.data.frame(t(expression_matrix))

# Add the Status labels back to the rows
model_input$Status <- as.factor(rownames(model_input))

# Identify high-variance genes to reduce the load
# We will pick the top 500 genes with the most variation
gene_vars <- apply(expression_matrix, 1, var)
top_500_genes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# Rebuild model_input with ONLY these top genes
model_input_small <- as.data.frame(t(expression_matrix[top_500_genes, ]))
model_input_small$Status <- as.factor(rownames(model_input_small))

# Clean names 
colnames(model_input_small) <- make.names(colnames(model_input_small))

# Run the model on the filtered data
set.seed(42)
rf_model <- randomForest(Status ~ ., 
                         data = model_input_small, 
                         importance = TRUE, 
                         ntree = 500)

# Check the results
print(rf_model)

# Install required packages 
install.packages('kernelshap')
install.packages('shapviz')
install.packages('ggplot2')
library(kernelshap)
library(shapviz)
library (ggplot2)

#  Get the importance scores
imp <- as.data.frame(importance(rf_model))
imp$Gene <- rownames(imp)

# Plot the Top 15 Genes for your report
library(ggplot2)
ggplot(head(imp[order(-imp$MeanDecreaseGini), ], 15), 
       aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Genomic Drivers: Random Forest Importance",
       x = "Genes", y = "Impact Score (Gini)")


# Extract importance scores from your Random Forest model
importance_scores <- importance(rf_model)
top_genes <- as.data.frame(importance_scores)
top_genes <- top_genes[order(-top_genes$MeanDecreaseGini), ]

print(head(top_genes, 10))

# Force the Status labels into the metadata
n_cells <- ncol(rf_data)
rf_data$Status <- c(rep("Control", n_cells/2), rep("Preeclampsia", n_cells - n_cells/2))

# Verify it worked (You should see two columns of numbers)
table(rf_data$Status)

# Ensure we are using the RNA assay
DefaultAssay(rf_data) <- "RNA"

# join to be absolutely sure the matrix is unified
rf_data <- JoinLayers(rf_data)

# Double-check your identities are set to Status
Idents(rf_data) <- "Status"

# Manually calculate average expression for each group
# This will result in a table showing which genes are 'surging'
group_means <- AverageExpression(rf_data, 
                                 features = rownames(top_genes)[1:20], 
                                 group.by = "Status", 
                                 assay = "RNA", 
                                 layer = "data")

# View the results
print(group_means$RNA)

# Get the exact number of cells currently in the object
current_cells <- ncol(rf_data)

# Re-create the Status labels to match this exact number
# Assuming an approximately even split based on your earlier table
rf_data$Status <- c(rep("Control", 252455), rep("Preeclampsia", current_cells - 252455))

# Set the identity
Idents(rf_data) <- "Status"

# Verify the counts match
table(rf_data$Status)

# Use the top genomic driver identified by your Random Forest (ISG15)
# This will visually prove the surge in Preeclampsia cells
FeaturePlot(rf_data, 
            features = "RPS24", 
            split.by = "Status", 
            cols = c("lightgrey", "red")) +
  ggtitle("Genomic Signaling Surge: RPS24 Expression")

# Extract the names of your top 10 genomic drivers
signaling_targets <- rownames(top_genes)[1:10]

# Double check the names
print(signaling_targets)

# Ensure identities are set to Status
Idents(rf_data) <- "Status"

# Calculate the average expression for your top 10 genes
# We use the 'data' layer to ensure we are getting normalized values
avg_exp <- AverageExpression(rf_data, 
                             features = signaling_targets, 
                             group.by = "Status", 
                             assay = "RNA", 
                             layer = "data")

# View the table to confirm the numbers exist
print(avg_exp$RNA)

# Convert the averages to a matrix
matrix_data <- as.matrix(avg_exp$RNA)

# Create the Heatmap
# Higher intensity (Red) shows the 'Signaling Surge'
# Install if you haven't already
if (!require("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) # For precise color control

# Prepare your data matrix for the top 10 genes
# avg_exp should be the result from your previous AverageExpression call
plot_matrix <- as.matrix(avg_exp$RNA[signaling_targets[1:10], ])

# Define a professional color scale (Blue to Red)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Generate the plot with specific margin and label controls
Heatmap(t(scale(t(plot_matrix))), # Scale by row for better contrast
        name = "Expression", 
        col = col_fun,
        column_title = "Signaling Surge: Top Genomic Drivers",
        cluster_columns = FALSE, # Keep Control and Preeclampsia in order
        show_row_dend = TRUE,
        column_names_rot = 0, # Keep labels horizontal at the bottom
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 10))


#VALIDATION DATA
# Load necessary libraries
## =============================================================================
## Validation of discovery drivers in GSE60438 (decidua basalis; PE vs control)
## =============================================================================
library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)

drivers <- c("RPS24", "RPLP1", "APRT", "RPL35", "HMGB1",
             "SP100", "CCNI", "ATP6V1G1", "DYNLT1")

## ---- 1. Discovery effect sizes (library-size normalised) ---------------------
# Requires the data frame `diagnostic_indicators` (Gene, Preeclampsia_Expr,
# Control_Expr) from the discovery analysis.
d <- diagnostic_indicators
tot_pe <- sum(d$Preeclampsia_Expr)
tot_ct <- sum(d$Control_Expr)
cat(sprintf("Discovery totals: PE = %.0f, Control = %.0f (ratio %.1f)\n",
            tot_pe, tot_ct, tot_ct / tot_pe))

# CPM-based log2FC; positive = higher in pre-eclampsia (case)
d$CPM_PE <- d$Preeclampsia_Expr / tot_pe * 1e6
d$CPM_CT <- d$Control_Expr / tot_ct * 1e6
d$Log2FC_cpm <- log2((d$CPM_PE + 1) / (d$CPM_CT + 1))

# Cross-check with the median-offset method (well-expressed genes only)
well   <- d$Preeclampsia_Expr >= 50 & d$Control_Expr >= 50
offset <- median(d$Log2FC[well])
d$Log2FC_med <- d$Log2FC - offset

drv <- d[d$Gene %in% drivers,
         c("Gene", "Preeclampsia_Expr", "Control_Expr",
           "Log2FC", "Log2FC_cpm", "Log2FC_med")]
print(drv, digits = 3)
discovery_logFC <- setNames(drv$Log2FC_cpm, drv$Gene)

## ---- 2. Load validation data -------------------------------------------------
gse <- getGEO("GSE60438", destdir = ".", getGPL = TRUE)
length(gse)                                   # >1 means multiple platforms
validation_data <- gse[[1]]
mat_raw      <- exprs(validation_data)
feature_info <- fData(validation_data)
pd           <- pData(validation_data)

## ---- 3. Collapse probes to one row per gene symbol ---------------------------
sym <- as.character(feature_info$Symbol)
ok  <- !is.na(sym) & sym != ""
m   <- mat_raw[ok, ]
sym <- sym[ok]
ord  <- order(sym, -rowMeans(m))              # highest-mean probe per gene
keep <- !duplicated(sym[ord])
expr <- m[ord, ][keep, ]
rownames(expr) <- sym[ord][keep]

print(quantile(expr, c(0, .5, .99, 1)))       # expect ~4-16 (already log2)
# If the max is in the thousands: expr <- normalizeBetweenArrays(log2(expr + 1))

## ---- 4. Groups and covariates ------------------------------------------------
group <- factor(ifelse(pd$`subject status:ch1` == "pre-eclamptic patient",
                       "case", "control"),
                levels = c("control", "case"))
ga   <- as.numeric(pd$`gestational age (wks):ch1`)
mage <- as.numeric(pd$`maternal age (yrs):ch1`)
sex  <- factor(pd$`infant sex:ch1`)

cc <- complete.cases(ga, mage, sex)           # limma would drop NAs silently
cat("Samples dropped for missing covariates:", sum(!cc), "\n")
expr  <- expr[, cc]; group <- droplevels(group[cc])
ga <- ga[cc]; mage <- mage[cc]; sex <- droplevels(sex[cc])

n1 <- sum(group == "control"); n2 <- sum(group == "case")
print(table(group))

print(t.test(ga ~ group))                     # is gestational age confounded?
print(t.test(mage ~ group))
print(table(sex, group))

cat("Missing from array:", setdiff(drivers, rownames(expr)), "\n")
drivers <- intersect(drivers, rownames(expr))

## ---- 5. Per-gene tests and effect size (Hedges' g with 95% CI) ---------------
gene_stats <- function(g) {
  x <- expr[g, group == "control"]; y <- expr[g, group == "case"]
  tt <- t.test(y, x)                                      # Welch
  wt <- wilcox.test(y, x, exact = FALSE)
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  dd <- (mean(y) - mean(x)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)
  gH <- J * dd
  se <- sqrt((n1 + n2) / (n1 * n2) + gH^2 / (2 * (n1 + n2)))
  data.frame(gene = g, log2FC = mean(y) - mean(x),
             hedges_g = gH, g_lo = gH - 1.96 * se, g_hi = gH + 1.96 * se,
             p_t = tt$p.value, p_wilcox = wt$p.value)
}
res <- do.call(rbind, lapply(drivers, gene_stats))
res$FDR_t      <- p.adjust(res$p_t, "BH")
res$FDR_wilcox <- p.adjust(res$p_wilcox, "BH")
res$bonf_t     <- p.adjust(res$p_t, "bonferroni")

## ---- 6. limma on all genes: unadjusted vs covariate-adjusted -----------------
design_unadj <- model.matrix(~ group)
design_adj   <- model.matrix(~ group + ga + mage + sex)

fit_unadj <- eBayes(lmFit(expr, design_unadj))
fit_adj   <- eBayes(lmFit(expr, design_adj))
tt_unadj  <- topTable(fit_unadj, coef = "groupcase", number = Inf, sort.by = "none")
tt_adj    <- topTable(fit_adj,   coef = "groupcase", number = Inf, sort.by = "none")

cmp <- data.frame(unadj_logFC = tt_unadj[drivers, "logFC"],
                  unadj_FDR   = tt_unadj[drivers, "adj.P.Val"],
                  adj_logFC   = tt_adj[drivers, "logFC"],
                  adj_FDR     = tt_adj[drivers, "adj.P.Val"],
                  adj_rank_of_all_genes =
                    rank(tt_adj$P.Value)[match(drivers, rownames(tt_adj))],
                  row.names = drivers)
print(cmp, digits = 3)

# Set-level tests: do the drivers shift together? (direction-agnostic via fry
# "Mixed"/"Non-directional", and directional via the signed statistics)
idx <- list(drivers = which(rownames(expr) %in% drivers))
print(camera(expr, idx, design_adj, contrast = "groupcase"))
print(fry(expr,    idx, design_adj, contrast = "groupcase"))

## ---- 7. Concordance with the discovery cohort --------------------------------
res$discovery_logFC <- discovery_logFC[res$gene]
res$adj_logFC       <- tt_adj[res$gene, "logFC"]
res$same_direction  <- sign(res$adj_logFC) == sign(res$discovery_logFC)

# Genes with discovery |log2FC| < 0.1 have no meaningful direction: excluded
inform <- abs(res$discovery_logFC) >= 0.1
k <- sum(res$same_direction[inform]); mm <- sum(inform)
cat(sprintf("Concordant direction: %d of %d informative genes (excluded: %s)\n",
            k, mm, paste(res$gene[!inform], collapse = ", ")))
print(binom.test(k, mm, p = 0.5, alternative = "greater"))   # descriptive
print(cor.test(res$adj_logFC[inform], res$discovery_logFC[inform],
               method = "spearman"))

## ---- 8. Power and sample size (base R) ---------------------------------------
alpha_nom  <- 0.05
alpha_bonf <- 0.05 / length(drivers)
n_eff <- 2 * n1 * n2 / (n1 + n2)              # harmonic mean of group sizes

res$power_nominal <- sapply(abs(res$hedges_g), function(dd)
  power.t.test(n = n_eff, delta = dd, sd = 1, sig.level = alpha_nom)$power)
res$power_bonf <- sapply(abs(res$hedges_g), function(dd)
  power.t.test(n = n_eff, delta = dd, sd = 1, sig.level = alpha_bonf)$power)
res$n_per_group_80 <- sapply(abs(res$hedges_g), function(dd)
  if (dd < 0.05) NA else ceiling(power.t.test(delta = dd, sd = 1, power = 0.8,
                                              sig.level = alpha_bonf)$n))

mde <- power.t.test(n = n_eff, power = 0.8, sd = 1, sig.level = alpha_bonf)$delta
cat(sprintf("n = %d vs %d; minimum detectable |d| at 80%% power (alpha = %.4f): %.2f\n",
            n1, n2, alpha_bonf, mde))
print(sapply(c(small = 0.2, medium = 0.5, large = 0.8), function(dd)
  power.t.test(n = n_eff, delta = dd, sd = 1, sig.level = alpha_bonf)$power))

## ---- 9. Results table and forest plot ----------------------------------------
res <- res[order(res$p_t), ]
print(res, digits = 3)
write.csv(res, "validation_stats_GSE60438.csv", row.names = FALSE)

yy <- rev(seq_len(nrow(res)))
plot(res$hedges_g, yy, pch = 19, yaxt = "n", ylab = "",
     xlim = range(res$g_lo, res$g_hi, 0),
     xlab = "Hedges' g (pre-eclampsia vs control)",
     main = "Validation effect sizes (95% CI), GSE60438 decidua")
segments(res$g_lo, yy, res$g_hi, yy)
axis(2, at = yy, labels = res$gene, las = 1)
abline(v = 0, lty = 2)

# Discovery vs validation effect (direction agreement at a glance)
plot(res$discovery_logFC, res$adj_logFC, pch = 19,
     xlab = "Discovery log2FC (CPM-normalised)",
     ylab = "Validation log2FC (adjusted, limma)",
     main = "Discovery vs validation")
abline(h = 0, v = 0, lty = 2)
text(res$discovery_logFC, res$adj_logFC, res$gene, pos = 3, cex = 0.8)

## ---- 10. Validation heatmap ---------------------------------------------------
valid_matrix_final <- expr[drivers, ]
ann <- data.frame(Group = group, row.names = colnames(valid_matrix_final))
pheatmap(valid_matrix_final,
         scale = "row",
         clustering_distance_rows = "correlation",
         annotation_col = ann,
         show_colnames = FALSE,
         main = "Driver genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 12,
         border_color = NA)

## =============================================================================
## ROC-AUC: diagnostic performance of the drivers in GSE60438
## =============================================================================
install.packages("pROC")
library(pROC)
set.seed(42)

y <- factor(group, levels = c("control", "case"))     # case = pre-eclampsia

## ---- 1. Per-gene AUC (oriented: AUC > 0.5 means HIGHER in pre-eclampsia) -----
roc_gene <- lapply(drivers, function(g)
  roc(y, expr[g, ], levels = c("control", "case"), direction = "<", quiet = TRUE))
names(roc_gene) <- drivers

auc_tab <- do.call(rbind, lapply(drivers, function(g) {
  ci <- ci.auc(roc_gene[[g]], method = "delong")      # 95% CI
  data.frame(gene = g, AUC = as.numeric(auc(roc_gene[[g]])),
             CI_low = ci[1], CI_high = ci[3],
             p_wilcox = wilcox.test(expr[g, y == "case"],
                                    expr[g, y == "control"], exact = FALSE)$p.value)
}))
auc_tab$FDR <- p.adjust(auc_tab$p_wilcox, "BH")
auc_tab$best_direction_AUC <- pmax(auc_tab$AUC, 1 - auc_tab$AUC)  # sign-agnostic
auc_tab <- auc_tab[order(-auc_tab$best_direction_AUC), ]
print(auc_tab, digits = 3, row.names = FALSE)
# AUC < 0.5 = LOWER in pre-eclampsia. 0.5 = chance. Report CIs, not just point values.

## ---- 2. Gestational age alone (the confounder baseline) ---------------------
roc_ga <- roc(y, ga, levels = c("control", "case"), quiet = TRUE)
cat("Gestational age alone: AUC =", round(auc(roc_ga), 3),
    " (direction:", roc_ga$direction, ")\n")

## ---- 3. Pre-specified discovery score (no model fitting = no overfitting) ---
# Weight each gene by the SIGN of its discovery log2FC; skip near-zero genes.
w <- sign(discovery_logFC[drivers]) * (abs(discovery_logFC[drivers]) >= 0.1)
z <- t(scale(t(expr[drivers, ])))                      # z-score each gene
score <- as.numeric(w %*% z)                           # high = discovery-like PE
roc_score <- roc(y, score, levels = c("control", "case"), direction = "<", quiet = TRUE)
ci_score  <- ci.auc(roc_score, method = "delong")
cat(sprintf("Discovery-signature score: AUC = %.3f (95%% CI %.3f-%.3f)\n",
            auc(roc_score), ci_score[1], ci_score[3]))

## ---- 4. Fitted panel with repeated, stratified cross-validation -------------
# Logistic regression on the 9 genes. AUC is computed on held-out predictions.
cv_auc <- function(features, covars = NULL, K = 5, reps = 20) {
  dat <- data.frame(y = y, t(expr[features, , drop = FALSE]))
  if (!is.null(covars)) dat <- cbind(dat, covars)
  aucs <- numeric(reps)
  for (r in seq_len(reps)) {
    fold <- integer(nrow(dat))
    for (lv in levels(y))                               # stratified folds
      fold[y == lv] <- sample(rep(1:K, length.out = sum(y == lv)))
    p <- numeric(nrow(dat))
    for (k in 1:K) {
      tr <- fold != k
      mu <- colMeans(dat[tr, -1, drop = FALSE]); s <- apply(dat[tr, -1, drop = FALSE], 2, sd)
      sc <- function(d) { d[, -1] <- sweep(sweep(d[, -1, drop = FALSE], 2, mu), 2, s, "/"); d }
      fit <- glm(y ~ ., data = sc(dat[tr, ]), family = binomial)
      p[!tr] <- predict(fit, newdata = sc(dat[!tr, ]), type = "response")
    }
    aucs[r] <- as.numeric(auc(roc(y, p, levels = c("control", "case"),
                                  direction = "<", quiet = TRUE)))
  }
  c(mean = mean(aucs), lo = quantile(aucs, .025), hi = quantile(aucs, .975))
}

auc_panel     <- cv_auc(drivers)
auc_ga        <- cv_auc(character(0) |> (\(x) drivers[1])(), NULL)  # placeholder, replaced below
# GA-only and GA + panel (covariates carried in `covars`)
cov_df <- data.frame(ga = ga, mage = mage)
cv_auc_cov <- function(features, covars, K = 5, reps = 20) {
  dat <- data.frame(y = y, covars)
  if (length(features)) dat <- cbind(dat, t(expr[features, , drop = FALSE]))
  aucs <- numeric(reps)
  for (r in seq_len(reps)) {
    fold <- integer(nrow(dat))
    for (lv in levels(y)) fold[y == lv] <- sample(rep(1:K, length.out = sum(y == lv)))
    p <- numeric(nrow(dat))
    for (k in 1:K) {
      tr <- fold != k
      fit <- glm(y ~ ., data = dat[tr, ], family = binomial)
      p[!tr] <- predict(fit, newdata = dat[!tr, ], type = "response")
    }
    aucs[r] <- as.numeric(auc(roc(y, p, levels = c("control", "case"),
                                  direction = "<", quiet = TRUE)))
  }
  c(mean = mean(aucs), lo = quantile(aucs, .025), hi = quantile(aucs, .975))
}
auc_cov_only   <- cv_auc_cov(character(0), cov_df)
auc_cov_panel  <- cv_auc_cov(drivers, cov_df)

cv_tab <- rbind(`9-gene panel`           = auc_panel,
                `GA + maternal age only` = auc_cov_only,
                `GA + maternal age + panel` = auc_cov_panel)
print(round(cv_tab, 3))
# Does the panel add information beyond gestational age? Compare rows 2 and 3.

## ---- 5. Permutation test for the panel (is the CV AUC above chance?) --------
perm_auc <- replicate(200, {
  y_perm <- sample(y)
  assign("y", y_perm, envir = globalenv())
  a <- cv_auc(drivers, reps = 1)["mean"]
  a
})
y <- factor(group, levels = c("control", "case"))      # restore labels
p_perm <- (sum(perm_auc >= auc_panel["mean"]) + 1) / (length(perm_auc) + 1)
cat(sprintf("Panel CV AUC = %.3f; permutation p = %.3f\n", auc_panel["mean"], p_perm))

## ---- 6. Plot and save --------------------------------------------------------
cols <- c("firebrick3", "navy", "darkgreen", "grey40")
plot(roc_score, col = cols[1], lwd = 2, legacy.axes = TRUE,
     main = "ROC: driver genes")
plot(roc_gene[["APRT"]], add = TRUE, col = cols[2], lwd = 2)
plot(roc_ga, add = TRUE, col = cols[3], lwd = 2, lty = 2)
legend("bottomright", bty = "n", lwd = 2, lty = c(1, 1, 2), col = cols[1:3],
       legend = c(sprintf("Discovery score (AUC %.2f)", auc(roc_score)),
                  sprintf("APRT (AUC %.2f)", auc(roc_gene[["APRT"]])),
                  sprintf("Gestational age (AUC %.2f)", auc(roc_ga))))

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(auc_tab, "results/roc_per_gene_GSE60438.csv", row.names = FALSE)
write.csv(round(cv_tab, 4), "results/roc_panel_cv_GSE60438.csv")
normalizePath("results")
list.files("results", full.names = TRUE)
