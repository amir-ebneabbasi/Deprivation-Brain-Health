# Load Libraries
library(GENESIS)
library(methods)
library(GWASTools)
library(SNPRelate)
library(BiocParallel)
library(SeqVarTools)

# Paths
# OLD
output_dir <- "path/to/output"
gdsfile <- file.path(output_dir, "genotype_autosomes.gds")
pruned_file <- file.path(output_dir, "pruned_snps.txt")
pc_rds_file <- file.path(output_dir, "processed_pc_matrix.rds")
unrelated_samples_file <- file.path(output_dir, "processed_unrelated_samples.rds")

# NEW
pcrelate_file <- file.path(output_dir, "pcrelate.rds")
pcrelate_sparse_file <- file.path(output_dir, "pcrelate_sparse.rds")

# Load PCs
cat("Loading PC matrix from RDS file...\n")
pc <- readRDS(pc_rds_file)
cat("PC matrix loaded from:", pc_rds_file, "\n")
cat("PC matrix dimensions:", dim(pc), "\n")
cat("PC matrix has", nrow(pc), "samples ×", ncol(pc), "PCs\n")

# Load unrelated samples
cat("Loading unrelated samples...\n")
unrel <- readRDS(unrelated_samples_file)
cat("Loaded", length(unrel), "unrelated samples\n")

# Load pruned SNPs
cat("Loading pruned SNP list...\n")
pruned <- scan(pruned_file, what = "")
cat("Loaded", length(pruned), "pruned SNPs\n")

# Load GDS genotype data
cat("Loading GDS file...\n")
gds_reader <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(gds_reader)
geno_ids <- getScanID(genoData)
cat("Loaded", length(geno_ids), "sample IDs from GDS\n")

# Run PC-Relate
cat("\nRunning PC-Relate...\n")
genoBlock <- GenotypeBlockIterator(genoData, snpInclude = pruned)

pcrelate <- pcrelate(
  genoBlock,
  training.set= unrel,
  pcs = pc,
  sample.block.size = 10000,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4)
)

saveRDS(pcrelate, pcrelate_file)
cat("PC-Relate saved to:", pcrelate_file, "\n")

# Create Sparse GRM
cat("\nCreating sparse GRM...\n")
pcrelate_sparse <- pcrelateToMatrix(
  pcrelate,
  thresh = 2^(-11/2)
)

# Save Results
saveRDS(pcrelate_sparse, pcrelate_sparse_file)
cat("Sparse GRM saved to:", pcrelate_sparse_file, "\n")

# Close GDS connection
close(genoData)
cat("GDS connection closed\n")
