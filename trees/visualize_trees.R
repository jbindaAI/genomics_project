suppressMessages(suppressWarnings(library(ape)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R <basename>")
}
basename <- args[1]

tree_files <- c(
  sprintf("trees/consensus_results/%s/boostrapped/consensus_tree.contree", basename),
  sprintf("trees/consensus_results/%s/vanilla/consensus_tree.contree", basename),
  sprintf("trees/super_tree_results/ortologs/%s/boostrapped/supertree.nwk", basename),
  sprintf("trees/super_tree_results/ortologs/%s/vanilla/supertree.nwk", basename),
  sprintf("trees/super_tree_results/paralogs/%s/boostrapped/supertree.nwk", basename),
  sprintf("trees/super_tree_results/paralogs/%s/vanilla/supertree.nwk", basename)
)

existing_tree_files <- tree_files[file.exists(tree_files)]
if (length(existing_tree_files) == 0) {
  stop("No valid tree files found for the provided basename: ", basename)
}

# Validate and load the tree files
trees <- list()
valid_files <- list()

for (file in existing_tree_files) {
  tree <- tryCatch(
    read.tree(file),
    error = function(e) {
      warning(sprintf("Failed to read tree file: %s. Skipping.", file))
      return(NULL)
    }
  )

  # Check if the tree object is valid
  if (!is.null(tree) && !is.null(tree$tip.label)) {
    trees <- append(trees, list(tree))
    valid_files <- append(valid_files, file)
  } else {
    warning(sprintf("Tree file %s is invalid or empty. Skipping.", file))
  }
}

if (length(trees) == 0) {
  stop("No valid trees could be loaded. Exiting.")
}


output_dir <- file.path("trees", "figures")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_files <- c(
  file.path(output_dir, "consensus_tree_boostrapped.pdf"),
  file.path(output_dir, "consensus_tree_vanilla.pdf"),
  file.path(output_dir, "supertree_ortologs_boostrapped.pdf"),
  file.path(output_dir, "supertree_ortologs_vanilla.pdf"),
  file.path(output_dir, "supertree_paralogs_boostrapped.pdf"),
  file.path(output_dir, "supertree_paralogs_vanilla.pdf")
)

output_png_files <- c(
  file.path(output_dir, "consensus_tree_boostrapped.png"),
  file.path(output_dir, "consensus_tree_vanilla.png"),
  file.path(output_dir, "supertree_ortologs_boostrapped.png"),
  file.path(output_dir, "supertree_ortologs_vanilla.png"),
  file.path(output_dir, "supertree_paralogs_boostrapped.png"),
  file.path(output_dir, "supertree_paralogs_vanilla.png")
)

# Ensure the number of output files matches the number of valid input files
if (length(trees) != length(output_files)) {
  stop("The number of valid trees does not match the hard-coded output file names.")
}

# Plot and save the trees
for (i in seq_along(trees)) {
  # Save to PDF
  pdf(output_files[i], width = 8, height = 6)
  plot(trees[[i]], main = gsub(".pdf", "", basename(output_files[i]), fixed = TRUE))
  dev.off()

  # Save to PNG
  png(output_png_files[i], width = 800, height = 600)
  plot(trees[[i]], main = gsub(".png", "", basename(output_png_files[i]), fixed = TRUE))
  dev.off()
}

cat("Plots have been saved to the directory: ", output_dir, "\n")
