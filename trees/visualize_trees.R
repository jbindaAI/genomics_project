suppressMessages(suppressWarnings(library(ape)))

tree_files <- c("trees/consensus_results/bacteria/vanilla/consensus_tree.contree",
                "trees/consensus_results/bacteria/boostrapped/consensus_tree.contree",
                "trees/super_tree_results/bacteria/vanilla/supertree.nwk",
                "trees/super_tree_results/bacteria/boostrapped/supertree.nwk")


trees <- lapply(tree_files, read.tree)

output_dir <- "trees/figures"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Output filenames for saving the plots
output_files <- file.path(output_dir, c("consensus_tree.pdf",
                                        "bootstraped_consensus_tree.pdf",
                                        "supertree.pdf",
                                        "bootstraped_supertree.pdf"))

output_png_files <- file.path(output_dir, c("consensus_tree.png",
                                            "bootstraped_consensus_tree.png",
                                            "supertree.png",
                                            "bootstraped_supertree.png"))

# Save the plots
for (i in seq_along(trees)) {
  pdf(output_files[i], width = 8, height = 6)
  plot(trees[[i]], main = gsub(".pdf", "", basename(output_files[i]), fixed = TRUE))
  dev.off()
  
  png(output_png_files[i], width = 800, height = 600)
  plot(trees[[i]], main = gsub(".png", "", basename(output_png_files[i]), fixed = TRUE))
  dev.off()
}

cat("Plots have been saved to the directory: ", output_dir, "\n")
