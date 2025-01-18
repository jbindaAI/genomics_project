# build_supertree.R
suppressMessages(suppressWarnings(library(phangorn)))
suppressMessages(suppressWarnings(library(ape)))

args <- commandArgs(trailingOnly = TRUE)

tree_file <- args[1]   # First argument: input tree file
output_file <- args[2] # Second argument: output file for the supertree
method <- args[3] 
multicore <- args[4]
n_cores <- args[5]

if (n_cores == 0) {
n_cores = NULL
}

trees <- read.tree(tree_file)

st_mrp <- superTree(trees, method = method, multicore=multicore, mc.cores=n_cores, minit = 25, maxit = 50)

write.tree(st_mrp, file = output_file)
