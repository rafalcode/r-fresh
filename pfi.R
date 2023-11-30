#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(pathfindR)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()

# example_pathfindR_input
# available immediately?

# Gene.symbol	logFC	adj.P.Val
# FAM110A	-0.6939359	0.0000034
# RNASE2	1.3535040	0.0000101
# S100A8	1.5448338	0.0000347
# S100A9	1.0280904	0.0002263
# TEX261	-0.3235994	0.0002263
# ARHGAP17	-0.6919330	0.0002708

# By default, run_pathfindR() creates a temporary directory for writing the output files, including active subnetwork search results and a HTML report. To set the output directory, use output_dir:
# output_df <- run_pathfindR(example_pathfindR_input, p_val_threshold = 0.01, output_dir = "pfidir")
# output_df <- run_pathfindR(example_pathfindR_input, p_val_threshold = 0.01)

input_processed <- input_processing(example_pathfindR_input)
visualize_terms(result_df = example_pathfindR_output,
                input_processed = input_processed,
                hsa_KEGG=T)
