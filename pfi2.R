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

# follwinf function copied from source
create_HTML_report2 <- function(input, input_processed, final_res, dir_for_report)
{
    message("## Creating HTML report")
    rmarkdown::render(input="pfr_results.Rmd", output_dir = dir_for_report)
    rmarkdown::render(input="pfr_enriched_terms.Rmd",
        params = list(df = final_res), output_dir = dir_for_report)
    rmarkdown::render(input="pfr_conversion_table.Rmd",
        params = list(df = input_processed, original_df = input), output_dir = dir_for_report)
}

example_processed <- input_processing(example_pathfindR_input)
#
kegg_list <- fetch_gene_set(
  gene_sets = "KEGG",
  min_gset_size = 10,
  max_gset_size = 300)
# you get a list of two from this.
# [[1]] are the genesets
# and [2]] are the desciption

n_iter <- 2 ## number of iterations

run_as <- F
if(run_as) {
    combined_res <- NULL ## to store the result of each iteration
    for (i in 1:n_iter) {
      ###### Active Subnetwork Search
      snws_file <- paste0("active_snws_", i) # Name of output file
      active_snws <- active_snw_search(
        input_for_search = example_processed,
        pin_name_path = "Biogrid",
        snws_file = snws_file,
        score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
        sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
        search_method = "GR", # we suggest using GR
        seedForRandom=i) # setting seed to ensure reproducibility per iteration
    
      ###### Enrichment Analyses
      current_res <- enrichment_analyses(
        snws = active_snws,
        sig_genes_vec = example_processed$GENE,
        pin_name_path = "Biogrid",
        genes_by_term = kegg_list$genes_by_term,
        term_descriptions = kegg_list$term_descriptions,
        adj_method="bonferroni",
        enrichment_threshold=.05,
        list_active_snw_genes=TRUE) # listing the non-input active snw genes in output
    
      ###### Combine results via `rbind`
      combined_res <- rbind(combined_res, current_res)
    }
}

###### Summarize Combined Enrichment Results
summarized_df <- summarize_enrichment_results(enrichment_res=combined_res, list_active_snw_genes=T)

###### Annotate Affected Genes Involved in Each Enriched Term
final_res <- annotate_term_genes(
  result_df = summarized_df,
  input_processed = example_processed,
  genes_by_term = kegg_list$genes_by_term)

# quite naff as it turns out, only tables.
# create_HTML_report2(example_pathfindR_input, example_processed, final_res, "pfi2outdir")

CairoPNG("pfivt.png", 800, 800)
visualize_terms(
  result_df = final_res,
  hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
  pin_name_path = "Biogrid")
dev.off()

CairoPNG("pfiec.png", 800, 800)
enrichment_chart(final_res[1:10, ])
dev.off()
