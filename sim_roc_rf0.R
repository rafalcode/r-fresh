library(genefilter)
library(MASS)
library(pROC)
library(docopt)
library(DGCA)
library(Cairo)

# sigma is the pairwise correlation parameter matrix. AMcK follows MASS in this, gets ready to use MASS::mvrnorm
# once set up rnbinom  (or whatever) can be run on it
# to geenrate the variates.
# goc is most likely gain of correlation, loc is loss of (most likely)

# they to understandig this is the mvrnorm function and the Sigma option. This is the covraiance marix and says
# how one variable (gene) varies with respect to the variation of another.
# so it is not at all like rnorm whihc only has one variable. Sigma is what defines the
# correlation bewteen two variables, so it's no surprisie that alot of this code
# is about setting it properly.
args = commandArgs(trailingOnly = TRUE)
print(args)

nsamps_per_cond <- 6 # actually it will be nsamps per conditions ther being two conditions so N=40.
seed = runif(1, 0, 100000)
seed = as.numeric(42)

ddcor_test = T
EBcoexpress_test = F
discordant_test = F
compare_ddcor_eb_discordant = F

compare_ddcor_classes = T
compare_eb_classes = F
compare_discordant_classes = F

#reproducibility
set.seed(seed)

#color palette
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette = c("#000000", "darkgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")

#################################
# functions

# Sets the superdiagonal and subdiagonal of a matrix to the number x.
# so he calls a matrix variable by reserved word matrix - ooh great. Changed to "mat"
subsup_diag <- function(mat, x)
{
    # this is the key func that matches up g1 with g2, then g2 with g3 and so on.
    # the subdiag just reverses the direction, though correlation is undirected so no meaning.
	diag(mat[,-1]) = x # -1 here means no first col.
	diag(mat[-1,]) = x # first row skipped.
    # in this way we see that superdiagnal is the diagonal above the principal diagonal. And the sub, well, you can take a guess.
	return(mat)
}

#changes the superdiagonal and subdiagonal of a particular submatrix of a larger matrix
add_sigma_to_matrix <- function(total_matrix, row_start, nrows, x)
{
	extent = row_start + nrows # it'll be a square submatrix, i.e. nrows x nrows.
	extract_mat = total_matrix[row_start:extent, row_start:extent]
	mat = subsup_diag(extract_mat, x)
	total_matrix[row_start:extent, row_start:extent] <- mat
	return(total_matrix)
}

#updates the corrrelation matrix and adds the gene pairs that were modified to a class for subsequent extraction
update_corr_structure <- function(sigma, n_genes, variance, corr_factor, dcor = NULL, dcor_specific = NULL)
{
    # this sets the sigma matrix, which will store sigma .. the variance
	# this next value is updated each time, so the changes are applied sequentially,
    # after the last batch so to speak
    # wacth letters_unique it's not in the arguments, no, it's picked up globally. Usual R confusion.
	total_dc_genes = sigma[["total_dc_genes"]]
	sigma[["sigma_tot"]] = add_sigma_to_matrix(sigma[["sigma_tot"]],
		row_start = total_dc_genes + 1,
		nrows = n_genes, x = variance * corr_factor)
    # that's the number all set, yep! THis is mostly a labelling function.
    # following is just a labelling of the pairs. 
	if(!is.null(dcor)){
        # should be if, else's.
		if(dcor == "goc"){
			sigma[["goc_gene_pairs"]] = c(sigma[["goc_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if(dcor == "loc") {
			sigma[["loc_gene_pairs"]] = c(sigma[["loc_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
	}
    # also a labelling ...
	if(!is.null(dcor_specific)){
        # should be if, else's.
		if(dcor_specific == "+/0"){
			sigma[["pos_none_gene_pairs"]] = c(sigma[["pos_none_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if(dcor_specific == "+/-"){
			sigma[["pos_neg_gene_pairs"]] = c(sigma[["pos_neg_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if(dcor_specific == "0/+"){
			sigma[["none_pos_gene_pairs"]] = c(sigma[["none_pos_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if(dcor_specific == "0/-"){
			sigma[["none_neg_gene_pairs"]] = c(sigma[["none_neg_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if (dcor_specific == "-/+"){
			sigma[["neg_pos_gene_pairs"]] = c(sigma[["neg_pos_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		} else if(dcor_specific == "-/0"){
			sigma[["neg_none_gene_pairs"]] = c(sigma[["neg_none_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
	}
	sigma[["total_dc_genes"]] = total_dc_genes + n_genes
	return(sigma)
}

#for gene name generation, i.e. "aa", "ab", "ac", etc
amck_combine <- function (vecs, number) # Andy's combine function
{
    utils::combn(vecs, number, paste, collapse = "") # seldom used function from utils package
}

letters_unique0 = amck_combine(c(letters[1:26], letters[1:26]), 2)
letters_unique = sort(unique(letters_unique0)) # just a series of 2-letter labels, 676 of them.

# parameters for simulation data
n_non_expr_genes = 20
n_housekeeping_genes = 10
n_activated_genes = 30
n_tot_genes = n_non_expr_genes + n_housekeeping_genes + n_activated_genes
ngenesubset <- 2 # subset size for the 9 different correlation combos (+/+, +/= etc.) (sually all eqaul in size)

low_mean = 2
high_mean = 50
low_var = 50
high_var = 100
poscorrvalue = 0.5
negcorrvalue = -0.5 #if this is set too low, then the matrix will be non-positive definite

#create overall correlation matrices and classes to store them
sigma_tot_a <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0
	), class = "sigma")

sigma_tot_b <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0,
	loc_gene_pairs = vector(),
	goc_gene_pairs = vector(),
	pos_none_gene_pairs = vector(),
	pos_neg_gene_pairs = vector(),
	none_pos_gene_pairs = vector(),
	none_neg_gene_pairs = vector(),
	neg_pos_gene_pairs = vector(),
	neg_none_gene_pairs = vector()
	), class = "sigma")

# So a and b are the two conditions, and they have a sigma matrix each
# that simple. 

#A+, B+ = no differential correlation
# for ngeneset=2
# 50 0 and 50 0
# 0 50 and 0 50
# so the sigma is the same .. that's the way no correlation is obtained.
n_apbp <- ngenesubset # number_apos_bpos
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbp,
	variance = high_var, corr_factor = poscorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbp,
	variance = high_var, corr_factor = poscorrvalue, dcor = NULL)
# what happens?
# super and subdiagonals get value of fifty.

write.csv(sigma_tot_a[[1]], "sigma_tot_a0.csv", quote=F)
write.csv(sigma_tot_b[[1]], "sigma_tot_b0.csv", quote=F)

#A+, B= = loss of correlation
# 50 0 and 0 0
# 0 50 and 0 0
n_apbe <- ngenesubset # number a pos b equal
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbe,
	variance = high_var, corr_factor = poscorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbe,
	variance = high_var, corr_factor = 0, dcor = "loc", dcor_specific = "+/0")

#A+, B- = loss of correlation
# 50 0 and -50 0
# 0 50 and 0 -50
# actually you work these out 
n_apbm <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbm,
	variance = high_var, corr_factor = poscorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbm,
	variance = high_var, corr_factor = negcorrvalue, dcor = "loc",
	dcor_specific = "+/-")

#A=, B+ = gain of correlation
n_aebp <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebp,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebp,
	variance = high_var, corr_factor = poscorrvalue, dcor = "goc",
	dcor_specific = "0/+")

write.csv(sigma_tot_a[[1]], "sigma_tot_a1.csv", quote=F)
write.csv(sigma_tot_b[[1]], "sigma_tot_b1.csv", quote=F)

#A=, B- = loss of correlation
n_aebm <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = negcorrvalue, dcor = "loc",
	dcor_specific = "0/-")

# note clear oversight, he sets the n_genes to n_aebm for the rest too,
# doesn't matter because they are all 20, corrected.
#A-, B+ = gain of correlation
n_ambp <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_ambp,
	variance = high_var, corr_factor = negcorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_ambp,
	variance = high_var, corr_factor = poscorrvalue, dcor = "goc",
	dcor_specific = "-/+")

#A-, B= = gain of correlation
n_ambe <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_ambe,
	variance = high_var, corr_factor = negcorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_ambe,
	variance = high_var, corr_factor = 0, dcor = "goc",
	dcor_specific = "-/0")

#A-, B- = no differential correlation
n_ambm <- ngenesubset
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_ambm,
	variance = high_var, corr_factor = negcorrvalue, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_ambm,
	variance = high_var, corr_factor = negcorrvalue, dcor = NULL)

#A=, B= = no differential correlation
n_aebe = n_activated_genes - sigma_tot_a[["total_dc_genes"]]
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebe,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebe,
	variance = high_var, corr_factor = 0, dcor = NULL)

sigma_tot_a_mat = sigma_tot_a[["sigma_tot"]]
sigma_tot_b_mat = sigma_tot_b[["sigma_tot"]]

write.csv(sigma_tot_a[[1]], "sigma_tot_a2.csv", quote=F)
write.csv(sigma_tot_b[[1]], "sigma_tot_b2.csv", quote=F)

########################################
#set the means and variances of the genes

#set the variance
# Note usage of Gene expression favourite probdist rnbinom, but it's for the mean
# rnorm used later .. rnbinom mean for the different geens .. I suppose it makes sense.
var_diag = c(rep(high_var, n_activated_genes),
	rep(low_var, n_housekeeping_genes),
	rep(high_var, n_non_expr_genes))

#set the mean for the high expression genes
mu = rnbinom(n = (n_activated_genes + n_housekeeping_genes), mu = high_mean, size = 0.5)

#set the mean for the low expression genes
mu = c(mu, rnbinom(n = n_non_expr_genes, mu = low_mean, size = 0.5))

#if mu is less than one, set it to one
mu[mu < 1] = 1

diag(sigma_tot_a_mat) = var_diag
diag(sigma_tot_b_mat) = var_diag

#generate the simulation matrix
sim_data_a = as.matrix(MASS::mvrnorm(nsamps_per_cond, mu = mu, Sigma = sigma_tot_a_mat))
sim_data_b = as.matrix(MASS::mvrnorm(nsamps_per_cond, mu = mu, Sigma = sigma_tot_b_mat))

true_cases = c(sigma_tot_b$loc_gene_pairs, sigma_tot_b$goc_gene_pairs)
stop("o")
#run the ddcor differential correlation pipeline on the simulated data

print("starting the ddcor pipeline")

sim_data_a = t(sim_data_a)
sim_data_b = t(sim_data_b)

# this is then the expression matrix
sim_data_merge = cbind(sim_data_a, sim_data_b)
rownames(sim_data_merge) = letters_unique[1:nrow(sim_data_merge)]

conditions = c(rep("cond_a", ncol(sim_data_a)), rep("cond_b", ncol(sim_data_b)))

design_mat = model.matrix(~ conditions + 0)
labels = c("cond_a", "cond_b")
colnames(design_mat) = labels
npairs = (nrow(sim_data_merge)^2)/2 - nrow(sim_data_merge)

ddcor_start = Sys.time()
# ddcor_res = DGCA::ddcorAll(inputMat=sim_data_merge, design=design_mat, compare=labels, corr_cutoff=.95,
#                           adjust="none", corrType="pearson", nPerm=0, heatmap=F, nPairs="all")
# heatmap option is now heatmapPlot
# modify take out defaults based on DGCA manual which is:
# ddcorAll(inputMat, design, compare, inputMatB = NULL, splitSet = NULL,
#   impute = FALSE, corrType = "pearson", nPairs = "all",
#   sortBy = "zScoreDiff", adjust = "perm", nPerms = 10, classify = TRUE,
#   sigThresh = 1, corSigThresh = 0.05, heatmapPlot = FALSE,
#   color_palette = NULL, verbose = FALSE, plotFdr = FALSE,
#   corr_cutoff = 0.99, signType = "none", getDCorAvg = FALSE,
#   dCorAvgType = "gene_average", dCorAvgMethod = "median",
#   oneSidedPVal = FALSE, customize_heatmap = FALSE, heatmapClassic = FALSE,
#   corPower = 2, ...)
CairoPNG("ddcorhm.png", 1000, 1000)
ddcor_res = DGCA::ddcorAll(inputMat=sim_data_merge, design=design_mat, compare=labels, corr_cutoff=.95,
                           adjust="none", nPerm=0, heatmapPlot=T)
dev.off()
# the bottom part of heat map (which uses heatmap.2() is not comign out. Top part seems ok.
# actually for the 60 not 600 gene run, the hm does come out ... with borders.
# it's prob the borders that are making it all black.

print(paste0("time taken to run the ddcor pipeline with ", nsamps_per_cond, " samples..."))
ddcor_time = difftime(ddcor_start, Sys.time(), units = "secs")
print(difftime(ddcor_start, Sys.time()))

ddcor_res$combos = paste(ddcor_res$Gene1, ddcor_res$Gene2, sep = " ")
ddcor_res$true = ddcor_res$combos %in% true_cases

# function for plotting sens vs. spec.
plot_groupwise_ROC_ddcor <- function(res_pval = ddcor_res$pValDiff,
	pair_names = ddcor_res$combos, group_true, all_true = true_cases)
{
	group_true_vec = pair_names %in% group_true
	to_remove = pair_names %in% setdiff(all_true, group_true)
	tmp_df = data.frame(pvals = res_pval, true = group_true_vec)
	tmp_df = tmp_df[!to_remove, ]
	roc_list = roc(tmp_df$true, tmp_df$pvals, smooth = FALSE)
	return(roc_list)
}

print("starting ddcor ROC curve calculation")
# following function is defined up above.
full_roc = plot_groupwise_ROC_ddcor(group_true = true_cases)

print("full AUC for ddcor with nsamps_per_cond, and two conds ...")
print(full_roc)

ddcor_line = c("ddcor", nsamps_per_cond, seed, ddcor_time, as.numeric(full_roc$auc))
write.table(ddcor_line, "auc_results.txt", sep = "\t", quote=F, append=T, col.names=F, row.names=F)

roc_loc_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$loc_gene_pairs)
roc_goc_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$goc_gene_pairs)
roc_p0_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$pos_none_gene_pairs)
roc_pn_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$pos_neg_gene_pairs)
roc_0p_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$none_pos_gene_pairs)
roc_0n_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$none_neg_gene_pairs)
roc_np_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$neg_pos_gene_pairs)
roc_n0_only = plot_groupwise_ROC_ddcor(group_true = sigma_tot_b$neg_none_gene_pairs)

tiff(file = paste0("ddcor_classes_", nsamps_per_cond, "_samples", seed, ".tiff"),
		width = 4, height = 4, units = 'in', res = 300)
plot(roc_goc_only, col = cbPalette[1])
plot(roc_loc_only, add = TRUE, col = cbPalette[2])
plot(roc_p0_only, add = TRUE, col = cbPalette[3])
plot(roc_pn_only, add = TRUE, col = cbPalette[4])
plot(roc_0p_only, add = TRUE, col = cbPalette[5])
plot(roc_0n_only, add = TRUE, col = cbPalette[6])
plot(roc_np_only, add = TRUE, col = cbPalette[7])
plot(roc_n0_only, add = TRUE, col = cbPalette[8])
dev.off()
