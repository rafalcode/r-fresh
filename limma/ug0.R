#!/usr/bin/env Rscript
# this script does what?
# from the user's guide
# page 44
# library(ggplot2)
library(Cairo)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()

# FileName  Strain Treatment
#     File1   WT          U
#     File2   WT           S
#         File3    Mu         U
#         File4    Mu          S
#             File5    Mu          S

phe <- data.frame(Strain=factor(c("WT", "WT", "Mu", "Mu", "Mu"), levels=c("WT", "Mu")),
                  Treat=factor(c("U", "S", "U", "S", "S"), levels=c("U", "S")))
des <- model.matrix(~Strain*Treat, data=phe)

# OK, so that pretty much reflects what the userguide says
# though you have to trust, because intercept isn't given the comparison, it's assumed you know. 
# Coefficient         Comparison              Interpretation
# Intercept           WT.U                    Baseline level of unstimulated WT
# StrainMu            Mu.U-WT.U               Difference between unstimulated strains
# TreatmentS          WT.S-WT.U               Stimulation effect for WT
# StrainMu:TreatmentS (Mu.S-Mu.U)-(WT.S-WT.U) Interaction

# To commentate, these are the coefficients of the linear model.Strain comes first, so that's the first difference, and it's done in terms of the first level of the second factor


