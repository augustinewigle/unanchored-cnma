#####################################
# Script for results in Application
# By Augustine Wigle
# March 15 2022
#####################################

# libraries
library(dplyr)
library(tidyr)
library(runjags)

source("make_jags_data.R")

# Read in mortality data

allmort <- read.csv("all_mortality.csv")
allmort <- allmort[with(allmort, order(study, -usual)),] # reorder so first arm in each study is always usual care

# Mortality data in format for Rucker model - combined arms 1 and 2 from study id 15
allmort2 <- read.csv("all_mortality2.csv")

# Create jags data
jags_data <- make_jags_data(allmort, studycol = "study", components = c("usual", "edu", "beh", "cog",  "relax","support"))

# Running the models with JAGS -------------------------------------------------------------------------------------------

# Anchored, arm-level - Welton model-------------------------------------------
set.seed(22)
inits_anchor <- list(list(alpha = rep(1,jags_data$arm_un$nstudy),
                          sigma = 0.1,
                          d1 = rep(0.1, jags_data$arm_un$ncomp-1)),
                     list(alpha = runif(jags_data$arm_un$nstudy, -1,1),
                          sigma = 1,
                          d1 = runif(jags_data$arm_un$ncomp-1, -0.5, 0.5)))

# using univariate conditional distributions - this is equivalent but faster than sampling from the multivariate error distribution
anchor_samples <- run.jags(model = "arm_anchored_decomposed.R", # if you want to use the multivariate errors, change to "arm_anchored.R"
                        data = jags_data$arm_anchor,
                        inits = inits_anchor,
                        monitor = c("d1", "sigma"),
                        n.chains = 2,
                        burnin = 10000,
                        sample = 50000) # warning for unused variable Sigma is ok, Sigma only needed if multivariate version used
anchor_dic <- extract(anchor_samples, "DIC")
results_anchor <- summary(anchor_samples)

# Unanchored, arm-level ------------------------------------------------

# initial values for two chains
inits_un <- list(list(alpha = rep(1,jags_data$arm_un$nstudy),
                      sigma = 0.1,
                      d = rep(0.1, jags_data$arm_un$ncomp)),
                 list(alpha = runif(jags_data$arm_un$nstudy, -1,1),
                      sigma = 1,
                      d = runif(jags_data$arm_un$ncomp, -0.5, 0.5)))

# Using the conditional distributions
arm_samples <- run.jags(model = "arm_unanchored_decomposed.R",
                       data = jags_data$arm_un,
                       inits = inits_un,
                       monitor = c("d", "sigma"),
                       n.chains = 2,
                       burnin = 20000,
                       sample = 50000) # warning for unused variable - ok

un_dic <- extract(arm_samples, "DIC")
results_un <- summary(arm_samples)

# Unanchored, contrast-based ----------------------------------------------

# initial values for 2 chains
inits_con <- list(list(sigma = 0.1,
                       d = rep(0.1, jags_data$arm_un$ncomp)),
                  list(sigma = 1,
                       d = runif(jags_data$arm_un$ncomp, -0.5, 0.5)))

# Using the conditional distributions
con_samples <- run.jags(model = "contrast_unanchored_decomposed.R",
                        data = jags_data$contrast_un,
                        inits = inits_con,
                        monitor = c("d", "sigma"),
                        n.chains = 2,
                        burnin = 10000,
                        sample = 40000) # warning for unused variable

con_dic <- extract(con_samples, "DIC")

results_con <- summary(con_samples)

# Getting frequentist results using netmeta package --------------------------------------

library(netmeta)

names_id <- c("usual", "edu", "beh", "cog",  "relax","support")
trt_names <- vector(mode = "character", length = nrow(allmort2))

for(i in 1:nrow(allmort2)) {
  
  trt_names[i] <- paste(names_id[which(allmort2[i,2:7] >0)], collapse = "+")
  
}

# Use pairwise function to prepare the data
p1 <- pairwise(treat = trt_names, event = allmort2$r, n = allmort2$n, studlab = allmort2$study, incr = 0.5, warn = T, allstudies = T)

# do NMA and then CNMA
nm <- netmeta(TE = p1$TE, seTE = p1$seTE, treat1 = p1$treat1, treat2 = p1$treat2, studlab = p1$studlab, sm = "OR", reference.group = "usual")
cnma <- netcomb(nm)

rucker_means <- cnma$TE.random[c("edu", "beh", "cog", "support"), "usual"]
rucker_lowers <- cnma$lower.random[c("edu", "beh", "cog", "support"), "usual"]
rucker_uppers <- cnma$upper.random[c("edu", "beh", "cog", "support"), "usual"]

# Plotting the results -------------------------------------------------------------------
library(forestplot)
# # Results for anchored Welton model - from Welton et al 2009
# welton_df <- data.frame(mean = c(0.29, -0.58, -0.01, -0.38, 0.21),
#                         lower = c(-0.27, -1.13, -0.52, -1.16, -0.66),
#                         upper = c(0.85, -0.05, 0.45, 0.37, 1.06),
#                         labeltext = c("Edu",
#                                       "Beh",
#                                       "Cog",
#                                       "Rel",
#                                       "Sup"),
#                         group = "2.1 (Welton)")

anchor_mcmc <- as.matrix(anchor_samples$mcmc)
anchor_df <- data.frame(mean = results_anchor[-6,"Mean"],
                        upper = results_anchor[-6,"Upper95"],
                        lower = results_anchor[-6,"Lower95"],
                        labeltext = c("Edu",
                                      "Beh",
                                      "Cog",
                                      "Rel",
                                      "Sup"))
anchor_df$group = "Bayes Arm Anchored"

# make d_usual,x for unanchored models
un_mcmc <- as.matrix(arm_samples$mcmc)
un_edu <- un_mcmc[,"d[2]"] - un_mcmc[,"d[1]"]
un_beh <- un_mcmc[,"d[3]"] - un_mcmc[,"d[1]"]
un_cog <- un_mcmc[,"d[4]"] - un_mcmc[,"d[1]"]
un_rel <- un_mcmc[,"d[5]"] - un_mcmc[,"d[1]"]
un_sup <- un_mcmc[,"d[6]"] - un_mcmc[,"d[1]"]
un_results <- cbind(un_edu, un_beh, un_cog, un_rel, un_sup)
un_df <- data.frame(mean = apply(un_results, 2, mean),
                    lower = apply(un_results, 2, quantile, probs = 0.025),
                    upper = apply(un_results, 2, quantile, probs = 0.975),
                    labeltext = c("Edu",
                                  "Beh",
                                  "Cog",
                                  "Rel",
                                  "Sup"))
un_df$group <- "Bayes Arm Unanchored"

con_mcmc <- as.matrix(con_samples$mcmc)
con_edu <- con_mcmc[,"d[2]"] - con_mcmc[,"d[1]"]
con_beh <- con_mcmc[,"d[3]"] - con_mcmc[,"d[1]"]
con_cog <- con_mcmc[,"d[4]"] - con_mcmc[,"d[1]"]
con_rel <- con_mcmc[,"d[5]"] - con_mcmc[,"d[1]"]
con_sup <- con_mcmc[,"d[6]"] - con_mcmc[,"d[1]"]
con_results <- cbind(con_edu, con_beh, con_cog, con_rel, con_sup)
con_df <- data.frame(mean = apply(con_results, 2, mean),
                     lower = apply(con_results, 2, quantile, probs = 0.025),
                     upper = apply(con_results, 2, quantile, probs = 0.975),
                     labeltext = c("Edu",
                                   "Beh",
                                   "Cog",
                                   "Rel",
                                   "Sup"))
con_df$group <- "Bayes Con Unanchored"

# Rucker model
rucker_df <- data.frame(mean = c(rucker_means, NA_real_), lower = c(rucker_lowers, NA_real_), upper = c(rucker_uppers, NA_real_))
rucker_df$labeltext <- c("Edu",
                         "Beh",
                         "Cog",
                         "Sup",
                         "Rel")
colnames(rucker_df) <- c("mean", "lower", "upper", "labeltext")
rucker_df$group <- "Freq Con Unanchored"
rucker_df <- rucker_df[c(1,2,3,5,4),]

# Plotting

# newdf <- as_tibble(rbind(con_df, un_df, welton_df, rucker_df))
newdf <- as_tibble(rbind(con_df, un_df, anchor_df, rucker_df))
newdf$labeltext <- as.character(newdf$labeltext)

newdf %>% group_by(group)%>%forestplot(shapes_gp = fpShapesGp(box = c("navyblue", "cornflowerblue", "darkred", "lightpink") %>% lapply(function(x) gpar(fill = x, col = "#555555")),
                                                              default = gpar(vertices = TRUE)),
                                       fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI, fpDrawNormalCI),
                                       legend_args = fpLegend(pos = list(x = .90, y = 0.5),
                                                              gp = gpar(col = "black", fill = "white"),
                                                              title = "Models"),
                                       ci.vertices = TRUE,
                                       ci.vertices.height = 0.05,
                                       boxsize = .135,
                                       title = "Log-odds ratio of components vs. Usual care",
                                       xlab = "Log odds ratio vs Usual care")

# compare the DICs
un_dic
anchor_dic
