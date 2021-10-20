# unanchored-cnma

all_mortality.csv - data on all mortality from the systematic literature review done by Welton et al

all_mortality2.csv - same dataset as before but arms 1 and 2 of trial id 15 are combined because the treatments are the same - for compatibility with netmeta

application.R - code to obtain the results found in the Application section of the paper

make_jags_data.R - function to format input data object and create the design matrices needed to run the models in JAGS

arm_unanchored.R - JAGS code to run the novel arm-based Bayesian model introduced in section 3.2

arm_unanchored_decomposed.R - JAGS code to run the model in section 3.2 written in a faster implementation

contrast_unanchored.R - JAGS code to run the novel arm-based Bayesian model introduced in section 3.1

contrast_unanchored_decomposed.R - JAGS code to run the model in section 3.1 written in a faster implementation
