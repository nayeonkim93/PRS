############################ INSTALL PACKAGE ####################################

# install.packages("dplyr")
# install.packages("data.table")
# install.packages("devtools", denpendencies = TRUE)
# library(devtools)
# install_github("tshmak/lassosum")

############################ LIBRARY UPLOAD #####################################

library(dplyr)
library(lassosum)
library(data.table)

################################## LASSOSUM #####################################

# Summary Statistic from PheWeb 
ss <- fread("/home/nayeonkim1/summary_statistics/UKBB_T2DB.txt")

# Define reference panel, test sample, and LDblock
ref.bfile <- "/home/lee7801/DATA/PRS/EURID_1000G_plink" # reference file, usually 1000 Genome
test.bfile <- "" # test.bfile is the plink file you used for model training.
LDblocks <- "EUR.hg19"  # LDblocks can be ASN, EUR, etc. Details can be found at https://github.com/tshmak/lassosum.

# Phenotype and covariate
pheno <- fread("/home/nayeonkim1/pheno/KBN_WHOLE_T2DB.txt") # phenotype
covar <- fread("/home/nayeonkim1/pheno/KBN_WHOLE_covar.txt") # covariate

# Split into validation set and test set
train_set <- fread("/home/nayeonkim1/pheno/train.txt")
test_set <- fread("/home/nayeonkim1/pheno/test.txt")

# Correlation calculation
cor <- p2cor(p = ss$pval, n = length(ss$pos), sign=ss$beta)

#' A1 Alternative allele 
#' A2 Reference allele 

# Lassosum 
out1 <- lassosum.pipeline(cor, chr=ss$`#chrom`, pos=ss$pos, 
                          A1=ss$alt, A2=ss$ref, 
                          ref.bfile=ref.bfile, test.bfile=test.bfile, 
                          remove.test = test_set,
                          LDblocks=LDblocks)

# Validation with sample 1
v1 <- validate(out1, pheno = pheno, covar = covar,  keep = train_set) 

# Validation with sample 2 (with the best lambda and s)
out2 <- subset(out1, s=v1$best.s, lambda=v1$best.lambda)
v2 <- validate(out2, pheno = pheno, covar = covar, keep= test_set) 

