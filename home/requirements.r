# Parallel Installation of Packages
options(Ncpus = parallel::detectCores())

# install all of the CRAN Packages
# can you comment out any of the packages you think are not necessary
install.packages('tidyverse', repos = "http://cran.us.r-project.org")
install.packages('devtools', repos = "http://cran.us.r-project.org")
install.packages('liqueueR', repos = "http://cran.us.r-project.org")

options(BioC_mirror = "http://bioconductor.org", repos = "http://cran.r-project.org")

# First Install all Bioconductor R Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("DESeq2", silent = TRUE)
BiocManager::install("EnhancedVolcano", silent = TRUE)
BiocManager::install("gprofiler2", silent = TRUE)
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene", silent = TRUE)

# Special Load-In ChIPseeker Dependencies
library(devtools)
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("YuLab-SMU/enrichplot")
BiocManager::install("ChIPseeker", silent = TRUE)