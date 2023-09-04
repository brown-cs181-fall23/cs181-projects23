# DO NOT MODIFY THIS FILE #

# List of packages to check
packages_to_check <- c(
  "tidyverse",
  "devtools",
  "optparse",
  "pheatmap",
  "liqueueR",
  "VennDiagram",
  "DESeq2",
  "EnhancedVolcano",
  "GenomicRanges",
  "gprofiler2",
  "ChIPseeker"
)
# Create a function to check if a package is installed and loadable
check_package <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    message(paste("Package ",
    package_name, " is not installed or cannot be loaded."))
    return(FALSE)
  } else {
    message(paste("Package ",
    package_name, " is installed and can be loaded."))
    return(TRUE)
  }
}

# Check each package in the list
for (package in packages_to_check) {
  if (!check_package(package)) {
    message(paste("ERROR: ",
    package, " NOT installed. Please contact course staff."))
  }
}