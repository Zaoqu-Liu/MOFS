##----------------------------------------------------------------------
## This script is designed for loading R packages in the main shiny app
##----------------------------------------------------------------------
## Maintainer: Mingjie Wang (huzai920621@126.com)
## Log: 
## 2022-02-04, initial build


#### Set packages resources ####
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  

message("[++] It will cost a while for the first time.")


# Packages from cran and bioconductor
pkgs.to.check <- c(
  "shiny",
  "shinyjs",
  "bslib",
  "readr",
  "stringr",
  "slickR",
  "markdown",
  "DT",
  "shinyWidgets",
  "shinyFeedback",
  "neuralnet"
)

for (pkg.name in pkgs.to.check) {
  if (!requireNamespace(pkg.name)) {
    message(paste("[++] Installing", pkg.name))
    tryCatch(
      pacman::p_install(pkg.name,character.only = T),
      error = function(e) {
        message(paste("[++] ERROR: R packages ["), pkg.name, "] can NOT be installed.")
      }
    )
  }
  library(pkg.name,character.only = T)
}

