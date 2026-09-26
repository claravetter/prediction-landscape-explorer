#!/usr/bin/env Rscript
# Restore only software dependencies; no study input is read and no inference runs.
args<-commandArgs(TRUE)
if("--help"%in%args){cat("Rscript scripts/restore_environment.R [LIBRARY_PATH]\nDefaults to .Rlibrary in this repository. Requires network access for missing packages.\n");quit(status=0)}
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
project<-normalizePath(file.path(dirname(script),".."))
library_path<-if(length(args))args[1] else file.path(project,".Rlibrary")
dir.create(library_path,recursive=TRUE,showWarnings=FALSE)
library_path<-normalizePath(library_path)
.libPaths(c(library_path,.Library))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",lib=library_path)
renv::restore(project=project,library=library_path,lockfile=file.path(project,"renv.lock"),prompt=FALSE)
lock<-renv::lockfile_read(file.path(project,"renv.lock"))
for(pkg in names(lock$Packages)){
 # Recommended packages may be provided by this R installation. Check their
 # exact DESCRIPTION versions as well, without falling back to user libraries.
 actual<-as.character(utils::packageDescription(pkg,lib.loc=c(library_path,.Library))$Version)
 if(!identical(actual,lock$Packages[[pkg]]$Version))stop("Dependency version mismatch: ",pkg)
}
cat("PASS: restored",length(lock$Packages),"locked R packages to",library_path,"\n")
