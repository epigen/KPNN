pkgs <- installed.packages()

for(pkgx in c("ggplot2", "data.table")){
  if(!pkgx %in% row.names(pkgs)){
    install.packages(pkgx, repos = "https://cloud.r-project.org")
  }
}
