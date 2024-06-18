.onLoad <- function(libname, pkgname){
  msg <- paste0(pkgname, " v", "0.1.0", "\n\n","For help: https://grswsci.top/","\n\n")
  citation <- paste0("If you use ", pkgname, " in published research, please thank:\n",
                      "Yuyao Liu, Bioinformatics R&D Department, Hefei GuangRe Biotechnology Co., Ltd, Hefei, China","\n\n")
  packageStartupMessage(paste0(msg, citation))
  options(timeout = 10000)
  browseURL("https://grswsci.top/", browser = getOption("browser"),encodeIfNeeded = FALSE)
}

