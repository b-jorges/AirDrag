libraries <- c("cowplot","tidyverse","grid","gridExtra","lme4","lmerTest","ggpubr","apastats","rstatix")
# Where is this script?
invisible(lapply(libraries, function(x) {
  if(!require(x, character.only = T, quietly = T)) {
    install.packages(x)
    require(x, character.only = T)
  }
}
))
rm(libraries)

Where_Am_I <- function(path=T){
  if (path == T){
    dirname(rstudioapi::getSourceEditorContext()$path)
  }
  else {
    rstudioapi::getSourceEditorContext()$path
  }
}
theme_set(theme_minimal(10))
scale_color_discrete <- function(...) {
  ggthemes::scale_color_colorblind(...)
}
scale_fill_discrete <- function(...) {
  ggthemes::scale_fill_colorblind(...)
}
source("Funs.R")