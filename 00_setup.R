# Shared setup: packages, paths, seed
pkgs <- c("tidyverse","janitor","naniar","mice","gtsummary","gt","survival","survminer",
          "broom","ranger","here","quarto")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(naniar); library(mice)
  library(gtsummary); library(survival); library(survminer); library(broom)
  library(ranger); library(here)
})
set.seed(2026)
dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("data","clean"), showWarnings = FALSE)

