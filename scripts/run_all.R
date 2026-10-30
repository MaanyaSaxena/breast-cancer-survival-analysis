for (s in c("01_clean","02_explore","03_survival","04_models"))
  source(here::here("scripts", paste0(s, ".R")), echo = FALSE)
quarto::quarto_render(here::here("reports","report.qmd"))
