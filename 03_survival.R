# 03_survival.R -- Kaplan-Meier + log-rank, Cox model, forest plot, PH check
source(here::here("scripts","00_setup.R"))
clean <- read_csv(here("data","clean","metabric_clean.csv"), show_col_types = FALSE) |>
  mutate(across(c(stage, grade), factor), across(where(is.character), factor))

km_plot <- function(var, title) {
  d <- clean |> drop_na(all_of(var))
  fit <- surv_fit(as.formula(paste("Surv(time, event) ~", var)), data = d)
  p <- ggsurvplot(fit, data = d, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                  xlab = "Months", ylab = "Overall survival", title = title,
                  legend.title = "", ggtheme = theme_minimal())
  png(here("figures", paste0("km_", var, ".png")), width = 8, height = 7, units = "in", res = 200, bg = "white")
  print(p); dev.off()
  tibble(variable = var, p_logrank = surv_pvalue(fit)$pval)
}
km_res <- bind_rows(
  km_plot("stage", "Survival by tumor stage"),
  km_plot("grade", "Survival by tumor grade"),
  km_plot("age_group", "Survival by age group"),
  km_plot("hr_status", "Survival by hormone receptor status"),
  km_plot("her2_status", "Survival by HER2 status"))
write_csv(km_res, here("data","clean","logrank_results.csv")); print(km_res)

# ---- Cox model: complete-case (primary) ----
cox_vars <- c("age","tumor_size","stage","grade","nodes_pos","er_status","her2_status","chemo","hormone_tx","radio_tx")
cc <- clean |> select(time, event, all_of(cox_vars)) |> drop_na()
cat("Complete cases for Cox:", nrow(cc), "\n")
cox_f <- as.formula(paste("Surv(time, event) ~", paste(cox_vars, collapse = " + ")))
cox_cc <- coxph(cox_f, data = cc)
hr_cc <- tidy(cox_cc, exponentiate = TRUE, conf.int = TRUE) |> mutate(model = "Complete-case")

# ---- Sensitivity: multiple imputation (m = 5), pooled with Rubin's rules ----
imp_dat <- clean |> select(time, event, all_of(cox_vars))
imp_dat$cumhaz <- nelsonaalen(imp_dat, time, event)   # standard: include outcome info when imputing for inference
imp <- mice(imp_dat |> select(-time), m = 5, printFlag = FALSE, seed = 2026)
imp_long <- complete(imp, "all")
fits <- lapply(imp_long, function(x) coxph(cox_f, data = bind_cols(imp_dat[c("time","event")], x[cox_vars])))
hr_mi <- summary(pool(as.mira(fits)), conf.int = TRUE, exponentiate = TRUE) |> as_tibble() |>
  transmute(term, estimate, conf.low = `2.5 %`, conf.high = `97.5 %`, p.value, model = "Multiple imputation")

hr_all <- bind_rows(hr_cc |> select(term, estimate, conf.low, conf.high, p.value, model), hr_mi)
write_csv(hr_all, here("data","clean","cox_hazard_ratios.csv"))

fp <- ggplot(hr_all, aes(estimate, term, xmin = conf.low, xmax = conf.high, colour = model)) +
  geom_pointrange(position = position_dodge(width = 0.6)) + geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() + theme_minimal() +
  labs(title = "Adjusted hazard ratios (95% CI)", x = "Hazard ratio (log scale)", y = NULL, colour = NULL)
ggsave(here("figures","forest_plot.png"), fp, width = 8, height = 6, dpi = 200)

# ---- Proportional hazards check (Schoenfeld residuals) ----
zph <- cox.zph(cox_cc); print(zph)
write_csv(as_tibble(zph$table, rownames = "term"), here("data","clean","ph_test.csv"))
png(here("figures","schoenfeld.png"), width = 2600, height = 2000, res = 150); print(ggcoxzph(zph)); dev.off()

# Remedy 1: stratify on ER status only
f_er <- update(cox_f, . ~ . - er_status + strata(er_status))
cox_er <- coxph(f_er, data = cc); zph_er <- cox.zph(cox_er); print(zph_er)
write_csv(as_tibble(zph_er$table, rownames = "term"), here("data","clean","ph_test_er_strat.csv"))
write_csv(tidy(cox_er, exponentiate = TRUE, conf.int = TRUE), here("data","clean","cox_stratified.csv"))

# Remedy 2: time-varying age effect
cox_tt <- coxph(update(f_er, . ~ . + tt(age)), data = cc, tt = function(x, t, ...) x * log(t))
print(summary(cox_tt))
write_csv(tidy(cox_tt, exponentiate = FALSE, conf.int = TRUE), here("data","clean","cox_tt_age.csv"))

# Sensitivity: truncate follow-up at 10 years (120 months)
cc10 <- cc |> mutate(event = if_else(time > 120, 0L, event), time = pmin(time, 120))
cox_10 <- coxph(cox_f, data = cc10); print(cox.zph(cox_10))
write_csv(tidy(cox_10, exponentiate = TRUE, conf.int = TRUE), here("data","clean","cox_10yr.csv"))
