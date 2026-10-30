# 04_models.R -- Cox vs random survival forest on held-out data + risk groups
source(here::here("scripts","00_setup.R"))
clean <- read_csv(here("data","clean","metabric_clean.csv"), show_col_types = FALSE) |>
  mutate(across(c(stage, grade), factor), across(where(is.character), factor))

vars <- c("age","tumor_size","stage","grade","nodes_pos","er_status","her2_status","chemo","hormone_tx","radio_tx")
d <- clean |> select(patient_id, time, event, all_of(vars))

set.seed(2026)
test_id <- d |> group_by(event) |> slice_sample(prop = 0.3) |> pull(patient_id)   # stratified 70/30
d <- d |> mutate(is_test = patient_id %in% test_id)

imp <- mice(d |> select(all_of(vars)), m = 1, maxit = 10, ignore = d$is_test, printFlag = FALSE, seed = 2026)
dd <- bind_cols(d |> select(patient_id, time, event, is_test), complete(imp, 1))
train <- filter(dd, !is_test); test <- filter(dd, is_test)
cat("Train:", nrow(train), " Test:", nrow(test), "\n")

num <- c("age","tumor_size","nodes_pos")
mu <- sapply(train[num], mean); sdv <- sapply(train[num], sd)
scale_df <- function(x) { for (v in num) x[[v]] <- (x[[v]] - mu[[v]]) / sdv[[v]]; x }
train_s <- scale_df(train); test_s <- scale_df(test)

f <- as.formula(paste("Surv(time, event) ~", paste(vars, collapse = " + ")))
cox <- coxph(f, data = train_s)
lp_train <- predict(cox, newdata = train_s, type = "lp")
lp_test  <- predict(cox, newdata = test_s,  type = "lp")

rsf <- ranger(f, data = train, num.trees = 500, mtry = 4, min.node.size = 15,
              importance = "permutation", seed = 2026)
rsf_risk_test  <- rowSums(predict(rsf, data = test)$chf)   # ensemble cumulative hazard = risk score
rsf_risk_train <- rowSums(rsf$chf)                          # out-of-bag for training rows
# Check train (OOB) and test score distributions are comparable
print(summary(rsf_risk_train)); print(summary(rsf_risk_test))

cidx <- function(risk, d) {
  cc <- concordance(Surv(time, event) ~ risk, data = d, reverse = TRUE)  # higher risk = earlier death
  tibble(C = cc$concordance, SE = sqrt(cc$var))
}
res <- bind_rows(bind_cols(model = "Cox PH", cidx(lp_test, test)),
                 bind_cols(model = "Random survival forest", cidx(rsf_risk_test, test))) |>
  mutate(lower = C - 1.96 * SE, upper = C + 1.96 * SE)
print(res); write_csv(res, here("data","clean","model_comparison.csv"))

ggsave(here("figures","model_comparison.png"),
  ggplot(res, aes(model, C, ymin = lower, ymax = upper, colour = model)) + geom_pointrange(size = 0.8) +
    geom_hline(yintercept = 0.5, linetype = 3) + coord_cartesian(ylim = c(0.5, 0.85)) +
    guides(colour = "none") + theme_minimal() + labs(title = "Held-out C-index (95% CI)", x = NULL, y = "C-index"),
  width = 6, height = 4, dpi = 200)

vi <- tibble(variable = names(rsf$variable.importance), importance = rsf$variable.importance)
ggsave(here("figures","rsf_importance.png"),
  ggplot(vi, aes(reorder(variable, importance), importance)) + geom_col(fill = "steelblue") +
    coord_flip() + theme_minimal() + labs(title = "RSF permutation importance", x = NULL, y = "Importance"),
  width = 6, height = 4, dpi = 200)

# ---- Risk groups for BOTH models (no model chosen using test results) ----

counts <- list()
make_groups <- function(r_tr, r_te, label, file) {
  cuts <- quantile(r_tr, c(1/3, 2/3))
  d <- test
  d$risk_group <- cut(r_te, c(-Inf, cuts, Inf), labels = c("Low","Medium","High"))
  stopifnot(!anyNA(d$risk_group))                       # every test patient classified
  counts[[label]] <<- d |> group_by(risk_group) |>
    summarise(patients = n(), deaths = sum(event), .groups = "drop") |> mutate(model = label)
  fit <- surv_fit(Surv(time, event) ~ risk_group, data = d)
  p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                  palette = c("#2a9d8f","#e9c46a","#e76f51"), xlab = "Months",
                  ylab = "Overall survival", title = paste("Test-set survival by risk group:", label),
                  ggtheme = theme_minimal())
  png(here("figures", file), width = 8, height = 7, units = "in", res = 200, bg = "white")
  print(p); dev.off()
  tidy(coxph(Surv(time, event) ~ risk_group, data = d), exponentiate = TRUE, conf.int = TRUE) |>
    mutate(model = label)
}
hr <- bind_rows(
  make_groups(lp_train,       lp_test,       "Cox PH",                 "risk_groups_cox.png"),
  make_groups(rsf_risk_train, rsf_risk_test, "Random survival forest", "risk_groups_rsf.png"))
write_csv(hr, here("data","clean","risk_group_hr.csv")); print(hr)
group_counts <- bind_rows(counts); print(group_counts)
write_csv(group_counts, here("data","clean","risk_group_counts.csv"))
