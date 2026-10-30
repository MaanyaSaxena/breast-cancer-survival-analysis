# 02_explore.R -- summary table + descriptive charts
source(here::here("scripts","00_setup.R"))
clean <- read_csv(here("data","clean","metabric_clean.csv"), show_col_types = FALSE) |>
  mutate(across(c(stage, grade), factor))

tbl <- clean |>
  select(age, age_group, tumor_size, stage, grade, nodes_pos, er_status, pr_status,
         her2_status, hr_status, chemo, hormone_tx, radio_tx, event) |>
  mutate(event = factor(event, labels = c("Censored","Died"))) |>
  tbl_summary(by = event, missing_text = "Missing") |> add_overall() |> add_p()
gt::gtsave(as_gt(tbl), here("figures","table1.html"))

p1 <- ggplot(clean, aes(age)) + geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  labs(title = "Age at diagnosis", x = "Age (years)", y = "Patients")
p2 <- ggplot(drop_na(clean, grade), aes(grade, fill = factor(event, labels = c("Censored","Died")))) +
  geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) +
  labs(title = "Outcome by tumor grade", x = "Grade", y = NULL, fill = NULL)
p3 <- ggplot(drop_na(clean, tumor_size), aes(tumor_size)) +
  geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
  labs(title = "Tumor size", x = "mm", y = "Patients")
ggsave(here("figures","eda_age.png"), p1, width = 6, height = 4, dpi = 200)
ggsave(here("figures","eda_grade.png"), p2, width = 6, height = 4, dpi = 200)
ggsave(here("figures","eda_size.png"), p3, width = 6, height = 4, dpi = 200)

# Median follow-up via reverse Kaplan-Meier
fu <- survfit(Surv(time, 1 - event) ~ 1, data = clean)
cat("Median follow-up (months):", surv_median(fu)$median, "\n")

