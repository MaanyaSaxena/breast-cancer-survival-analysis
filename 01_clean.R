# 01_clean.R -- import, join, validate, clean, document missingness
source(here::here("scripts","00_setup.R"))


pat_file <- here("data","raw","data_clinical_patient.txt")
smp_file <- here("data","raw","data_clinical_sample.txt")
stopifnot(file.exists(pat_file), file.exists(smp_file))

# cBioPortal files start with 4 metadata lines beginning with '#'
patient <- read_tsv(pat_file, comment = "#", show_col_types = FALSE) |> clean_names()
sample  <- read_tsv(smp_file, comment = "#", show_col_types = FALSE) |> clean_names()

log <- character()
add_log <- function(x) { log <<- c(log, x); message(x) }
add_log(sprintf("Raw patients: %d rows; raw samples: %d rows", nrow(patient), nrow(sample)))

stopifnot(!anyDuplicated(patient$patient_id))
if (anyDuplicated(sample$patient_id))
  add_log("Samples had duplicate patient IDs; keeping first sample per patient")
sample <- distinct(sample, patient_id, .keep_all = TRUE)

# Unmatched IDs before joining
add_log(sprintf("Unmatched IDs: %d in patient file only, %d in sample file only",
                nrow(anti_join(patient, sample, by = "patient_id")),
                nrow(anti_join(sample, patient, by = "patient_id"))))

dat <- inner_join(patient, sample, by = "patient_id", suffix = c("", "_smp"))
add_log(sprintf("After join: %d rows", nrow(dat)))
stopifnot(nrow(dat) == n_distinct(dat$patient_id))

# ---- Validation of labels and numeric ranges (fail loudly) ----
chk <- function(ok, msg) if (!isTRUE(ok)) stop("Validation failed: ", msg, call. = FALSE)
in_rng <- function(x, lo, hi) all(is.na(x) | (x >= lo & x <= hi))
chk(all(dat$os_status %in% c("0:LIVING", "1:DECEASED", NA)), "unexpected os_status labels")
chk(all(dat$grade %in% c(1, 2, 3, NA)), "grade outside 1-3")
chk(all(as.numeric(dat$tumor_stage) %in% c(0, 1, 2, 3, 4, NA)), "tumor_stage outside 0-4")
chk(all(dat$er_status  %in% c("Positive", "Negative", NA)), "unexpected ER status labels")
chk(all(dat$pr_status  %in% c("Positive", "Negative", NA)), "unexpected PR status labels")
chk(all(dat$her2_status %in% c("Positive", "Negative", NA)), "unexpected HER2 status labels")
for (v in c("chemotherapy", "hormone_therapy", "radio_therapy"))
  chk(all(toupper(dat[[v]]) %in% c("YES", "NO", NA)), paste("unexpected labels in", v))
chk(in_rng(as.numeric(dat$age_at_diagnosis), 0, 120), "implausible age")
chk(in_rng(as.numeric(dat$tumor_size), 0, 300), "implausible tumor size (mm)")
chk(in_rng(as.numeric(dat$lymph_nodes_examined_positive), 0, 100), "implausible positive node count")
chk(in_rng(as.numeric(dat$os_months), -Inf, 1000), "implausible survival time")
add_log("Validation passed: labels, stage/grade ranges, and numeric ranges")

# Outcome: all-cause overall survival. os_status is "0:LIVING" / "1:DECEASED"
dat <- dat |> mutate(event = as.integer(str_sub(os_status, 1, 1)), time = as.numeric(os_months))
bad <- dat |> filter(is.na(time) | time <= 0 | is.na(event))
add_log(sprintf("Dropped %d patients with missing/non-positive survival time or status (missing: %d; non-positive time: %d)",
                nrow(bad), sum(is.na(bad$time) | is.na(bad$event)), sum(!is.na(bad$time) & bad$time <= 0)))
add_log("Selection rule: the outcome cannot be analysed without a valid positive survival time, so these patients are excluded (not imputed)")
dat <- anti_join(dat, bad, by = "patient_id")

yn <- function(x) factor(str_to_title(x), levels = c("No","Yes"))
clean <- dat |>
  transmute(
    patient_id, time, event,
    age         = as.numeric(age_at_diagnosis),
    tumor_size  = as.numeric(tumor_size),
    stage       = factor(case_when(as.numeric(tumor_stage) <= 1 ~ "0-1",
                                   as.numeric(tumor_stage) == 2 ~ "2",
                                   as.numeric(tumor_stage) >= 3 ~ "3-4"),
                         levels = c("0-1", "2", "3-4")),
    grade       = factor(grade, levels = c(1, 2, 3)),
    nodes_pos   = as.numeric(lymph_nodes_examined_positive),
    er_status   = factor(er_status),
    pr_status   = factor(pr_status),
    her2_status = factor(her2_status),
    chemo       = yn(chemotherapy),
    hormone_tx  = yn(hormone_therapy),
    radio_tx    = yn(radio_therapy)
  ) |>
  mutate(
    age_group = cut(age, c(0, 50, 65, Inf), labels = c("<50","50-64","65+"), right = FALSE),
    hr_status = factor(if_else(er_status == "Positive" | pr_status == "Positive",
                               "HR-positive", "HR-negative"),
                       levels = c("HR-positive","HR-negative"))
  )
# If a column is not found, inspect names(patient) and names(sample): cBioPortal exports vary slightly.

# Missingness
write_csv(miss_var_summary(clean), here("data","clean","missingness_summary.csv"))
ggsave(here("figures","missingness.png"),
       gg_miss_var(clean, show_pct = TRUE) + labs(title = "Missing values by variable"),
       width = 7, height = 5, dpi = 200)

key <- c("age","tumor_size","stage","grade","nodes_pos","er_status","pr_status","her2_status")
miss_check <- clean |>
  mutate(any_missing = if_any(all_of(key), is.na)) |>
  group_by(any_missing) |>
  summarise(n = n(), event_rate = mean(event), median_followup = median(time), .groups = "drop")
print(miss_check); write_csv(miss_check, here("data","clean","missingness_vs_outcome.csv"))
add_log(sprintf("Complete cases on key variables: %d of %d", sum(complete.cases(clean[key])), nrow(clean)))

write_csv(clean, here("data","clean","metabric_clean.csv"))
writeLines(log, here("data","clean","cleaning_log.txt"))
