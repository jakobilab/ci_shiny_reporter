# List of required packages
packages <- c(
  "shiny", "bslib", "dplyr", "ggplot2", "ggExtra", "readr",
  "tidyr", "patchwork", "e1071", "forcats", "ggpubr",
  "fmsb", "scales", "shinydashboard", "shinythemes",
  "stringr", "broom", "RColorBrewer", "tibble"
)

missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)
lapply(packages, library, character.only = TRUE)


# ── Data Loading ──────────────────────────────────────────────────────────────
read_and_fix <- function(path, type_label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(across(is.character, ~na_if(.x, ""))) %>%
    mutate(type = type_label) %>%
    mutate(across(where(is.logical) | where(is.character), as.character))
}

open_source  <- rbind(read_and_fix("open_source_repos_quality.csv", "open_source"),
                      read_and_fix("open_source_repos_2__2.csv",    "open_source"))
bioinfo      <- rbind(read_and_fix("bioinformatics_repos_quality.csv", "bioinformatics"),
                      read_and_fix("bioinformatics_repos_2__2.csv",    "bioinformatics"))
astro        <- rbind(read_and_fix("astrophysics_repos_quality.csv", "astrophysics"),
                      read_and_fix("astrophysics_repos_2__2.csv",    "astrophysics"))
swdev        <- rbind(read_and_fix("software_engineering_repos_quality.csv", "software_engineering"),
                      read_and_fix("software_engineering_repos_2__2.csv",    "software_engineering"))
imagerec     <- rbind(read_and_fix("image_recognition_repos_quality.csv", "image_recognition"),
                      read_and_fix("image_recognition_repos_2__2.csv",    "image_recognition"))

repo <- rbind(open_source, bioinfo, astro, swdev, imagerec)
repo$watchers             <- as.numeric(repo$watchers)
repo$stars                <- as.numeric(repo$stars)
repo$forks                <- as.numeric(repo$forks)
repo$days_since_last_commit <- as.numeric(repo$days_since_last_commit)
repo$total_issues         <- as.numeric(repo$total_issues)
repo$closed_issues        <- as.numeric(repo$closed_issues)
repo$commit_count         <- as.numeric(repo$commit_count)
repo$open_issues          <- as.numeric(repo$open_issues)

repo <- repo %>%
  filter(days_since_last_commit >= 0) %>%
  filter(commit_count > 4) %>%
  filter(total_issues > 4)


# ── RQI Computation (source of truth: Rmd) ───────────────────────────────────
repo <- repo %>%
  filter(is.na(error) | error == "") %>%
  filter(!(is.na(has_tests) | has_tests == "")) %>%
  mutate(
    backlog_health = 1 / (1 + as.numeric(median_issue_days_open)),
    days_since_last_commit = ifelse(is.na(days_since_last_commit),
                                    max(days_since_last_commit, na.rm = TRUE),
                                    days_since_last_commit),
    popularity    = log1p(stars + forks + watchers),
    ci_present    = ifelse(has_CI    %in% c(TRUE, "True", "Has CI"), 1, 0),
    tests_present = ifelse(has_tests %in% c(TRUE, "True"), 1, 0)
  ) %>%
  mutate(
    recency_z    = scale(-days_since_last_commit),
    issue_z      = scale(backlog_health),
    pop_z        = scale(popularity),
    z_score_mean = rowMeans(cbind(recency_z, issue_z, pop_z), na.rm = TRUE),
    z_winsor     = pmax(pmin(z_score_mean,
                             quantile(z_score_mean, 0.99, na.rm = TRUE)),
                        quantile(z_score_mean, 0.01, na.rm = TRUE)),
    z_normalized_1to5 = scales::rescale(z_winsor, to = c(1, 5))
  ) %>%
  mutate(
    recency_score = scales::rescale(-days_since_last_commit, to = c(0, 100)),
    repo_age_days = as.numeric(repo_age_days),
    commit_count  = as.numeric(commit_count),
    is_active     = as.integer(days_since_last_commit <= 365),
    age_cohort    = factor(case_when(
      repo_age_days <= 730  ~ "0–2 years",
      repo_age_days <= 1825 ~ "2–5 years",
      repo_age_days <= 3650 ~ "5–10 years",
      TRUE                  ~ "10+ years"
    ), levels = c("0–2 years", "2–5 years", "5–10 years", "10+ years"))
  )


# ── Quality groups ────────────────────────────────────────────────────────────
coverage_data <- repo %>%
  mutate(
    quality_group = factor(case_when(
      z_normalized_1to5 >= 4 ~ "High (>=4*)",
      z_normalized_1to5 <= 2 ~ "Low (<=2*)",
      TRUE ~ "Mid (2–4*)"
    ), levels = c("Low (<=2*)", "Mid (2–4*)", "High (>=4*)"))
  )


# ── Radar globals (Rmd rescale_metrics logic) ─────────────────────────────────
global_minmax <- repo %>%
  summarise(
    imin = min(backlog_health, na.rm = TRUE),
    imax = max(backlog_health, na.rm = TRUE),
    pmin = min(popularity, na.rm = TRUE),
    pmax = max(popularity, na.rm = TRUE)
  )
global_max <- c(100, 100, 100)
global_min <- c(0,   0,   0)

rescale_metrics <- function(df) {
  df %>% mutate(
    backlog_health = scales::rescale(backlog_health, to = c(0,100),
                                     from = c(global_minmax$imin, global_minmax$imax)),
    popularity     = scales::rescale(popularity,     to = c(0,100),
                                     from = c(global_minmax$pmin, global_minmax$pmax))
  )
}

top100_avg <- repo %>%
  group_by(type) %>%
  arrange(desc(z_normalized_1to5)) %>%
  slice_head(n = 100) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

bottom100_avg <- repo %>%
  group_by(type) %>%
  arrange(z_normalized_1to5) %>%
  slice_head(n = 100) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

median100_avg <- repo %>%
  group_by(type) %>%
  arrange(z_normalized_1to5) %>%
  slice(round(n()/2 - 49) : round(n()/2 + 50)) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

# ── Radar: by-domain & top20/bottom20-by-domain (Rmd Figs S10 / S11) ──────────
all_by_type_radar <- repo %>%
  group_by(type) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

top20_by_type_radar <- repo %>%
  group_by(type) %>%
  arrange(desc(z_normalized_1to5)) %>%
  slice_head(n = 20) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

bottom20_by_type_radar <- repo %>%
  group_by(type) %>%
  arrange(z_normalized_1to5) %>%
  slice_head(n = 20) %>%
  summarise(recency = mean(recency_score, na.rm=TRUE),
            backlog_health = mean(backlog_health, na.rm=TRUE),
            popularity = mean(popularity, na.rm=TRUE), .groups="drop") %>%
  rescale_metrics()

build_single_radar_df <- function(row_df) {
  df <- rbind(global_max, global_min,
              as.numeric(row_df[1, c("recency", "backlog_health", "popularity")]))
  df <- as.data.frame(df)
  colnames(df) <- c("Recent\nActivity", "Backlog\nHealth", "Popularity")
  df
}

# ── Generic helper: pairwise Wilcoxon significance brackets ──────────────────
build_sig_brackets <- function(data, group_col, value_col) {
  groups    <- levels(droplevels(factor(data[[group_col]])))
  group_pos <- setNames(seq_along(groups), groups)
  out <- list()
  if (length(groups) >= 2) {
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        a <- groups[i]; b <- groups[j]
        sub <- data[data[[group_col]] %in% c(a, b), ]
        sub[[group_col]] <- droplevels(factor(sub[[group_col]]))
        if (nlevels(sub[[group_col]]) < 2) next
        form <- stats::as.formula(paste(value_col, "~", group_col))
        p <- tryCatch(wilcox.test(form, data = sub)$p.value, error = function(e) NA)
        if (is.na(p)) next
        lbl <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
        if (lbl != "ns") {
          out[[length(out) + 1]] <- data.frame(g1 = a, g2 = b,
                                               x = group_pos[[a]], xend = group_pos[[b]],
                                               label = lbl, stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(out) == 0) {
    return(data.frame(g1 = character(), g2 = character(),
                      x = numeric(), xend = numeric(), label = character()))
  }
  do.call(rbind, out)
}

# ── Generic helper: rescale a z-score vector to 1–5 (Rmd scale_1to5) ──────────
scale_1to5 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  1 + 4 * (x - rng[1]) / (rng[2] - rng[1])
}


require_two_groups <- function(g, label = "group") {
  validate(need(length(g) > 0,
                paste0("No rows matched for this ", label, " comparison — check that the ",
                       "filter (e.g. domain/type spelling) matches the data exactly.")))
  g <- droplevels(factor(g))
  validate(need(nlevels(g) == 2,
                paste0("Only ", nlevels(g), " ", label, " categor",
                       if (nlevels(g) == 1) "y is" else "ies are",
                       " present in this subset, so no comparison can be made.")))
  g
}


# ── ANOVA / Tukey ─────────────────────────────────────────────────────────────
anova_result <- aov(z_normalized_1to5 ~ type, data = repo)
Tukey_result <- TukeyHSD(anova_result)

tukey_df <- as.data.frame(Tukey_result$type)
tukey_df$pair <- rownames(tukey_df)

tukey_sig <- tukey_df %>%
  filter(`p adj` < 0.05) %>%
  separate(pair, into = c("group1", "group2"), sep = "-")

repo_summary <- repo %>%
  group_by(type) %>%
  summarise(mean = mean(z_normalized_1to5, na.rm=TRUE),
            sd   = sd(z_normalized_1to5, na.rm=TRUE),
            n    = n(),
            se   = sd / sqrt(n),
            ci_lower = mean - qt(0.975, df=n-1)*se,
            ci_upper = mean + qt(0.975, df=n-1)*se)


# ── Within-type CI/test t-tests ───────────────────────────────────────────────
within_type_tests <- repo %>%
  group_by(type) %>%
  summarise(
    ci_p = ifelse(length(unique(ci_present)) > 1 &
                    sum(ci_present==0,na.rm=TRUE)>=2 &
                    sum(ci_present==1,na.rm=TRUE)>=2,
                  t.test(z_normalized_1to5 ~ ci_present)$p.value, NA),
    tests_p = ifelse(length(unique(tests_present)) > 1 &
                       sum(tests_present==0,na.rm=TRUE)>=2 &
                       sum(tests_present==1,na.rm=TRUE)>=2,
                     t.test(z_normalized_1to5 ~ tests_present)$p.value, NA),
    mean_ci_yes   = mean(z_normalized_1to5[ci_present==1],    na.rm=TRUE),
    mean_ci_no    = mean(z_normalized_1to5[ci_present==0],    na.rm=TRUE),
    mean_tests_yes = mean(z_normalized_1to5[tests_present==1], na.rm=TRUE),
    mean_tests_no  = mean(z_normalized_1to5[tests_present==0], na.rm=TRUE),
    .groups = "drop"
  )


# ── Language processing ───────────────────────────────────────────────────────
lang_long <- repo %>%
  select(repo, type, z_normalized_1to5, tests_present, ci_present, languages) %>%
  filter(!is.na(languages), languages != "") %>%
  mutate(languages = str_replace_all(languages, "\\s+", " ")) %>%
  separate_rows(languages, sep = ",\\s*") %>%
  mutate(
    language = str_trim(str_extract(languages, "^[^\\(]+")),
    bytes    = as.numeric(str_extract(languages, "(?<=\\()[0-9]+"))
  ) %>%
  filter(!is.na(language), language != "", !is.na(bytes))

lang_primary <- lang_long %>%
  group_by(repo) %>%
  slice_max(order_by = bytes, n = 1, with_ties = FALSE) %>%
  ungroup()

min_n <- 20
lang_summary_filt <- lang_primary %>%
  group_by(language) %>%
  summarise(n = n(), mean_rating = mean(z_normalized_1to5, na.rm=TRUE), .groups="drop") %>%
  filter(n >= min_n)

# Bioinformatics-only language summary (Rmd Figs S15–S17 use a lower n threshold)
min_n_bio <- 10
bio_lang_summary_filt <- lang_primary %>%
  filter(type == "bioinformatics") %>%
  group_by(language) %>%
  summarise(n = n(), mean_rating = mean(z_normalized_1to5, na.rm=TRUE), .groups="drop") %>%
  filter(n >= min_n_bio)


# ── Cohort summaries ──────────────────────────────────────────────────────────
cohort_ci_summary <- repo %>%
  group_by(age_cohort, ci_present) %>%
  summarise(mean_z = mean(z_normalized_1to5,na.rm=TRUE),
            se_z   = sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),
            n=n(), .groups="drop") %>%
  mutate(ci_label = ifelse(ci_present==1,"Has CI","No CI"))

cohort_test_summary <- repo %>%
  group_by(age_cohort, tests_present) %>%
  summarise(mean_z = mean(z_normalized_1to5,na.rm=TRUE),
            se_z   = sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),
            n=n(), .groups="drop") %>%
  mutate(test_label = ifelse(tests_present==1,"Has Tests","No Tests"))


# ── Summary tables (Rmd Table 1 / Table S1) ───────────────────────────────────
summary_tbl_type <- repo %>%
  group_by(type) %>%
  summarise(mean_rqi        = mean(z_normalized_1to5, na.rm = TRUE),
            mean_recency    = mean(recency_score,     na.rm = TRUE),
            mean_backlog    = mean(backlog_health,    na.rm = TRUE),
            mean_popularity = mean(popularity,        na.rm = TRUE),
            n               = n(),
            .groups = "drop") %>%
  arrange(desc(mean_rqi))

top_languages_tbl <- lang_long %>% count(language, sort = TRUE) %>% head(15)


# ── Shared colour palette ─────────────────────────────────────────────────────
type_colors <- c(
  "bioinformatics"     = "#4C72B0",
  "astrophysics"       = "#DD8452",
  "image_recognition"  = "#937860",
  "open_source"        = "#64B5CD",
  "software_engineering" = "#9467BD"
)

unique_type <- unique(repo$type)


# ── Citation data (loaded if file exists; mirrors Rmd's publication logic) ────
cit_file <- "repos_with_citations_gemma_v2.csv"
has_citations <- file.exists(cit_file)
if (has_citations) {
  citations_raw <- read_csv(cit_file, show_col_types = FALSE) %>%
    distinct(repo, .keep_all = TRUE) %>%
    mutate(
      
      doi_is_paper = !is.na(doi) &
        str_detect(doi, "^10\\.") &
        !str_detect(doi, "^10\\.(5281|6084|17605|32614|5438)/"),
      published = (citation_verified %in% TRUE) | doi_is_paper,
      pub_label = factor(if_else(published, "Published", "Not Published"),
                         levels = c("Not Published", "Published"))
    )
}

# ── Funding data (loaded if file exists; mirrors Rmd's currency conversion) ───
nih_file   <- "oa_publication_summary_gemma.csv"
has_funding <- file.exists(nih_file)
if (has_funding) {
  nih <- read_csv(nih_file, show_col_types = FALSE)
  
  usd_per_unit <- c(
    USD = 1.00, GBP = 1.34, EUR = 1.16, CHF = 1.25, AUD = 0.67,
    CAD = 0.72, CNY = 0.14, HKD = 0.128, SEK = 0.108, JPY = 0.0063, CLP = 0.00105
  )
  
  convert_cost_string <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    total <- 0
    matched_any <- FALSE
    for (part in strsplit(x, ";\\s*")[[1]]) {
      kv <- strsplit(part, ":")[[1]]
      if (length(kv) != 2) next
      cur  <- trimws(kv[1])
      amt  <- suppressWarnings(as.numeric(kv[2]))
      rate <- usd_per_unit[cur]
      if (!is.na(amt) && !is.na(rate)) {
        total <- total + amt * rate
        matched_any <- TRUE
      }
    }
    if (matched_any) total else NA_real_
  }
  
  nih <- nih %>%
    mutate(
      developer_cost_usd = sapply(developer_cost_by_currency, convert_cost_string),
      total_cost_usd     = sapply(total_cost_by_currency,     convert_cost_string)
    )
  
  funding_df <- nih %>%
    left_join(repo %>% select(repo, z_normalized_1to5, recency_z, issue_z, pop_z, type),
              by = "repo") %>%
    mutate(
      funding_group = factor(
        ifelse(is.na(funding_sources) | funding_sources == "",
               "Not Grant Funded", "Grant Funded"),
        levels = c("Not Grant Funded", "Grant Funded")
      )
    )
  
  funded_sources_list <- funding_df %>%
    filter(!is.na(funding_sources) & funding_sources != "") %>%
    pull(funding_sources) %>%
    strsplit("; ") %>% unlist() %>% trimws() %>% unique() %>% sort()
}


# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel(
    tags$div(
      style = "line-height: 1.35; padding: 6px 0;",
      tags$h2(
        "Continuous Integration and Software Quality in Scientific Software: A Large-Scale Empirical Analysis of GitHub Repositories",
        style = "font-weight: 700; font-size: 26px; margin-bottom: 6px;"
      ),
      tags$div(
        style = "font-size: 17px; color: #555;",
        tags$a(
          href = "https://doi.org/10.XXXX/XXXXXXX", 
          target = "_blank",
          style = "color: #555; text-decoration: underline;",
          "DOI: 10.XXXX/XXXXXXX" 
        ),
        tags$span(" \u2003|\u2003 ", style = "color:#aaa;"),
        tags$a(
          href = "https://paper-link.example.com",  
          target = "_blank",
          icon("file-text"), "Read the paper"
        ),
        tags$span(" \u2003|\u2003 ", style = "color:#aaa;"),
        tags$a(
          href = "https://jakobilab.org",
          target = "_blank",
          icon("globe"), "jakobilab.org"
        )
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 # Repo-type filter (most tabs)
                 conditionalPanel(
                   condition = paste0("['distribution_tab','coverage_tab','impact_tab',",
                                      "'impact_by_type_tab','radar_tab','cohort_tab',",
                                      "'lang_tab'].indexOf(input.main_tabs) >= 0"),
                   selectInput("repo_type", "Repository Type",
                               choices = c("All", unique_type), selected = "All")
                 ),
                 
                 # Metric filter (coverage tab)
                 conditionalPanel(
                   condition = "input.main_tabs == 'coverage_tab'",
                   selectInput("test_type", "Metric",
                               choices = c("All", "ci_coverage", "test_coverage"), selected = "All")
                 ),
                 
                 # Pair selector (Tukey comparison tab)
                 conditionalPanel(
                   condition = "input.main_tabs == 'comparison_tab'",
                   selectInput("repo_1", "First Repo Type",  choices = unique_type),
                   selectInput("repo_2", "Second Repo Type", choices = unique_type, selected = unique_type[2])
                 ),
                 
                 # Radar view selector
                 conditionalPanel(
                   condition = "input.main_tabs == 'radar_tab'",
                   selectInput("radar_view", "Radar View",
                               choices = c("Top / Mid / Bottom 100"        = "tmb",
                                           "Average Profile by Domain"      = "by_domain",
                                           "Top 20 vs Bottom 20 by Domain"  = "top_bottom_20"),
                               selected = "tmb"),
                   helpText("\"Repository Type\" above only applies to the Top/Mid/Bottom view.")
                 ),
                 
                 # Language filters
                 conditionalPanel(
                   condition = "input.main_tabs == 'lang_tab'",
                   selectInput("lang_scope", "Domain Scope",
                               choices = c("All Domains" = "all", "Bioinformatics Only" = "bio"),
                               selected = "all"),
                   selectInput("lang_split", "Split by",
                               choices = c("None" = "none", "Tests" = "tests", "CI" = "ci", "QA Bucket" = "qa"),
                               selected = "tests")
                 ),
                 
                 # Citation sub-tab
                 conditionalPanel(
                   condition = "input.main_tabs == 'citation_tab'",
                   selectInput("cit_plot", "Citation Plot",
                               choices = c("Mean vs RQI (loess)"                    = "loess",
                                           "Median vs RQI (loess)"                  = "median_loess",
                                           "Violin by RQI bin"                      = "violin",
                                           "Mean bar by RQI bin"                    = "mean_bar",
                                           "Median bar (Low/High)"                  = "bar_lh",
                                           "Median bar (RQI coarse bins)"           = "bar_bins",
                                           "Median bar (RQI fine bins)"             = "bar_bins_fine",
                                           "Median bar, published only (Fig 6)"     = "fig6",
                                           "Citations vs CI"                        = "bar_ci",
                                           "Median RQI by publication status"       = "pub_bar",
                                           "RQI distribution by publication status" = "pub_violin",
                                           "RQI components by publication (Fig 5A)" = "fig5a",
                                           "Citations vs RQI + Age"                 = "age_scatter",
                                           "Scatter + LM (published)"               = "scatter_lm"),
                               selected = "loess")
                 ),
                 
                 # Funding sub-tab
                 conditionalPanel(
                   condition = "input.main_tabs == 'funding_tab'",
                   selectInput("fund_plot", "Funding Plot",
                               choices = c("Median RQI by Funding (Bioinformatics)"              = "bar_bio",
                                           "RQI Components by Funding, Fig 5B (Bioinformatics)"  = "fig5b_bio",
                                           "RQI vs Developer Award (Bioinformatics)"             = "lm_dev_bio",
                                           "RQI vs Organization Award (Bioinformatics)"          = "lm_org_bio",
                                           "RQI Components by Funding (All Domains)"             = "fig_all",
                                           "RQI vs Developer Award (All Domains)"                = "lm_dev_all",
                                           "RQI vs Organization Award (All Domains)"             = "lm_org_all"),
                               selected = "bar_bio")
                 )
    ),
    
    mainPanel(
      tabsetPanel(id = "main_tabs",
                  
                  # 1 ── Coverage
                  tabPanel("CI & Test Coverage", value = "coverage_tab",
                           h4("CI and Testing Prevalence by Quality Tier"),
                           plotOutput("test_plot", height = "400px"),
                           hr(),
                           h4("Prevalence Regression Across RQI"),
                           plotOutput("prev_reg_plot", height = "400px")
                  ),
                  
                  # 2 ── Distribution
                  tabPanel("Distribution", value = "distribution_tab",
                           h4("Quality Distribution by Type (stacked)"),
                           plotOutput("dist_plot", height = "380px"),
                           hr(),
                           h4("Overall Distribution + Q-Q"),
                           plotOutput("dist_qq_plot", height = "350px"),
                           hr(),
                           h4("Density by Type (faceted)"),
                           plotOutput("dist_facet_plot", height = "420px")
                  ),
                  
                  # 3 ── ANOVA / Tukey overall
                  tabPanel("ANOVA / Tukey", value = "anova_tab",
                           h4("Bioinformatics vs All — Tukey Significant Pairs"),
                           plotOutput("bio_tukey_plot", height = "480px")
                  ),
                  
                  # 4 ── Tukey pairwise comparison
                  tabPanel("Pairwise Comparison", value = "comparison_tab",
                           h4("Tukey HSD Pairwise Comparison"),
                           plotOutput("tukey_plot", height = "420px")
                  ),
                  
                  # 5 ── CI & Test impact (overall)
                  tabPanel("CI & Test Impact", value = "impact_tab",
                           fluidRow(
                             column(6, h4("CI Impact"), plotOutput("ci_impact_plot", height = "380px")),
                             column(6, h4("Test Impact"), plotOutput("test_impact_plot", height = "380px"))
                           )
                  ),
                  
                  # 6 ── CI & Test impact by type (faceted)
                  tabPanel("Impact by Type", value = "impact_by_type_tab",
                           h4("CI Presence Effect per Repository Type"),
                           plotOutput("ci_by_type_plot", height = "500px"),
                           hr(),
                           h4("Testing Presence Effect per Repository Type"),
                           plotOutput("test_by_type_plot", height = "500px")
                  ),
                  
                  # 7 ── Radar
                  tabPanel("Radar", value = "radar_tab",
                           h4("Repository Quality Radar Profiles"),
                           uiOutput("radar_plot_ui")
                  ),
                  
                  # 8 ── Age cohort
                  tabPanel("Age Cohorts", value = "cohort_tab",
                           h4("Mean RQI by Age Cohort — CI"),
                           plotOutput("cohort_ci_plot", height = "360px"),
                           hr(),
                           h4("Mean RQI by Age Cohort — Testing"),
                           plotOutput("cohort_test_plot", height = "360px"),
                           hr(),
                           h4("Bioinformatics CI & Test Coverage by Cohort"),
                           plotOutput("bio_cohort_plot", height = "360px"),
                           hr(),
                           h4("Combined View (Journal Figure Style)"),
                           plotOutput("cohort_combined_plot", height = "420px")
                  ),
                  
                  # 9 ── Language
                  tabPanel("Language", value = "lang_tab",
                           h4("Mean RQI by Primary Language"),
                           plotOutput("lang_plot", height = "520px")
                  ),
                  
                  # 10 ── Activity (survival)
                  tabPanel("Activity / Survival", value = "survival_tab",
                           h4("Odds Ratios: Predictors of Repo Activity"),
                           plotOutput("survival_plot", height = "480px"),
                           hr(),
                           h4("Activity Rate by CI/Test Presence and Domain"),
                           plotOutput("activity_plot", height = "480px")
                  ),
                  
                  # 11 ── Citations (optional)
                  tabPanel("Citations", value = "citation_tab",
                           uiOutput("citation_ui")
                  ),
                  
                  # 12 ── Funding (optional)
                  tabPanel("Funding", value = "funding_tab",
                           uiOutput("funding_ui")
                  ),
                  
                  # 13 ── Summary tables
                  tabPanel("Summary Tables", value = "summary_tab",
                           h4("Table 1 — Mean Scores by Repository Type"),
                           tableOutput("table1_summary"),
                           hr(),
                           h4("Table S1 — Most Frequent Languages"),
                           tableOutput("tableS1_languages")
                  )
      )
    )
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output) {
  
  # ── Reactive filtered repo ─────────────────────────────────────────────────
  filtered_repo <- reactive({
    if (input$repo_type == "All") repo else filter(repo, type == input$repo_type)
  })
  
  # ── 1a. Coverage bar ────────────────────────────────────────────────────────
  output$test_plot <- renderPlot({
    if (input$repo_type == "All") {
      fs <- coverage_data %>%
        group_by(quality_group) %>%
        summarise(ci_coverage   = mean(ci_present==1,    na.rm=TRUE),
                  test_coverage = mean(tests_present==1, na.rm=TRUE),
                  .groups="drop") %>%
        pivot_longer(c(ci_coverage,test_coverage), names_to="metric", values_to="coverage")
    } else {
      fs <- coverage_data %>%
        filter(type == input$repo_type) %>%
        group_by(quality_group) %>%
        summarise(ci_coverage   = mean(ci_present==1,    na.rm=TRUE),
                  test_coverage = mean(tests_present==1, na.rm=TRUE),
                  .groups="drop") %>%
        pivot_longer(c(ci_coverage,test_coverage), names_to="metric", values_to="coverage")
    }
    if (input$test_type != "All") fs <- filter(fs, metric == input$test_type)
    
    ggplot(fs, aes(x=quality_group, y=coverage*100, fill=metric)) +
      geom_col(position=position_dodge(0.6), width=0.6, color="white") +
      geom_text(aes(label=paste0(round(coverage*100,1),"%")),
                position=position_dodge(0.6), vjust=-0.5, size=3.5) +
      scale_fill_manual(values=c("ci_coverage"="#1F78B4","test_coverage"="#33A02C"),
                        labels=c("Continuous Integration","Testing")) +
      labs(title="CI and Testing Prevalence by Quality Tier",
           x="Quality Group by RQI", y="Prevalence (%)", fill="Metric") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(plot.title=element_text(face="bold",hjust=0.5))
  })
  
  # ── 1b. Prevalence regression ───────────────────────────────────────────────
  output$prev_reg_plot <- renderPlot({
    df <- filtered_repo()
    prev_reg <- df %>%
      mutate(rqi_bin = round(z_normalized_1to5/0.15)*0.15) %>%
      group_by(rqi_bin) %>%
      summarise(ci_prev   = mean(ci_present==1,    na.rm=TRUE),
                test_prev = mean(tests_present==1, na.rm=TRUE),
                n=n(), .groups="drop") %>%
      filter(n >= 5)
    
    lm_ci   <- lm(ci_prev   ~ rqi_bin, data=prev_reg)
    lm_test <- lm(test_prev ~ rqi_bin, data=prev_reg)
    
    subtitle_str <- sprintf(
      "CI: slope=%.2f, R²=%.2f, p%s   |   Testing: slope=%.2f, R²=%.2f, p%s",
      coef(lm_ci)[2],   summary(lm_ci)$r.squared,
      ifelse(summary(lm_ci)$coefficients[2,4]<0.001,"<0.001",sprintf("=%.3f",summary(lm_ci)$coefficients[2,4])),
      coef(lm_test)[2], summary(lm_test)$r.squared,
      ifelse(summary(lm_test)$coefficients[2,4]<0.001,"<0.001",sprintf("=%.3f",summary(lm_test)$coefficients[2,4]))
    )
    
    prev_reg %>%
      pivot_longer(c(ci_prev,test_prev), names_to="metric", values_to="prevalence") %>%
      ggplot(aes(x=rqi_bin, y=prevalence*100, color=metric, fill=metric)) +
      geom_point(aes(size=n), alpha=0.5) +
      geom_smooth(method="lm", se=TRUE, alpha=0.15, linewidth=1.2) +
      scale_color_manual(values=c("ci_prev"="#1F78B4","test_prev"="#33A02C"),
                         labels=c("Continuous Integration","Testing")) +
      scale_fill_manual(values=c("ci_prev"="#1F78B4","test_prev"="#33A02C"),
                        labels=c("Continuous Integration","Testing")) +
      scale_size_continuous(range=c(1.5,5), guide="none") +
      labs(title="CI and Testing Prevalence Across Repository Quality",
           subtitle=subtitle_str, x="RQI (1–5)", y="Prevalence (%)",
           color="Metric", fill="Metric") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(plot.title=element_text(face="bold",hjust=0.5),
            plot.subtitle=element_text(hjust=0.5,color="gray40",size=16),
            legend.position="bottom")
  })
  
  # ── 2a. Distribution stacked ────────────────────────────────────────────────
  output$dist_plot <- renderPlot({
    df <- filtered_repo()
    mean_df <- df %>% group_by(type) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),.groups="drop")
    
    ggplot(df, aes(x=z_normalized_1to5, fill=type)) +
      geom_histogram(aes(y=after_stat(count/sum(count))), position="stack",
                     binwidth=0.15, color="white", alpha=0.7) +
      geom_density(aes(y=after_stat(..scaled..), color=type), size=1.1, alpha=0.9) +
      geom_vline(data=mean_df, aes(xintercept=mean_z, color=type),
                 linetype="dashed", linewidth=1) +
      scale_fill_manual(values=type_colors)  +
      scale_color_manual(values=type_colors) +
      labs(title="Comparative Quality Distribution Across Repository Types",
           x="Normalized (1–5)", y="Density / Proportion",
           fill="Repository Type", color="Repository Type") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(legend.position="top", plot.title=element_text(face="bold",hjust=0.5))
  })
  
  # ── 2b. Overall dist + Q-Q ──────────────────────────────────────────────────
  output$dist_qq_plot <- renderPlot({
    df     <- filtered_repo()
    z_all  <- df$z_normalized_1to5[!is.na(df$z_normalized_1to5)]
    z_samp <- if (length(z_all) > 10000) sample(z_all,10000) else z_all
    
    p_dist <- ggplot(df, aes(x=z_normalized_1to5)) +
      geom_histogram(aes(y=after_stat(density)), bins=40,
                     fill="#56B4E9", color="white", alpha=0.8) +
      geom_density(color="red", size=1.1) +
      geom_vline(aes(xintercept=mean(z_normalized_1to5,na.rm=TRUE)),
                 color="darkblue", linetype="dashed", linewidth=1) +
      labs(title="Overall Distribution of Normalized RQI",
           x="RQI (1–5)", y="Density") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
    
    p_qq <- ggplot(data.frame(z=z_samp), aes(sample=z)) +
      stat_qq(color="#1F78B4") +
      stat_qq_line(color="red", linewidth=1) +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
    
    p_dist + p_qq
  })
  
  # ── 2c. Density faceted ─────────────────────────────────────────────────────
  output$dist_facet_plot <- renderPlot({
    df <- filtered_repo()
    ggplot(df, aes(x=z_normalized_1to5, fill=type)) +
      geom_histogram(aes(y=after_stat(density)), bins=30, color="white", alpha=0.7) +
      geom_density(aes(color=type), size=1.1) +
      facet_wrap(~type, ncol=3) +
      scale_fill_manual(values=type_colors) +
      scale_color_manual(values=type_colors) +
      labs(title="Distribution and Density by Repository Type",
           x="RQI (1–5)", y="Density") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(legend.position="none",
            plot.title=element_text(face="bold",hjust=0.5))
  })
  
  # ── 3. Bio Tukey ────────────────────────────────────────────────────────────
  output$bio_tukey_plot <- renderPlot({
    bio_comp <- as.data.frame(Tukey_result$type) %>%
      tibble::rownames_to_column("comparison") %>%
      separate(comparison, into=c("group1","group2"), sep="-") %>%
      filter(group1=="bioinformatics" | group2=="bioinformatics") %>%
      mutate(other = ifelse(group1=="bioinformatics",group2,group1),
             sig = case_when(`p adj`<0.001~"***",`p adj`<0.01~"**",
                             `p adj`<0.05~"*",TRUE~"ns"))
    
    all_types   <- c("bioinformatics", bio_comp$other)
    plot_data   <- repo_summary %>% filter(type %in% all_types) %>%
      mutate(type=factor(type, levels=c("bioinformatics",sort(bio_comp$other))))
    type_levels <- levels(plot_data$type)
    bio_x       <- which(type_levels=="bioinformatics")
    bio_ci_up   <- repo_summary %>% filter(type=="bioinformatics") %>% pull(ci_upper)
    col_pal     <- setNames(c("#1F78B4",
                              RColorBrewer::brewer.pal(length(type_levels)-1,"Set2")),
                            type_levels)
    
    sig_bars <- bio_comp %>%
      filter(sig!="ns") %>%
      mutate(other_x=match(other,type_levels), bio_x=bio_x) %>%
      left_join(repo_summary %>% select(type,ci_upper), by=c("other"="type")) %>%
      rename(ci_upper_other=ci_upper) %>%
      mutate(y=pmax(ci_upper_other,bio_ci_up)+0.05+row_number()*0.12)
    
    ggplot(plot_data, aes(x=type, y=mean, fill=type)) +
      geom_col(alpha=0.85, width=0.65) +
      geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper), width=0.2) +
      geom_segment(data=sig_bars,
                   aes(x=bio_x,xend=other_x,y=y,yend=y),
                   inherit.aes=FALSE, linewidth=0.6) +
      geom_text(data=sig_bars,
                aes(x=(bio_x+other_x)/2,y=y+0.05,label=sig),
                inherit.aes=FALSE, size=5, fontface="bold") +
      scale_fill_manual(values=col_pal) +
      labs(x=NULL, y="Mean RQI", title="Bioinformatics vs All Other Categories") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=30,hjust=1,size=18),
            plot.title=element_text(face="bold",hjust=0.5))
  })
  
  # ── 4. Pairwise Tukey ───────────────────────────────────────────────────────
  output$tukey_plot <- renderPlot({
    req(input$repo_1, input$repo_2)
    validate(need(input$repo_1 != input$repo_2, "Select two different types."))
    
    pair_data <- repo_summary %>% filter(type %in% c(input$repo_1,input$repo_2))
    sig_check <- tukey_sig %>%
      filter((group1==input$repo_1 & group2==input$repo_2) |
               (group1==input$repo_2 & group2==input$repo_1))
    sig_label <- if (nrow(sig_check)>=1) paste0("p = ",signif(as.numeric(sig_check$`p adj`[1]),3)) else "p = ns"
    
    y_bar  <- max(pair_data$ci_upper,na.rm=TRUE)+0.05
    y_text <- y_bar+0.03
    
    ggplot(pair_data, aes(x=type, y=mean, fill=type)) +
      geom_col(alpha=0.8, width=0.6) +
      geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper), width=0.15) +
      geom_segment(aes(x=1,xend=2,y=y_bar,yend=y_bar), linewidth=1) +
      geom_text(aes(x=1.5,y=y_text,label=sig_label), size=5, fontface="bold") +
      labs(title=paste(input$repo_1,"vs",input$repo_2),
           subtitle="Tukey HSD Comparison", x=NULL, y="Mean RQI (±95% CI)") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=30,hjust=1,size=18),
            plot.title=element_text(face="bold",hjust=0.5))
  })
  
  # ── 5. CI / test impact overall ─────────────────────────────────────────────
  make_impact_plot <- function(df, var, colors, xlab, title_prefix) {
    label_vals <- if (var=="ci_present") c("No CI","Has CI") else c("No Tests","Has Tests")
    smry <- df %>%
      mutate(label=ifelse(.data[[var]]==1, label_vals[2], label_vals[1])) %>%
      group_by(label) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                sd_z=sd(z_normalized_1to5,na.rm=TRUE), n=n(),
                se_z=sd_z/sqrt(n), .groups="drop")
    
    p_val <- if (nrow(smry)==2) t.test(df$z_normalized_1to5 ~ df[[var]])$p.value else NA
    sig   <- if (!is.na(p_val)) paste0("p = ",signif(p_val,3)) else ""
    
    ggplot(smry, aes(x=label,y=mean_z,fill=label)) +
      geom_col(width=0.6,color="white") +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),width=0.2) +
      geom_segment(aes(x=1,xend=2,y=max(mean_z+se_z)+0.05,yend=max(mean_z+se_z)+0.05)) +
      geom_text(aes(x=1.5,y=max(mean_z+se_z)+0.1,label=sig),
                size=5, fontface="bold") +
      scale_fill_manual(values=colors) +
      labs(title=paste(title_prefix,"Impact on Quality"), x=xlab, y="Mean RQI (1–5)") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(legend.position="none")
  }
  
  output$ci_impact_plot   <- renderPlot(make_impact_plot(filtered_repo(),"ci_present",
                                                         c("No CI"="#A6CEE3","Has CI"="#1F78B4"),"CI","CI"))
  output$test_impact_plot <- renderPlot(make_impact_plot(filtered_repo(),"tests_present",
                                                         c("No Tests"="#B2DF8A","Has Tests"="#33A02C"),"Testing","Testing"))
  
  # ── 6. CI / test impact by type (faceted) ───────────────────────────────────
  output$ci_by_type_plot <- renderPlot({
    ci_group <- repo %>%
      group_by(type,ci_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),
                n=n(),.groups="drop") %>%
      left_join(within_type_tests %>% select(type,ci_p), by="type") %>%
      mutate(sig_label=case_when(ci_p<0.001~"***",ci_p<0.01~"**",ci_p<0.05~"*",TRUE~"ns"))
    
    ggplot(ci_group, aes(x=factor(ci_present),y=mean_z,fill=factor(ci_present))) +
      geom_col(position="dodge",width=0.7,color="white") +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),width=0.2,color="black") +
      facet_wrap(~type, scales="free_y") +
      geom_segment(data=ci_group %>% distinct(type,.keep_all=TRUE),
                   aes(x=1,xend=2,y=max(mean_z+se_z,na.rm=TRUE)+0.05,
                       yend=max(mean_z+se_z,na.rm=TRUE)+0.05),color="black") +
      geom_text(data=ci_group %>% distinct(type,.keep_all=TRUE),
                aes(x=1.5,y=max(mean_z+se_z,na.rm=TRUE)+0.1,label=sig_label),
                size=5,fontface="bold") +
      scale_fill_manual(values=c("0"="#A6CEE3","1"="#1F78B4")) +
      scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1)) +
      labs(title="CI Presence Effect per Repository Type",
           x="CI (0=No, 1=Yes)", y="RQI (1–5)") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(strip.text=element_text(face="bold"), legend.position="none")
  })
  
  output$test_by_type_plot <- renderPlot({
    test_group <- repo %>%
      group_by(type,tests_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),
                n=n(),.groups="drop") %>%
      left_join(within_type_tests %>% select(type,tests_p), by="type") %>%
      mutate(sig_label=case_when(tests_p<0.001~"***",tests_p<0.01~"**",tests_p<0.05~"*",TRUE~"ns"))
    
    ggplot(test_group, aes(x=factor(tests_present),y=mean_z,fill=factor(tests_present))) +
      geom_col(position="dodge",width=0.7,color="white") +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),width=0.2,color="black") +
      facet_wrap(~type, scales="free_y") +
      geom_segment(data=test_group %>% distinct(type,.keep_all=TRUE),
                   aes(x=1,xend=2,y=max(mean_z+se_z,na.rm=TRUE)+0.05,
                       yend=max(mean_z+se_z,na.rm=TRUE)+0.05),color="black") +
      geom_text(data=test_group %>% distinct(type,.keep_all=TRUE),
                aes(x=1.5,y=max(mean_z+se_z,na.rm=TRUE)+0.1,label=sig_label),
                size=5,fontface="bold") +
      scale_fill_manual(values=c("0"="#B2DF8A","1"="#33A02C")) +
      scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1)) +
      labs(title="Testing Presence Effect per Repository Type",
           x="Tests (0=No, 1=Yes)", y="RQI (1–5)") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(strip.text=element_text(face="bold"), legend.position="none")
  })
  
  # ── 7. Radar ────────────────────────────────────────────────────────────────
  output$radar_plot_ui <- renderUI({
    h <- switch(input$radar_view,
                "tmb"           = "580px",
                "by_domain"     = "820px",
                "top_bottom_20" = "1150px",
                "580px")
    plotOutput("radar_plot", height = h)
  })
  
  output$radar_plot <- renderPlot({
    if (input$radar_view == "tmb") {
      
      if (input$repo_type == "All") {
        top_df    <- top100_avg    %>% summarise(across(c(recency,backlog_health,popularity),mean))
        mid_df    <- median100_avg %>% summarise(across(c(recency,backlog_health,popularity),mean))
        bottom_df <- bottom100_avg %>% summarise(across(c(recency,backlog_health,popularity),mean))
      } else {
        top_df    <- top100_avg    %>% filter(type==input$repo_type)
        mid_df    <- median100_avg %>% filter(type==input$repo_type)
        bottom_df <- bottom100_avg %>% filter(type==input$repo_type)
      }
      req(nrow(top_df)==1, nrow(mid_df)==1, nrow(bottom_df)==1)
      
      radar_data <- as.data.frame(rbind(
        global_max, global_min,
        as.numeric(top_df[,c("recency","backlog_health","popularity")]),
        as.numeric(mid_df[,c("recency","backlog_health","popularity")]),
        as.numeric(bottom_df[,c("recency","backlog_health","popularity")])
      ))
      colnames(radar_data) <- c("Recent Activity\n(outer=more recent)",
                                "Backlog Health\n(outer=faster)",
                                "Popularity\n(outer=more popular)")
      rownames(radar_data) <- c("Max","Min","Top 100","Mid 100","Bottom 100")
      
      par(mar=c(2,2,4,2))
      fmsb::radarchart(radar_data, axistype=1,
                       pcol=c("#1B9E77","#7570B3","#D95F02"),
                       pfcol=scales::alpha(c("#1B9E77","#7570B3","#D95F02"),0.35),
                       plwd=3, cglcol="grey", cglty=1, axislabcol="grey40",
                       caxislabels=c("0","25","50","75","100"), vlcex=1.6,
                       title=paste("Top vs Mid vs Bottom —", input$repo_type))
      legend("topright", legend=c("Top 100","Mid 100","Bottom 100"),
             col=c("#1B9E77","#7570B3","#D95F02"), lty=1, lwd=3, bty="n")
      
    } else if (input$radar_view == "by_domain") {
      
      categories   <- sort(unique(repo$type))
      n_cats_local <- length(categories)
      ncol_use <- min(3, n_cats_local)
      nrow_use <- ceiling(n_cats_local / ncol_use)
      
      par(mfrow = c(nrow_use, ncol_use), mar = c(1,2,4,2), oma = c(0,0,3,0), bg = "white")
      for (cat in categories) {
        row_df <- all_by_type_radar %>% filter(type == cat)
        if (nrow(row_df) > 0) {
          fmsb::radarchart(build_single_radar_df(row_df), axistype = 1,
                           pcol = type_colors[[cat]], pfcol = scales::alpha(type_colors[[cat]], 0.35),
                           plwd = 3, cglcol = "gray80", cglty = 1, cglwd = 0.8,
                           axislabcol = "gray40", vlcex = 1.4,
                           caxislabels = c("0","25","50","75","100"), title = "")
          mtext(cat, side = 3, line = 1.4, cex = 1.0, font = 2)
        } else {
          plot.new(); text(0.5, 0.5, paste0(cat, "\n(insufficient data)"), col = "gray50")
        }
      }
      mtext("Average Metric Profile by Domain", outer = TRUE, cex = 1.4, font = 2, line = 0.8)
      par(mfrow = c(1,1))
      
    } else { # top_bottom_20
      
      categories   <- sort(unique(repo$type))
      n_cats_local <- length(categories)
      
      par(mfrow = c(n_cats_local, 2), mar = c(1,2,4,2), oma = c(0,0,3,0), bg = "white")
      for (cat in categories) {
        top_row <- top20_by_type_radar    %>% filter(type == cat)
        bot_row <- bottom20_by_type_radar %>% filter(type == cat)
        
        if (nrow(top_row) > 0) {
          fmsb::radarchart(build_single_radar_df(top_row), axistype = 1,
                           pcol = "#1F78B4", pfcol = scales::alpha("#1F78B4", 0.35),
                           plwd = 3, cglcol = "gray80", cglty = 1, cglwd = 0.8,
                           axislabcol = "gray40", vlcex = 1.3,
                           caxislabels = c("0","25","50","75","100"), title = "")
          mtext(paste0(cat, " — Top 20"), side = 3, line = 1.4, cex = 0.9, font = 2)
        } else {
          plot.new(); text(0.5, 0.5, "insufficient data", col = "gray50")
        }
        
        if (nrow(bot_row) > 0) {
          fmsb::radarchart(build_single_radar_df(bot_row), axistype = 1,
                           pcol = "#E31A1C", pfcol = scales::alpha("#E31A1C", 0.35),
                           plwd = 3, cglcol = "gray80", cglty = 1, cglwd = 0.8,
                           axislabcol = "gray40", vlcex = 1.3,
                           caxislabels = c("0","25","50","75","100"), title = "")
          mtext(paste0(cat, " — Bottom 20"), side = 3, line = 1.4, cex = 0.9, font = 2)
        } else {
          plot.new(); text(0.5, 0.5, "insufficient data", col = "gray50")
        }
      }
      mtext("Top 20 vs Bottom 20 by Domain", outer = TRUE, cex = 1.4, font = 2, line = 0.8)
      par(mfrow = c(1,1))
    }
  })
  
  # ── 8. Age cohorts ──────────────────────────────────────────────────────────
  output$cohort_ci_plot <- renderPlot({
    df <- filtered_repo()
    smry <- df %>%
      group_by(age_cohort,ci_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),n=n(),.groups="drop") %>%
      mutate(ci_label=ifelse(ci_present==1,"Has CI","No CI"))
    
    ggplot(smry,aes(x=age_cohort,y=mean_z,color=ci_label,group=ci_label)) +
      geom_line(linewidth=1.1) + geom_point(size=3) +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),width=0.15) +
      scale_color_manual(values=c("Has CI"="#1F78B4","No CI"="#A6CEE3")) +
      labs(title="Mean RQI by Repo Age Cohort — CI vs No CI",
           x="Repository Age Cohort",y="RQI (1–5)",color=NULL) +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(plot.title=element_text(face="bold",hjust=0.5))
  })
  
  output$cohort_test_plot <- renderPlot({
    df <- filtered_repo()
    smry <- df %>%
      group_by(age_cohort,tests_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()),n=n(),.groups="drop") %>%
      mutate(test_label=ifelse(tests_present==1,"Has Tests","No Tests"))
    
    ggplot(smry,aes(x=age_cohort,y=mean_z,color=test_label,group=test_label)) +
      geom_line(linewidth=1.1) + geom_point(size=3) +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z),width=0.15) +
      scale_color_manual(values=c("Has Tests"="#33A02C","No Tests"="#B2DF8A")) +
      labs(title="Mean RQI by Repo Age Cohort — Tests vs No Tests",
           x="Repository Age Cohort",y="RQI (1–5)",color=NULL) +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(plot.title=element_text(face="bold",hjust=0.5))
  })
  
  output$bio_cohort_plot <- renderPlot({
    bio_cohort <- repo %>%
      filter(type=="bioinformatics") %>%
      group_by(age_cohort) %>%
      summarise(ci_pct=mean(ci_present==1,na.rm=TRUE)*100,
                test_pct=mean(tests_present==1,na.rm=TRUE)*100,
                n=n(),.groups="drop") %>%
      pivot_longer(c(ci_pct,test_pct),names_to="metric",values_to="pct") %>%
      mutate(metric=recode(metric,"ci_pct"="CI Coverage","test_pct"="Test Coverage"))
    
    ggplot(bio_cohort,aes(x=age_cohort,y=pct,fill=metric)) +
      geom_col(position=position_dodge(0.6),width=0.55,color="white") +
      geom_text(aes(label=paste0(round(pct,1),"%")),
                position=position_dodge(0.6),vjust=-0.4,size=3.2) +
      scale_fill_manual(values=c("CI Coverage"="#1F78B4","Test Coverage"="#33A02C")) +
      labs(title="Bioinformatics: CI and Test Coverage by Age Cohort",
           x="Repository Age Cohort",y="Coverage (%)",fill=NULL) +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(plot.title=element_text(face="bold",hjust=0.5))
  })
  
  output$cohort_combined_plot <- renderPlot({
    df <- filtered_repo()
    
    cohort_ci_s <- df %>%
      group_by(age_cohort, ci_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()), n=n(), .groups="drop") %>%
      mutate(ci_label=ifelse(ci_present==1,"Has CI","No CI"))
    
    cohort_test_s <- df %>%
      group_by(age_cohort, tests_present) %>%
      summarise(mean_z=mean(z_normalized_1to5,na.rm=TRUE),
                se_z=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n()), n=n(), .groups="drop") %>%
      mutate(test_label=ifelse(tests_present==1,"Has Tests","No Tests"))
    
    base_theme <- theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(plot.title=element_text(face="bold",hjust=0.5),
            axis.text.x=element_text(angle=20,hjust=1,size=18),
            legend.title=element_blank())
    
    p_ci <- ggplot(cohort_ci_s, aes(x=age_cohort,y=mean_z,color=ci_label,group=ci_label)) +
      geom_line(linewidth=1.1) + geom_point(size=3) +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z), width=0.15) +
      scale_color_manual(values=c("Has CI"="#1F78B4","No CI"="#A6CEE3")) +
      labs(title="CI vs No CI", x="Repository Age Cohort", y="RQI (1–5)") +
      base_theme
    
    p_tests <- ggplot(cohort_test_s, aes(x=age_cohort,y=mean_z,color=test_label,group=test_label)) +
      geom_line(linewidth=1.1) + geom_point(size=3) +
      geom_errorbar(aes(ymin=mean_z-se_z,ymax=mean_z+se_z), width=0.15) +
      scale_color_manual(values=c("Has Tests"="#33A02C","No Tests"="#B2DF8A")) +
      labs(title="Tests vs No Tests", x="Repository Age Cohort", y="RQI (1–5)") +
      base_theme
    
    (p_ci | p_tests) +
      plot_annotation(title="Mean RQI by Repository Age Cohort", tag_levels="A",
                      theme=theme(plot.title=element_text(face="bold",hjust=0.5,size=24))) +
      plot_layout(guides="collect") &
      theme(legend.position="bottom")
  })
  
  # ── 9. Language ─────────────────────────────────────────────────────────────
  output$lang_plot <- renderPlot({
    if (input$lang_scope == "bio") {
      df           <- filter(lang_primary, type == "bioinformatics")
      min_n_use    <- min_n_bio
      valid_langs  <- bio_lang_summary_filt$language
      title_prefix <- "Bioinformatics: "
    } else {
      df           <- if (input$repo_type=="All") lang_primary else filter(lang_primary,type==input$repo_type)
      min_n_use    <- min_n
      valid_langs  <- lang_summary_filt$language
      title_prefix <- ""
    }
    df <- filter(df, language %in% valid_langs)
    
    if (input$lang_split == "none") {
      smry <- df %>% group_by(language) %>%
        summarise(n=n(), mean_rating=mean(z_normalized_1to5,na.rm=TRUE), .groups="drop")
      ggplot(smry,aes(x=fct_reorder(language,mean_rating),y=mean_rating)) +
        geom_col(fill="#4C72B0", width=0.7) +
        coord_flip() +
        labs(title=paste0(title_prefix,"Mean RQI by Primary Language (n≥",min_n_use,")"),
             x=NULL,y="Mean RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
      
    } else if (input$lang_split == "tests") {
      smry <- df %>% group_by(language,tests_present) %>%
        summarise(n=n(),mean_rating=mean(z_normalized_1to5,na.rm=TRUE),
                  se=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n),.groups="drop")
      ggplot(smry,aes(x=fct_reorder(language,mean_rating),y=mean_rating,
                      fill=factor(tests_present))) +
        geom_col(position=position_dodge(0.7),width=0.65,color="white") +
        geom_errorbar(aes(ymin=mean_rating-se,ymax=mean_rating+se),
                      position=position_dodge(0.7),width=0.2) +
        coord_flip() +
        scale_fill_manual(values=c("0"="#B2DF8A","1"="#33A02C"),
                          labels=c("No Tests","Has Tests")) +
        labs(title=paste0(title_prefix,"Mean RQI by Primary Language (split by Testing) — n≥", min_n_use),
             x=NULL,y="Mean RQI (1–5)",fill="Tests Present") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
      
    } else if (input$lang_split == "ci") {
      smry <- df %>% group_by(language,ci_present) %>%
        summarise(n=n(),mean_rating=mean(z_normalized_1to5,na.rm=TRUE),
                  se=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n),.groups="drop")
      ggplot(smry,aes(x=fct_reorder(language,mean_rating),y=mean_rating,
                      fill=factor(ci_present))) +
        geom_col(position=position_dodge(0.7),width=0.65,color="white") +
        geom_errorbar(aes(ymin=mean_rating-se,ymax=mean_rating+se),
                      position=position_dodge(0.7),width=0.2) +
        coord_flip() +
        scale_fill_manual(values=c("0"="#A6CEE3","1"="#1F78B4"),
                          labels=c("No CI","Has CI")) +
        labs(title=paste0(title_prefix,"Mean RQI by Primary Language (split by CI) — n≥", min_n_use),
             x=NULL,y="Mean RQI (1–5)",fill="CI Present") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
      
    } else {
      smry <- df %>%
        mutate(qa_bucket=factor(case_when(
          ci_present==1 & tests_present==1 ~ "Has CI + Tests",
          ci_present==1 & tests_present==0 ~ "CI only",
          ci_present==0 & tests_present==1 ~ "Tests only",
          TRUE ~ "Neither"
        ), levels=c("Neither","Tests only","CI only","Has CI + Tests"))) %>%
        group_by(language,qa_bucket) %>%
        summarise(n=n(),mean_rating=mean(z_normalized_1to5,na.rm=TRUE),
                  se=sd(z_normalized_1to5,na.rm=TRUE)/sqrt(n),.groups="drop")
      ggplot(smry,aes(x=fct_reorder(language,mean_rating),y=mean_rating,fill=qa_bucket)) +
        geom_col(position=position_dodge(0.8),width=0.7,color="white") +
        geom_errorbar(aes(ymin=mean_rating-se,ymax=mean_rating+se),
                      position=position_dodge(0.8),width=0.2) +
        coord_flip() +
        labs(title=paste0(title_prefix,"Mean RQI by Primary Language (CI/Tests buckets) — n≥", min_n_use),
             x=NULL,y="Mean RQI (1–5)",fill="QA Bucket") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))
    }
  })
  
  # ── 10. Activity / survival ─────────────────────────────────────────────────
  output$survival_plot <- renderPlot({
    surv_model <- glm(
      is_active ~ ci_present + tests_present +
        log1p(repo_age_days) + log1p(commit_count) + type,
      data=repo, family=binomial
    )
    tidy_surv <- broom::tidy(surv_model, conf.int=TRUE, exponentiate=TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(significant = p.value < 0.05,
             term = recode(term,
                           "ci_present"="CI Present","tests_present"="Tests Present",
                           "log1p(repo_age_days)"="log(Repo Age)","log1p(commit_count)"="log(Commit Count)",
                           "typebioinformatics"="Type: Bioinformatics",
                           "typeastrophysics"="Type: Astrophysics",
                           "typeimage_recognition"="Type: Image Recognition",
                           "typeopen_source"="Type: Open Source",
                           "typesoftware_engineering"="Type: Software Engineering"))
    
    ggplot(tidy_surv,aes(x=estimate,y=fct_reorder(term,estimate),color=significant)) +
      geom_point(size=3) +
      geom_errorbarh(aes(xmin=conf.low,xmax=conf.high),height=0.2) +
      geom_vline(xintercept=1,linetype="dashed",color="gray40") +
      scale_color_manual(values=c("FALSE"="gray60","TRUE"="#E31A1C"),
                         labels=c("p ≥ 0.05","p < 0.05")) +
      labs(title="Odds Ratios: Predictors of Repository Activity",
           subtitle="Outcome: committed within last 365 days  |  OR > 1 = higher odds of being active",
           x="Odds Ratio (±95% CI)",y=NULL,color="Significance") +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(plot.title=element_text(face="bold",hjust=0.5),
            plot.subtitle=element_text(hjust=0.5,color="gray40"))
  })
  
  output$activity_plot <- renderPlot({
    act_smry <- repo %>%
      group_by(type,ci_present,tests_present) %>%
      summarise(activity_rate=mean(is_active,na.rm=TRUE)*100,n=n(),.groups="drop") %>%
      mutate(qa_bucket=factor(case_when(
        ci_present==1&tests_present==1~"CI + Tests",
        ci_present==1&tests_present==0~"CI Only",
        ci_present==0&tests_present==1~"Tests Only",
        TRUE~"Neither"
      ),levels=c("Neither","Tests Only","CI Only","CI + Tests")))
    
    ggplot(act_smry,aes(x=qa_bucket,y=activity_rate,fill=qa_bucket)) +
      geom_col(width=0.65,color="white") +
      geom_text(aes(label=paste0(round(activity_rate,1),"%")),vjust=-0.4,size=3.2) +
      facet_wrap(~type) +
      scale_fill_manual(values=c("Neither"="#FDBF6F","Tests Only"="#B2DF8A",
                                 "CI Only"="#A6CEE3","CI + Tests"="#1F78B4")) +
      labs(title="Repository Activity Rate by CI/Test Presence and Domain",
           x=NULL,y="Active Repos (%)",fill=NULL) +
      theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
      theme(plot.title=element_text(face="bold",hjust=0.5),
            axis.text.x=element_text(angle=30,hjust=1,size=18),
            legend.position="none")
  })
  
  # ── 11. Citations (conditional) ─────────────────────────────────────────────
  output$citation_ui <- renderUI({
    if (!has_citations) {
      p("Citation data not found (repos_with_citations_gemma_v2.csv).", style="color:gray;")
    } else {
      tagList(
        plotOutput("cit_plot_out", height="520px"),
        hr(),
        h4("Table S2 — Top 20 Most-Cited Repositories with RQI ≤ 2.5"),
        tableOutput("citation_table_s2")
      )
    }
  })
  
  output$cit_plot_out <- renderPlot({
    req(has_citations)
    craw <- citations_raw %>% filter(!is.na(citation_count), citation_count >= 0)
    
    if (input$cit_plot == "loess") {
      dat <- craw %>% mutate(rqi_bin=round(z_normalized_1to5*4)/4) %>%
        group_by(rqi_bin) %>% summarise(mean_citations=mean(citation_count),.groups="drop")
      ggplot(dat,aes(x=rqi_bin,y=mean_citations)) +
        geom_smooth(method="loess",span=0.2,color="#1F78B4",fill="#1F78B4",
                    alpha=0.15,linewidth=1.2,se=FALSE) +
        geom_point(color="#1F78B4",size=2) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(title="Mean Citation Count vs Repository Quality",x="RQI (1–5)",y="Mean Citations") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "median_loess") {
      dat <- craw %>% mutate(rqi_bin=round(z_normalized_1to5*2)/2) %>%
        group_by(rqi_bin) %>% summarise(median_citations=median(citation_count),.groups="drop")
      ggplot(dat,aes(x=rqi_bin,y=median_citations)) +
        geom_smooth(method="loess",span=0.5,color="#1F78B4",fill="#1F78B4",
                    alpha=0.15,linewidth=1.2,se=FALSE) +
        geom_point(color="#1F78B4",size=2) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(title="Median Citation Count vs Repository Quality",x="RQI (1–5)",y="Median Citations") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "violin") {
      dat <- craw %>% mutate(rqi_bin=factor(round(z_normalized_1to5*2)/2))
      ggplot(dat,aes(x=rqi_bin,y=log1p(citation_count),fill=rqi_bin)) +
        geom_violin(alpha=0.7,trim=TRUE) +
        geom_boxplot(width=0.05,fill="white",outlier.shape=NA,color="grey30") +
        geom_jitter(shape=16,position=position_jitter(0.2)) +
        scale_fill_viridis_d(option="plasma") +
        labs(title="Citation Distribution by RQI Bin",
             x="RQI (0.5 intervals)",y="log(Citations + 1)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(legend.position="none",plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "mean_bar") {
      dat <- craw %>% mutate(rqi_bin = round(z_normalized_1to5*2)/2)
      ggplot(dat, aes(x=rqi_bin, y=citation_count)) +
        stat_summary(fun="mean", geom="bar", fill="#1F78B4", alpha=0.8, width=0.4) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(title="Mean Citation Count by RQI Bin", x="RQI (1–5)", y="Mean Citation Count") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "bar_lh") {
      grp <- craw %>% filter(!is.na(z_normalized_1to5)) %>%
        mutate(rqi_group=factor(ifelse(z_normalized_1to5<=2.5,"Low","High"),levels=c("Low","High")))
      grp$rqi_group <- require_two_groups(grp$rqi_group, "RQI group")
      smry  <- grp %>% group_by(rqi_group) %>% summarise(median_cit=median(citation_count),.groups="drop")
      p_val <- wilcox.test(citation_count~rqi_group,data=grp)$p.value
      sig   <- ifelse(p_val<0.001,"***",ifelse(p_val<0.01,"**",ifelse(p_val<0.05,"*","ns")))
      max_b <- max(smry$median_cit); step <- max_b*0.12
      ann   <- data.frame(x=1,xend=2,y=max_b+step,label=sig)
      ggplot(smry,aes(x=rqi_group,y=median_cit,fill=rqi_group)) +
        geom_col(alpha=0.85,width=0.5) +
        geom_text(aes(label=round(median_cit,1)),vjust=-0.5,size=4) +
        geom_segment(data=ann,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,linewidth=0.5) +
        geom_text(data=ann,aes(x=(x+xend)/2,y=y+step*0.1,label=label),inherit.aes=FALSE,size=4) +
        scale_fill_manual(values=c("Low"="#D73027","High"="#1A9850")) +
        labs(title="Median Citations by Quality Group",x="RQI Group",y="Median Citations") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(legend.position="none",plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "bar_bins") {
      grp <- craw %>% filter(!is.na(z_normalized_1to5)) %>%
        mutate(rqi_group=factor(case_when(
          z_normalized_1to5<2~"1-2",z_normalized_1to5<3~"2-3",
          z_normalized_1to5<4~"3-4",TRUE~"4-5"),levels=c("1-2","2-3","3-4","4-5")))
      smry <- grp %>% group_by(rqi_group) %>% summarise(median_cit=median(citation_count),.groups="drop")
      all_pairs <- build_sig_brackets(grp, "rqi_group", "citation_count")
      max_bar <- max(smry$median_cit); step <- max_bar*0.10
      if (nrow(all_pairs) > 0) all_pairs$y <- max_bar + step*seq_len(nrow(all_pairs))
      p <- ggplot(smry,aes(x=rqi_group,y=median_cit,fill=rqi_group)) +
        geom_col(alpha=0.85,width=0.6) +
        geom_text(aes(label=round(median_cit,1)),vjust=-0.5,size=4) +
        scale_fill_manual(values=c("1-2"="#FC8D59","2-3"="#FEE08B","3-4"="#91CF60","4-5"="#1A9850")) +
        labs(title="Median Citations by RQI Bin",x="RQI Score Bin",y="Median Citations") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(legend.position="none",plot.title=element_text(face="bold",hjust=0.5))
      if (nrow(all_pairs) > 0) {
        p <- p +
          geom_segment(data=all_pairs,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,linewidth=0.5) +
          geom_text(data=all_pairs,aes(x=(x+xend)/2,y=y+step*0.15,label=label),inherit.aes=FALSE,size=3.5) +
          coord_cartesian(ylim=c(0, max_bar + step*(nrow(all_pairs)+2)))
      }
      p
      
    } else if (input$cit_plot == "bar_bins_fine") {
      grp <- craw %>% filter(!is.na(z_normalized_1to5)) %>%
        mutate(
          rqi_group = case_when(
            z_normalized_1to5 < 1.5 ~ "1-1.5", z_normalized_1to5 < 2.0 ~ "1.5-2",
            z_normalized_1to5 < 2.5 ~ "2-2.5", z_normalized_1to5 < 3.0 ~ "2.5-3",
            z_normalized_1to5 < 3.5 ~ "3-3.5", z_normalized_1to5 < 4.0 ~ "3.5-4",
            z_normalized_1to5 < 4.5 ~ "4-4.5", TRUE ~ "4.5-5"
          ),
          rqi_group = factor(rqi_group, levels=c("1-1.5","1.5-2","2-2.5","2.5-3",
                                                 "3-3.5","3.5-4","4-4.5","4.5-5")),
          quality = factor(ifelse(z_normalized_1to5<2.5,"Low","High"), levels=c("Low","High"))
        )
      smry <- grp %>% group_by(rqi_group,quality) %>% summarise(median_cit=median(citation_count),.groups="drop")
      all_pairs <- build_sig_brackets(grp,"rqi_group","citation_count")
      max_bar <- max(smry$median_cit); step <- max_bar*0.10
      if (nrow(all_pairs) > 0) all_pairs$y <- max_bar + step*seq_len(nrow(all_pairs))
      grp$quality <- require_two_groups(grp$quality, "quality group")
      p_lh <- wilcox.test(citation_count~quality, data=grp)$p.value
      med_low <- median(grp$citation_count[grp$quality=="Low"])
      med_high <- median(grp$citation_count[grp$quality=="High"])
      sig_lbl <- if (p_lh<0.001) "***" else if(p_lh<0.01) "**" else if(p_lh<0.05) "*" else "ns"
      subtitle_str <- sprintf("Low (≤2.5) median = %.1f vs High (>2.5) median = %.1f | Wilcoxon %s (p%s)",
                              med_low, med_high, sig_lbl, ifelse(p_lh<0.001,"<0.001",sprintf("=%.3f",p_lh)))
      p <- ggplot(smry, aes(x=rqi_group,y=median_cit,fill=quality)) +
        geom_col(alpha=0.85,width=0.6) +
        geom_text(aes(label=round(median_cit,1)), vjust=-0.5, size=4) +
        scale_fill_manual(values=c("Low"="#D73027","High"="#1A9850")) +
        labs(title="Median Citations by RQI Bin", subtitle=subtitle_str,
             x="RQI Score Bin", y="Median Citations", fill="Quality Group") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5),
              plot.subtitle=element_text(hjust=0.5,color="gray40",size=16),
              legend.position="bottom")
      if (nrow(all_pairs) > 0) {
        p <- p +
          geom_segment(data=all_pairs, aes(x=x,xend=xend,y=y,yend=y), inherit.aes=FALSE, linewidth=0.5) +
          geom_text(data=all_pairs, aes(x=(x+xend)/2, y=y+step*0.15, label=label), inherit.aes=FALSE, size=3.5) +
          coord_cartesian(ylim=c(0, max_bar + step*(nrow(all_pairs)+2)))
      }
      p
      
    } else if (input$cit_plot == "fig6") {
      grp <- craw %>% filter(citation_count>0, !is.na(z_normalized_1to5), published) %>%
        mutate(
          rqi_group = case_when(
            z_normalized_1to5<2~"1-2", z_normalized_1to5<2.5~"2-2.5",
            z_normalized_1to5<3.5~"2.5-3.5", z_normalized_1to5<4.5~"3.5-4.5", TRUE~"4.5-5"
          ),
          rqi_group = factor(rqi_group, levels=c("1-2","2-2.5","2.5-3.5","3.5-4.5","4.5-5")),
          quality = factor(ifelse(z_normalized_1to5<=2.5,"Low","High"), levels=c("Low","High"))
        )
      smry <- grp %>% group_by(rqi_group,quality) %>% summarise(median_cit=median(citation_count),.groups="drop")
      all_pairs <- build_sig_brackets(grp,"rqi_group","citation_count")
      max_bar <- max(smry$median_cit); step <- max_bar*0.10
      if (nrow(all_pairs) > 0) all_pairs$y <- max_bar + step*seq_len(nrow(all_pairs))
      grp$quality <- require_two_groups(grp$quality, "quality group")
      p_lh <- wilcox.test(citation_count~quality, data=grp)$p.value
      med_low <- median(grp$citation_count[grp$quality=="Low"])
      med_high <- median(grp$citation_count[grp$quality=="High"])
      sig_lbl <- if (p_lh<0.001) "***" else if(p_lh<0.01) "**" else if(p_lh<0.05) "*" else "ns"
      subtitle_str <- sprintf("Low (≤2.5) median = %.1f vs High (>2.5) median = %.1f | Wilcoxon %s (p%s)",
                              med_low, med_high, sig_lbl, ifelse(p_lh<0.001,"<0.001",sprintf("=%.3f",p_lh)))
      p <- ggplot(smry, aes(x=rqi_group,y=median_cit,fill=quality)) +
        geom_col(alpha=0.85,width=0.6) +
        geom_text(aes(label=round(median_cit,1)), vjust=-0.5, size=4) +
        scale_fill_manual(values=c("Low"="#D73027","High"="#1A9850")) +
        labs(title="Median Citation Count by RQI Bin (Published Repositories)",
             subtitle=subtitle_str, x="RQI Score Bin", y="Median Citation Count", fill="Quality Group") +
        theme_minimal(base_size=22) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5),
              plot.subtitle=element_text(hjust=0.5,color="gray40",size=16),
              legend.position="bottom")
      if (nrow(all_pairs) > 0) {
        p <- p +
          geom_segment(data=all_pairs, aes(x=x,xend=xend,y=y,yend=y), inherit.aes=FALSE, linewidth=0.5) +
          geom_text(data=all_pairs, aes(x=(x+xend)/2, y=y+step*0.15, label=label), inherit.aes=FALSE, size=3.5) +
          coord_cartesian(ylim=c(0, max_bar + step*(nrow(all_pairs)+2)))
      }
      p
      
    } else if (input$cit_plot == "bar_ci") {
      dat <- craw %>% mutate(ci_label=if_else(has_CI==TRUE,"Has CI","No CI"))
      t_r <- t.test(citation_count~ci_label,data=dat)
      ggplot(dat,aes(x=ci_label,y=citation_count,fill=ci_label)) +
        geom_bar(stat="summary",fun="mean",alpha=0.8,width=0.5) +
        geom_jitter(alpha=0.2,width=0.15,size=1.5) +
        scale_fill_manual(values=c("Has CI"="#1F78B4","No CI"="#A6CEE3")) +
        coord_cartesian(ylim=c(0,quantile(dat$citation_count,0.95))) +
        labs(title="Mean Citations by CI Presence",
             subtitle=paste("t-test p =",round(t_r$p.value,4)),x=NULL,y="Citation Count") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(legend.position="none",plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "pub_bar") {
      # Uses citations_raw (not craw) — citation_count is only populated for
      # WoS-matched repos, so filtering on it here would collapse pub_label
      # down to "Published" only, discarding almost every "Not Published" row.
      bio_pub <- citations_raw %>% filter(tolower(type)=="bioinformatics")
      bio_pub$pub_label <- require_two_groups(bio_pub$pub_label, "publication-status")
      wt <- wilcox.test(z_normalized_1to5 ~ pub_label, data=bio_pub)
      meds <- bio_pub %>% group_by(pub_label) %>%
        summarise(median_rqi=median(z_normalized_1to5,na.rm=TRUE),.groups="drop")
      y_max <- max(meds$median_rqi); sig_y <- y_max*1.1
      ggplot(bio_pub, aes(x=pub_label,y=z_normalized_1to5,fill=pub_label)) +
        geom_bar(stat="summary", fun="median", alpha=0.8, width=0.5) +
        geom_text(data=meds, aes(x=pub_label,y=median_rqi,label=round(median_rqi,2)),
                  vjust=-0.5, fontface="bold", inherit.aes=FALSE) +
        annotate("segment", x=1, xend=2, y=sig_y, yend=sig_y, linewidth=0.5) +
        annotate("text", x=1.5, y=sig_y*1.05,
                 label=ifelse(wt$p.value<0.0001,"p < 0.0001",paste("p =",round(wt$p.value,4))),
                 hjust=0.5, fontface="bold") +
        scale_fill_manual(values=c("Published"="#1F78B4","Not Published"="#A6CEE3")) +
        coord_cartesian(ylim=c(0,sig_y*1.15)) +
        labs(title="Median RQI by Publication Status (Bioinformatics)", x=NULL, y="Median RQI") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(legend.position="none", plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "pub_violin") {
      # Uses citations_raw (not craw) — see note in pub_bar above.
      bio_pub <- citations_raw %>% filter(tolower(type)=="bioinformatics")
      bio_pub$pub_label <- require_two_groups(bio_pub$pub_label, "publication-status")
      wt <- wilcox.test(z_normalized_1to5 ~ pub_label, data=bio_pub)
      meds <- bio_pub %>% group_by(pub_label) %>%
        summarise(median_rqi=median(z_normalized_1to5,na.rm=TRUE), n=n(), .groups="drop")
      p_label <- ifelse(wt$p.value<0.0001,"p < 0.0001", paste("p =",round(wt$p.value,4)))
      y_max <- max(bio_pub$z_normalized_1to5, na.rm=TRUE); sig_y <- y_max*1.08
      ggplot(bio_pub, aes(x=pub_label,y=z_normalized_1to5,fill=pub_label)) +
        geom_violin(alpha=0.7, trim=TRUE, linewidth=0.4) +
        geom_jitter(aes(color=pub_label), alpha=0.25, width=0.15, size=1.2, show.legend=FALSE) +
        geom_boxplot(width=0.12, fill="white", outlier.shape=NA, linewidth=0.5, coef=0) +
        geom_point(data=meds, aes(x=pub_label,y=median_rqi), shape=18, size=3.5, color="grey20", inherit.aes=FALSE) +
        geom_text(data=meds, aes(x=pub_label,y=median_rqi,label=paste0(round(median_rqi,2),"\n(n=",n,")")),
                  vjust=-0.6, fontface="bold", size=3.5, inherit.aes=FALSE) +
        annotate("segment", x=1,xend=2,y=sig_y,yend=sig_y, linewidth=0.5) +
        annotate("text", x=1.5,y=sig_y*1.03,label=p_label, hjust=0.5, fontface="bold", size=3.5) +
        scale_fill_manual(values=c("Published"="#1F78B4","Not Published"="#A6CEE3")) +
        scale_color_manual(values=c("Published"="#1F78B4","Not Published"="#A6CEE3")) +
        coord_cartesian(ylim=c(0,sig_y*1.12)) +
        labs(title="RQI Distribution by Publication Status (Bioinformatics)", x=NULL, y="RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(legend.position="none", plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "fig5a") {
      # Uses citations_raw (not craw) — see note in pub_bar above.
      base_data <- citations_raw %>%
        filter(tolower(type)=="bioinformatics", !is.na(z_normalized_1to5)) %>%
        mutate(
          pub_label    = factor(pub_label, levels=c("Published","Not Published")),
          recency_1to5 = scale_1to5(recency_z),
          issue_1to5   = scale_1to5(issue_z),
          pop_1to5     = scale_1to5(pop_z)
        )
      comp_long <- base_data %>%
        select(pub_label, recency_1to5, issue_1to5, pop_1to5) %>%
        pivot_longer(-pub_label, names_to="component", values_to="score") %>%
        mutate(component=recode(component, recency_1to5="Recency", issue_1to5="Activity", pop_1to5="Popularity"))
      comp_stats <- comp_long %>% group_by(pub_label,component) %>% summarise(med=median(score,na.rm=TRUE),.groups="drop")
      rqi_stats  <- base_data %>% group_by(pub_label) %>% summarise(med=median(z_normalized_1to5,na.rm=TRUE), n=n(), .groups="drop")
      invisible(require_two_groups(base_data$pub_label, "publication-status"))
      wt <- wilcox.test(z_normalized_1to5~pub_label, data=base_data, exact=FALSE)
      p_label <- ifelse(wt$p.value<0.0001,"p < 0.0001", paste0("p = ",round(wt$p.value,4)))
      comp_colors <- c(Recency="#E31A1C", Activity="#33A02C", Popularity="#FF7F00")
      deep_blue <- "#08519C"
      fill_values <- c(Published=deep_blue, `Not Published`=deep_blue, RQI=deep_blue, comp_colors)
      
      ggplot() +
        geom_violin(data=base_data, aes(x=pub_label,y=z_normalized_1to5,fill=pub_label),
                    trim=TRUE, width=0.75, alpha=0.20, linewidth=0.4, color="grey60") +
        geom_violin(data=comp_long %>% filter(component=="Popularity"), aes(x=pub_label,y=score),
                    fill=comp_colors[["Popularity"]], color=NA, trim=TRUE, width=0.55, alpha=0.55) +
        geom_violin(data=comp_long %>% filter(component=="Activity"), aes(x=pub_label,y=score),
                    fill=comp_colors[["Activity"]], color=NA, trim=TRUE, width=0.40, alpha=0.65) +
        geom_violin(data=comp_long %>% filter(component=="Recency"), aes(x=pub_label,y=score),
                    fill=comp_colors[["Recency"]], color=NA, trim=TRUE, width=0.25, alpha=0.75) +
        geom_point(data=comp_stats, aes(x=pub_label,y=med,fill=component), shape=23, size=3.5, color="black", stroke=0.6) +
        geom_text(data=comp_stats, aes(x=pub_label,y=med,label=round(med,2)),
                  color="black", hjust=-0.45, fontface="bold", size=4, inherit.aes=FALSE, show.legend=FALSE) +
        geom_text(data=rqi_stats, aes(x=pub_label,y=5.3,label=paste0("n=",n)),
                  vjust=0, fontface="bold", size=4.5, color="grey40", inherit.aes=FALSE) +
        geom_point(data=rqi_stats, aes(x=pub_label,y=med,fill="RQI"), shape=21, size=4.5, color="black", stroke=0.7) +
        geom_text(data=rqi_stats, aes(x=pub_label,y=med,label=round(med,2)),
                  hjust=1.45, fontface="bold", size=4, color="black", inherit.aes=FALSE, show.legend=FALSE) +
        annotate("text", x=1.5,y=5.15, label=paste("Wilcoxon (RQI)",p_label), hjust=0.5, fontface="italic", size=3.6, color="grey35") +
        scale_fill_manual(values=fill_values, breaks=c("RQI","Recency","Activity","Popularity"), name=NULL) +
        scale_y_continuous(limits=c(1,5.3), breaks=1:5) +
        labs(title="RQI and Components by Publication Status (Bioinformatics)",
             subtitle="Components scaled 1–5", x=NULL, y="Score (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5,size=24),
              plot.subtitle=element_text(hjust=0.5,color="grey40",size=16),
              panel.grid.major.x=element_blank(), legend.position="bottom", legend.title=element_blank()) +
        guides(fill=guide_legend(override.aes=list(shape=c(21,22,22,22), size=4, color="black", stroke=0.6)))
      
    } else if (input$cit_plot == "age_scatter") {
      model_data <- craw %>%
        filter(!is.na(z_normalized_1to5), !is.na(repo_age_days), !is.na(citation_count),
               published, citation_count>0) %>%
        mutate(log_citations=log1p(citation_count), age_years=repo_age_days/365)
      ggplot(model_data, aes(x=z_normalized_1to5,y=log_citations,color=age_years)) +
        geom_point(alpha=0.5,size=2) +
        geom_smooth(method="lm", color="grey20", linewidth=1, se=TRUE) +
        scale_color_viridis_c(name="Repo Age (years)", option="plasma") +
        labs(title="Citation Count Explained by RQI and Repository Age",
             x="RQI (1–5)", y="log(Citation Count + 1)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) + theme(plot.title=element_text(face="bold",hjust=0.5))
      
    } else if (input$cit_plot == "scatter_lm") {
      dat <- craw %>%
        filter(!is.na(z_normalized_1to5), published, citation_count>0) %>%
        mutate(log_citations=log1p(citation_count))
      ggplot(dat,aes(x=z_normalized_1to5,y=log_citations)) +
        geom_point(alpha=0.3,size=1.8,color="#1F78B4") +
        geom_smooth(method="lm",color="#1F78B4",fill="#1F78B4",
                    alpha=0.15,linewidth=1.2,se=TRUE) +
        labs(title="Citation Count vs Repository Quality Index",
             x="RQI (1–5)",y="log(Citations + 1)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16))+theme(plot.title=element_text(face="bold",hjust=0.5))
    }
  })
  
  output$citation_table_s2 <- renderTable({
    req(has_citations)
    citations_raw %>%
      filter(!is.na(citation_count), !is.na(z_normalized_1to5)) %>%
      arrange(desc(citation_count), z_normalized_1to5) %>%
      filter(z_normalized_1to5 <= 2.5) %>%
      select(repo, type, citation_count, z_normalized_1to5, days_since_last_commit, doi) %>%
      head(20)
  }, striped = TRUE, hover = TRUE, digits = 2)
  
  # ── 12. Funding (conditional) ────────────────────────────────────────────────
  output$funding_ui <- renderUI({
    if (!has_funding) {
      p("Funding data not found (oa_publication_summary_gemma.csv).", style="color:gray;")
    } else {
      plotOutput("funding_plot_out", height="540px")
    }
  })
  
  output$funding_plot_out <- renderPlot({
    req(has_funding)
    
    if (input$fund_plot == "bar_bio") {
      bar_data <- funding_df %>%
        filter(tolower(type)=="bioinformatics",!is.na(z_normalized_1to5),!is.na(funding_group))
      bar_data$funding_group <- require_two_groups(bar_data$funding_group, "funding-status")
      w  <- wilcox.test(z_normalized_1to5~funding_group,data=bar_data,exact=FALSE)
      pl <- ifelse(w$p.value<0.0001,"p < 0.0001",paste("p =",round(w$p.value,4)))
      meds <- bar_data %>% group_by(funding_group) %>%
        summarise(median_rqi=median(z_normalized_1to5,na.rm=TRUE),n=n(),.groups="drop")
      y_max <- max(meds$median_rqi); sig_y <- y_max*1.12
      
      ggplot(bar_data,aes(x=funding_group,y=z_normalized_1to5,fill=funding_group)) +
        geom_bar(stat="summary",fun="median",alpha=0.85,width=0.5) +
        geom_jitter(alpha=0.12,width=0.18,size=1.2,color="grey30") +
        geom_text(data=meds,aes(x=funding_group,y=median_rqi,
                                label=paste0(round(median_rqi,2),"\n(n=",n,")")),
                  vjust=-0.4,fontface="bold",size=3.8,inherit.aes=FALSE) +
        annotate("segment",x=1,xend=2,y=sig_y,yend=sig_y,linewidth=0.5) +
        annotate("text",x=1.5,y=sig_y*1.04,label=pl,hjust=0.5,fontface="bold",size=3.5) +
        scale_fill_manual(values=c("Grant Funded"="#1F78B4","Not Grant Funded"="#A6CEE3")) +
        scale_y_continuous(limits=c(0,sig_y*1.12)) +
        labs(title="Median RQI by Funding Status (Bioinformatics)",
             subtitle="Wilcoxon rank-sum test",x=NULL,y="Median RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(legend.position="none",plot.title=element_text(face="bold",hjust=0.5),
              plot.subtitle=element_text(hjust=0.5,color="grey40"))
      
    } else if (input$fund_plot == "lm_dev_bio") {
      lm_data <- funding_df %>%
        filter(tolower(type)=="bioinformatics", !is.na(developer_cost_usd),
               !is.na(z_normalized_1to5), developer_cost_usd>0) %>%
        mutate(log_award=log10(developer_cost_usd))
      mod <- lm(z_normalized_1to5~log_award, data=lm_data); ms <- summary(mod)
      beta <- round(coef(mod)[["log_award"]],3); r2 <- round(ms$r.squared,3)
      pval <- round(coef(ms)[2,"Pr(>|t|)"],4)
      sub <- paste0("β=",beta," | R²=",r2," | ", ifelse(pval<0.0001,"p<0.0001",paste("p=",pval))," | n=",nrow(lm_data))
      ggplot(lm_data, aes(x=log_award,y=z_normalized_1to5)) +
        geom_point(alpha=0.55,size=2.2,color="#1F78B4") +
        geom_smooth(method="lm", color="#1F78B4", fill="#1F78B4", alpha=0.15, linewidth=1.2, se=TRUE) +
        scale_x_continuous(name="Developer Award Amount (log₁₀ USD)",
                           labels=function(x) paste0("$",format(10^x,big.mark=",",scientific=FALSE))) +
        labs(title="RQI vs Developer Award Amount (Bioinformatics)", subtitle=sub, y="RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5), plot.subtitle=element_text(hjust=0.5,color="grey40",size=16))
      
    } else if (input$fund_plot == "lm_org_bio") {
      lm_data <- funding_df %>%
        filter(tolower(type)=="bioinformatics", !is.na(total_cost_usd),
               !is.na(z_normalized_1to5), total_cost_usd>0) %>%
        mutate(log_award=log10(total_cost_usd))
      mod <- lm(z_normalized_1to5~log_award, data=lm_data); ms <- summary(mod)
      beta <- round(coef(mod)[["log_award"]],3); r2 <- round(ms$r.squared,3)
      pval <- round(coef(ms)[2,"Pr(>|t|)"],4)
      sub <- paste0("β=",beta," | R²=",r2," | ", ifelse(pval<0.0001,"p<0.0001",paste("p=",pval))," | n=",nrow(lm_data))
      ggplot(lm_data, aes(x=log_award,y=z_normalized_1to5)) +
        geom_point(alpha=0.55,size=2.2,color="#1F78B4") +
        geom_smooth(method="lm", color="#1F78B4", fill="#1F78B4", alpha=0.15, linewidth=1.2, se=TRUE) +
        scale_x_continuous(name="Total Organization Award Amount (log₁₀ USD)",
                           labels=function(x) paste0("$",format(10^x,big.mark=",",scientific=FALSE))) +
        labs(title="RQI vs Organization Award Amount (Bioinformatics)", subtitle=sub, y="RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5), plot.subtitle=element_text(hjust=0.5,color="grey40",size=16))
      
    } else if (input$fund_plot == "lm_dev_all") {
      lm_data <- funding_df %>%
        filter(!is.na(developer_cost_usd), !is.na(z_normalized_1to5), developer_cost_usd>0) %>%
        mutate(log_award=log10(developer_cost_usd))
      mod <- lm(z_normalized_1to5~log_award, data=lm_data); ms <- summary(mod)
      beta <- round(coef(mod)[["log_award"]],3); r2 <- round(ms$r.squared,3)
      pval <- round(coef(ms)[2,"Pr(>|t|)"],4)
      sub <- paste0("β=",beta," | R²=",r2," | ", ifelse(pval<0.0001,"p<0.0001",paste("p=",pval))," | n=",nrow(lm_data))
      ggplot(lm_data, aes(x=log_award,y=z_normalized_1to5)) +
        geom_point(alpha=0.55,size=2.2,color="#1F78B4") +
        geom_smooth(method="lm", color="#1F78B4", fill="#1F78B4", alpha=0.15, linewidth=1.2, se=TRUE) +
        scale_x_continuous(name="Developer Award Amount (log₁₀ USD)",
                           labels=function(x) paste0("$",format(10^x,big.mark=",",scientific=FALSE))) +
        labs(title="RQI vs Developer Award Amount (All Domains)", subtitle=sub, y="RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5), plot.subtitle=element_text(hjust=0.5,color="grey40",size=16))
      
    } else if (input$fund_plot == "lm_org_all") {
      lm_data <- funding_df %>%
        filter(!is.na(total_cost_usd), !is.na(z_normalized_1to5), total_cost_usd>0) %>%
        mutate(log_award=log10(total_cost_usd))
      mod <- lm(z_normalized_1to5~log_award, data=lm_data); ms <- summary(mod)
      beta <- round(coef(mod)[["log_award"]],3); r2 <- round(ms$r.squared,3)
      pval <- round(coef(ms)[2,"Pr(>|t|)"],4)
      sub <- paste0("β=",beta," | R²=",r2," | ", ifelse(pval<0.0001,"p<0.0001",paste("p=",pval))," | n=",nrow(lm_data))
      ggplot(lm_data, aes(x=log_award,y=z_normalized_1to5)) +
        geom_point(alpha=0.55,size=2.2,color="#1F78B4") +
        geom_smooth(method="lm", color="#1F78B4", fill="#1F78B4", alpha=0.15, linewidth=1.2, se=TRUE) +
        scale_x_continuous(name="Total Organization Award Amount (log₁₀ USD)",
                           labels=function(x) paste0("$",format(10^x,big.mark=",",scientific=FALSE))) +
        labs(title="RQI vs Organization Award Amount (All Domains)", subtitle=sub, y="RQI (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5), plot.subtitle=element_text(hjust=0.5,color="grey40",size=16))
      
    } else {
      # fig5b_bio (bioinformatics) or fig_all (all domains): layered RQI+component violins
      scope_bio <- (input$fund_plot == "fig5b_bio")
      
      base_data <- funding_df %>% filter(!is.na(z_normalized_1to5))
      if (scope_bio) base_data <- base_data %>% filter(tolower(type) == "bioinformatics")
      
      base_data <- base_data %>%
        mutate(
          recency_1to5  = scale_1to5(recency_z),
          issue_1to5    = scale_1to5(issue_z),
          pop_1to5      = scale_1to5(pop_z),
          funding_group = factor(funding_group, levels=c("Grant Funded","Not Grant Funded"))
        )
      
      comp_long <- base_data %>%
        select(funding_group, recency_1to5, issue_1to5, pop_1to5) %>%
        pivot_longer(-funding_group, names_to="component", values_to="score") %>%
        mutate(component=recode(component, recency_1to5="Recency", issue_1to5="Activity", pop_1to5="Popularity"))
      
      comp_stats <- comp_long %>% group_by(funding_group,component) %>%
        summarise(med=median(score,na.rm=TRUE),.groups="drop")
      rqi_stats  <- base_data %>% group_by(funding_group) %>%
        summarise(med=median(z_normalized_1to5,na.rm=TRUE), n=n(), .groups="drop")
      
      invisible(require_two_groups(base_data$funding_group, "funding-status"))
      wt <- wilcox.test(z_normalized_1to5~funding_group, data=base_data, exact=FALSE)
      p_label <- ifelse(wt$p.value<0.0001,"p < 0.0001", paste0("p = ",round(wt$p.value,4)))
      
      comp_colors <- c(Recency="#E31A1C", Activity="#33A02C", Popularity="#FF7F00")
      deep_blue <- "#08519C"
      fill_values <- c(`Grant Funded`=deep_blue, `Not Grant Funded`=deep_blue, RQI=deep_blue, comp_colors)
      
      title_str <- if (scope_bio) "RQI and Components by Funding Status (Bioinformatics)" else
        "RQI and Components by Funding Status (All Domains)"
      subtitle_str <- paste0("Components scaled 1–5  |  Grant sources: ",
                             paste(funded_sources_list, collapse=", "))
      
      ggplot() +
        geom_violin(data=base_data, aes(x=funding_group,y=z_normalized_1to5,fill=funding_group),
                    trim=TRUE, width=0.75, alpha=0.20, linewidth=0.4, color="grey60") +
        geom_violin(data=comp_long %>% filter(component=="Popularity"), aes(x=funding_group,y=score),
                    fill=comp_colors[["Popularity"]], color=NA, trim=TRUE, width=0.55, alpha=0.55) +
        geom_violin(data=comp_long %>% filter(component=="Activity"), aes(x=funding_group,y=score),
                    fill=comp_colors[["Activity"]], color=NA, trim=TRUE, width=0.40, alpha=0.65) +
        geom_violin(data=comp_long %>% filter(component=="Recency"), aes(x=funding_group,y=score),
                    fill=comp_colors[["Recency"]], color=NA, trim=TRUE, width=0.25, alpha=0.75) +
        geom_point(data=comp_stats, aes(x=funding_group,y=med,fill=component), shape=23, size=3.5, color="black", stroke=0.6) +
        geom_text(data=comp_stats, aes(x=funding_group,y=med,label=round(med,2)),
                  color="black", hjust=-0.45, fontface="bold", size=4, inherit.aes=FALSE, show.legend=FALSE) +
        geom_text(data=rqi_stats, aes(x=funding_group,y=5.3,label=paste0("n=",n)),
                  vjust=0, fontface="bold", size=4.5, color="grey40", inherit.aes=FALSE) +
        geom_point(data=rqi_stats, aes(x=funding_group,y=med,fill="RQI"), shape=21, size=4.5, color="black", stroke=0.7) +
        geom_text(data=rqi_stats, aes(x=funding_group,y=med,label=round(med,2)),
                  hjust=1.45, fontface="bold", size=4, color="black", inherit.aes=FALSE, show.legend=FALSE) +
        annotate("text", x=1.5,y=5.15, label=paste("Wilcoxon (RQI)",p_label), hjust=0.5, fontface="italic", size=3.6, color="grey35") +
        scale_fill_manual(values=fill_values, breaks=c("RQI","Recency","Activity","Popularity"), name=NULL) +
        scale_y_continuous(limits=c(1,5.3), breaks=1:5) +
        labs(title=title_str, subtitle=subtitle_str, x=NULL, y="Score (1–5)") +
        theme_minimal(base_size=18) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=24), plot.subtitle=element_text(size=16)) +
        theme(plot.title=element_text(face="bold",hjust=0.5,size=24),
              plot.subtitle=element_text(hjust=0.5,color="grey40",size=9),
              panel.grid.major.x=element_blank(), legend.position="bottom", legend.title=element_blank()) +
        guides(fill=guide_legend(override.aes=list(shape=c(21,22,22,22), size=4, color="black", stroke=0.6)))
    }
  })
  
  # ── 13. Summary tables ───────────────────────────────────────────────────────
  output$table1_summary <- renderTable(summary_tbl_type, digits = 2)
  output$tableS1_languages <- renderTable(top_languages_tbl)
}


shinyApp(
  ui = ui,
  server = server,
  options = list(host = "0.0.0.0", port = 3838)
)