#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr)) 
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(logger, warn.conflicts = FALSE))

option_list <- list(
  make_option(
    c("-f", "--file_variant_table"),
    default = "variants_with_phenotypes.tsv",
    help = "File corresponding to the variants",
  ),
  make_option(
    c("-m", "--fixed_table"),
    default = "",
    help = "File containing the fixed mutations for each season",
  ),
  make_option(
    c("--season"),
    default = "23/24",
    help = "Season of the fixed mutations to consider",
  ),
  make_option(
    c("-a", "--file_annotations"),
    default = "",
    help = "[OPTIONAL] File containing metadata on the samples",
  ),
  make_option(
    c("-s", "--fluwarnsystem_alert_samples"), 
    default = "fluwarnsystem_alerts_samples_all.csv",
    help = "Output name for alert samples csv",
  ),
  make_option(
    c("-c", "--fluwarnsystem_alert_clusters"), 
    default = "fluwarnsystem_alerts_clusters_summaries_all.csv",
    help = "Output name for alert clusters csv",
  ),
  make_option(
    c("--verbose"), 
    default = FALSE, 
    action = "store_true",
    help = "Write more output to the console",
  ),
  make_option(
    c("--strict"), 
    default = "n", 
    help = "VirusWarn-Flu runs in strict mode",
  )
)

parser = OptionParser(option_list = option_list)
exit <- function() { invokeRestart("abort") }  
args <- parse_args(parser)

# Test if there is at least one argument and file_variant_file is supplied
args_tester <- function(args){
  if (length(args)==0) {
    print_help(parser)
    stop("At least one argument must be supplied.")
    exit()
  } else if (is.null(args$file_variant_table)) {
    print_help(parser)
    stop("Require --file_variant_table")
    exit()
  }
  return (args)
}
args = args_tester(args)

# Setup logging
if (args$verbose) {
  log_threshold(TRACE)
} else {
  log_threshold(WARN)
}

# ----- Functions ----- #
hamming <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}

hamming_tibble = function(tab) {
  mat = as.matrix(tab[, -1])
  hamming(mat)
}

get_sequence_clusters = function(tab, max_dist = 1) {
  matdists = hamming_tibble(tab)
  g = graph.adjacency(matdists <= max_dist, mode = "undirected")
  clus = components(g)
  tibble(
    ID = tab$ID,
    cluster_ID = clus$membership,
    cluster_size = clus$csize[clus$membership],
  )
}

compute_alert_levels_v1 <- function(pheno_table_wide){
  if (args$strict == "y") {
    pheno_table_wide_with_alert = pheno_table_wide %>%
      mutate(
        s_moc_roi_tot = (s_moc_M + s_moc_D + s_roi_M + s_roi_D + s_moc_I + s_roi_I + s_moc_F + s_roi_F),
        alert_level =
          case_when(
            (((s_moc_F >= 1 | s_roi_F >= 1) &
              (s_moc_M + s_moc_D + s_moc_I + s_moc_F) >= 1 & s_moc_roi_tot >= 3)) ~ "pink", # known variant
            (((s_moc_M + s_moc_D + s_moc_I) >= 1 & s_moc_roi_tot >= 5)) ~ "red", # alert
            (((s_moc_M + s_moc_D + s_moc_I) == 0) & 
              ((s_roi_M + s_roi_D + s_roi_I + s_roi_F) >= 8) & 
              ((s_pm_M + s_pm_D + s_pm_I + s_pm_F) < 25)) ~ "yellow", # accumulation alert roi
            (((s_moc_M + s_moc_D + s_moc_I) == 0) & 
              ((s_pm_M + s_pm_D + s_pm_I + s_pm_F) >= 25) & 
              ((s_roi_M + s_roi_D + s_roi_I + s_roi_F) < 8)) ~ "yellow", # accumulation alert pm
            TRUE ~ "grey"
          ))
  } else {
    pheno_table_wide_with_alert = pheno_table_wide %>%
      mutate(
        s_moc_roi_tot = (s_moc_M + s_moc_D + s_roi_M + s_roi_D + s_moc_I + s_roi_I + s_moc_F + s_roi_F),
        alert_level =
          case_when(
            (((s_moc_F >= 1 | s_roi_F >= 1) &
              (s_moc_M + s_moc_D + s_moc_I + s_moc_F) >= 1 & s_moc_roi_tot >= 3)) ~ "pink", # known variant
            (((s_moc_M + s_moc_D + s_moc_I) >= 1 & s_moc_roi_tot >= 5)) ~ "red", # alert
            (((s_moc_M + s_moc_D + s_moc_I) >= 1 & s_moc_roi_tot >= 3)) ~ "orange", # weak alert
            (((s_moc_M + s_moc_D + s_moc_I) == 0) & 
              ((s_roi_M + s_roi_D + s_roi_I + s_roi_F) >= 8) & 
              ((s_pm_M + s_pm_D + s_pm_I + s_pm_F) < 25)) ~ "yellow", # accumulation alert roi
            (((s_moc_M + s_moc_D + s_moc_I) == 0) & 
              ((s_pm_M + s_pm_D + s_pm_I + s_pm_F) >= 25) & 
              ((s_roi_M + s_roi_D + s_roi_I + s_roi_F) < 8)) ~ "yellow", # accumulation alert pm
            TRUE ~ "grey"
          ))
  }
  return(pheno_table_wide_with_alert)
}

# Setup variables
file_variant_table = args$file_variant_table
fixed_table = args$fixed_table
season = args$season
file_annotations = args$file_annotations

VARIANT_TYPE_LEVELS = c(
  #"PositiveSelection",
  "MutationOfConcern", 
  "RegionOfInterest",
  "NotAnnotated" 
)

indel_max_size = 5 # larger indel will not be considered
ID_COL = "ID"
DATE_COL = "SAMPLING_DATE"
MUT_TYPE_COL = "variant_type"
VARIANT_CLASS_COL = "VariantType"

var_pheno = read_tsv(file_variant_table, show_col_types = FALSE) %>%
  group_by(ID, aa_pattern) %>%
  filter(
    (n() == 1) |
      (type == "MutationOfConcern" & n() > 1) |
      (type == "RegionOfInterest" & n() > 1 & !any(type == "MutationOfConcern"))
  ) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = VARIANT_TYPE_LEVELS))

fixed = read_csv(fixed_table, show_col_types = FALSE) %>%
  filter(Season == season) %>%
  select(-...1)

var_pheno_fixed = var_pheno %>%
  mutate(
    type = ifelse(aa_pos_ref_start %in% fixed$Position, paste("Fixed", type, sep = ""), as.character(type)),
    variant_type = ifelse(aa_pos_ref_start %in% fixed$Position, "F", as.character(variant_type))
  ) 

var_pheno_final = var_pheno_fixed %>% 
  group_by(ID) %>%
  summarise(
    ListMutationsSelected = paste(aa_pattern, collapse = ","),
    ListMutationsSelected_M = paste(aa_pattern[variant_type == "M"], collapse = ","),
    ListMutationsSelected_D = paste(aa_pattern[variant_type == "D"], collapse = ","),
    ListMutationsSelected_I = paste(aa_pattern[variant_type == "I"], collapse = ","),
    ListMutationsSelected_F = paste(aa_pattern[variant_type == "F"], collapse = ","),
    nMutationsTotal = n(),
    nMutationsTotal_M = sum(variant_type == "M"),
    nMutationsTotal_D = sum(variant_type == "D"),
    nMutationsTotal_I = sum(variant_type == "I"),
    nMutationsTotal_F = sum(variant_type == "F"),
    s_moc_M = sum(type == "MutationOfConcern" & variant_type == "M"),
    s_moc_D = sum(type == "MutationOfConcern"  & variant_type == "D"),
    s_moc_I = sum(type == "MutationOfConcern" & variant_type == "I"),
    s_moc_F = sum(type == "FixedMutationOfConcern"),
    s_roi_M = sum(type == "RegionOfInterest" & variant_type == "M"),
    s_roi_D = sum(type == "RegionOfInterest" & variant_type == "D"),
    s_roi_I = sum(type == "RegionOfInterest" & variant_type == "I"),
    s_roi_F = sum(type == "FixedRegionOfInterest"),
    s_pm_M = sum(type == "NotAnnotated" & variant_type == "M"),
    s_pm_D = sum(type == "NotAnnotated" & variant_type == "D"),
    s_pm_I = sum(type == "NotAnnotated" & variant_type == "I"),
    s_pm_F = sum(type == "FixedNotAnnotated"),
    .groups = 'drop'
  )

if (args$strict == "y"){
  alert_colors = c(
    "pink" = "deeppink",    # known variant
    "red" = "#FF2400",      # alert 
    "yellow" = "#FFEA00",   # accumulation alert
    "grey" = "slategrey"    # no alert
  )
} else {
  alert_colors = c(
    "pink" = "deeppink",    # known variant
    "red" = "#FF2400",      # alert 
    "orange" = "orange",    # weak alert
    "yellow" = "#FFEA00",   # accumulation alert
    "grey" = "slategrey"    # no alert
  )
}
alert_level <- factor(names(alert_colors), ordered = TRUE)

var_pheno_final_with_alert <- compute_alert_levels_v1(var_pheno_final)
prediction_overview <- var_pheno_final_with_alert %>% group_by(alert_level) %>% count()

# Write output with prediction overview
write.table(
  prediction_overview,
  file = file.path("prediction_overview.txt"),
  sep = "\t",
)

filtered_df <- var_pheno_final_with_alert %>%
  filter(alert_level != "grey")
unique_colors <- unique(filtered_df$alert_level)

if (length(unique_colors) != 1) {
  mutations_per_alert_level <- suppressMessages(
    var_pheno_final_with_alert  %>%
    ungroup() %>%
    filter(alert_level != "grey") %>%
    select(ID, alert_level) %>%
    distinct() %>% ## WARNING THIS IS A HOTFIX
    inner_join(var_pheno_fixed %>% ungroup() %>% select(ID, aa_pattern)) %>%
    mutate(value = 1L) %>%
    group_by(alert_level) %>%
    nest()) %>% 
    mutate(data = map(data, ~ {
      pivot_wider(
        .x,
        names_from = aa_pattern,
        values_from = value,
        values_fill = list(value = 0)
      )
    }))
} else {
  mutations_per_alert_level <- suppressMessages(
    var_pheno_final_with_alert  %>%
    ungroup() %>%
    filter(alert_level != "grey") %>%
    select(ID, alert_level) %>%
    distinct() %>% ## WARNING THIS IS A HOTFIX
    inner_join(var_pheno_fixed %>% ungroup() %>% select(ID, aa_pattern)) %>%
    mutate(value = 1L) %>%
    group_by(alert_level) %>%
    nest())
    
  mutations_per_alert_level[[2]][[1]] <- mutations_per_alert_level[[2]][[1]] %>%
    group_by(ID, aa_pattern) %>%
    summarize(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(
      names_from = aa_pattern,
      values_from = value,
      values_fill = list(value = 0)
    )
}

alert_level_groups_with_clusters = mutations_per_alert_level %>%
  mutate(clusters = map(data, get_sequence_clusters))

if(nrow(alert_level_groups_with_clusters) != 0){
  alerts_with_clusters_ID = alert_level_groups_with_clusters %>% 
    select(-data) %>%
    unnest(c(alert_level, clusters)) %>%
    rename(cluster_ID_in_alert_level = cluster_ID)
} else {
  alerts_with_clusters_ID <- data.frame(
    alert_level=character(), 
    ID=character(), 
    cluster_ID_in_alert_level=character(),
    cluster_size=numeric())
}

samples_with_alert = suppressMessages(
  var_pheno_final_with_alert %>%
  left_join(alerts_with_clusters_ID) %>%
  arrange(desc(alert_level),
          desc(s_moc_roi_tot),
          desc(cluster_size),
          desc(DATE_COL))
)

log_info("*** Writing Results ***")

common_mutations_in_clusters = suppressMessages(
  samples_with_alert %>%
  filter(alert_level != "grey") %>%
  inner_join(var_pheno_fixed) %>%
  group_by(alert_level, cluster_ID_in_alert_level) %>%
  summarise(
    FreqMutThr = round( length( ID %>% unique() ) * 0.3, d = 0 ) + 1,
    aa_pattern_common = fct_lump_min(aa_pattern, FreqMutThr),
  ) %>%
  summarise(
    ListFrequentMutations_gt30perc = str_c(aa_pattern_common %>%
                                            unique())
  ) %>%
  ungroup() %>% 
    arrange(desc(alert_level), cluster_ID_in_alert_level),
)

list_clusters_properties = suppressMessages(
  samples_with_alert %>%
  filter(!is.na(cluster_ID_in_alert_level)) %>%
  group_by(alert_level, cluster_ID_in_alert_level, cluster_size) %>%
    summarise(
      across(starts_with("s_"), ~ round(mean(.), d = 1), .names = "{.col}.avg")
    )
)

list_clusters_properties_with_mutations = suppressMessages(
  common_mutations_in_clusters %>%
  inner_join(list_clusters_properties) %>%
  arrange(desc(alert_level), desc(s_moc_roi_tot.avg),) %>%
  select(!c(ends_with("_T.avg"))) %>%
  relocate(cluster_ID_in_alert_level, .after = last_col()) %>%
  rename(n_samples = cluster_size) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x , 0))),
)

# Write output with alert clusters
write_csv(
  list_clusters_properties_with_mutations,
  file = file.path(args$fluwarnsystem_alert_clusters),
)

samples_out = suppressMessages(
  samples_with_alert %>%
    mutate(
      across(where(is.numeric), ~ replace_na(.x , 0)),
      across(where(is.character), ~ replace_na(.x , "")),
    ) %>%
    select(
      !c(ends_with("_T.avg"),
         ends_with("MutationsSelected_T")),
    ) %>%
    relocate(
      alert_level,
      s_moc_roi_tot,
      s_moc_M, s_moc_D, s_moc_I, s_moc_F,
      s_roi_M,	s_roi_D,	s_roi_I, s_roi_F,
      s_pm_M, s_pm_D, s_pm_I, s_pm_F,
      starts_with("ListMutationsSelected"),
      any_of("LINEAGE.LATEST"),
      cluster_size,
      cluster_ID_in_alert_level),
    ) %>%
    arrange(
      across(c(
        alert_level, 
        s_moc_roi_tot,
        s_moc_D, s_moc_M,
        s_roi_D, s_roi_M))
)

# Write output with alert per sample
write_csv(
  samples_out,
  file = file.path(args$fluwarnsystem_alert_samples),
)

log_info("*** Success ***")
