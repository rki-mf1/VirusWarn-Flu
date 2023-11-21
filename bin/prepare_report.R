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
  }else if (is.null(args$file_variant_table)){
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

# Setup variables
file_variant_table = args$file_variant_table
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
GEOLOC_COL = "PRIMARY_DIAGNOSTIC_LAB_PLZ"
MUT_TYPE_COL = "variant_type"
VARIANT_CLASS_COL = "VariantType"

var_pheno = read_tsv(file_variant_table, col_type = "ccccicicc") %>%
  mutate(type = factor(type, levels = VARIANT_TYPE_LEVELS))

var_pheno_wide_filter = var_pheno %>%
  filter(variant_size <= indel_max_size) %>%
  mutate(Test = 1L) %>%
  pivot_wider(id_cols = ID:variant_size,
              names_from = "type",
              values_from = "Test",
              values_fn = max,
              values_fill = 0)

score_mutation <- function(pheno_table){
  if ((all(VARIANT_TYPE_LEVELS %in% colnames(pheno_table)))){
    pheno_table_infomask = pheno_table %>%
      mutate(infomask = str_c(
        #PositiveSelection,
        MutationOfConcern,
        RegionOfInterest,
        NotAnnotated, 
        sep = ".")
      )
    
    scores_combination = expand_grid(
      #PositiveSelection = c(TRUE,FALSE),
      MutationOfConcern = c(TRUE,FALSE),
      RegionOfInterest = c(TRUE,FALSE),
      NotAnnotated = c(TRUE,FALSE),
    ) %>%
    mutate(
      s_pm = case_when(
        NotAnnotated & !RegionOfInterest ~ 1,
        NotAnnotated & !MutationOfConcern ~ 1,
        !NotAnnotated & !MutationOfConcern & !RegionOfInterest ~ 1,
        TRUE ~ 0
      ),
      s_moc = case_when(
        !NotAnnotated & MutationOfConcern ~ 1,
        NotAnnotated & MutationOfConcern ~ 1,
        TRUE ~ 0
      ),
      s_roi = case_when(
        !NotAnnotated & RegionOfInterest ~ 1,
        NotAnnotated & RegionOfInterest ~ 1,
        !NotAnnotated & !MutationOfConcern & RegionOfInterest ~ 1,
        TRUE ~ 0
      )
    )
  } else {
    pheno_table_infomask = pheno_table %>%
      mutate(infomask = str_c(
        #PositiveSelection,
        RegionOfInterest,
        NotAnnotated, 
        sep = ".")
      )
    
    scores_combination = expand_grid(
      #PositiveSelection = c(TRUE,FALSE),
      RegionOfInterest = c(TRUE,FALSE),
      NotAnnotated = c(TRUE,FALSE)
    ) %>%
    mutate(
      s_pm = case_when(
        NotAnnotated & !RegionOfInterest ~ 1,
        !NotAnnotated & !RegionOfInterest ~ 1,
        TRUE ~ 0
      ),
      s_moc = case_when(
        TRUE ~ 0
      ),
      s_roi = case_when(
        NotAnnotated & RegionOfInterest ~ 1,
        !NotAnnotated & RegionOfInterest ~ 1,
        TRUE ~ 0
      )
    )
  }
  
  pheno_table_with_score = suppressMessages(inner_join(pheno_table_infomask, scores_combination))
  
  return(pheno_table_with_score)
}

var_pheno_wide_filter_with_score = score_mutation(var_pheno_wide_filter)

var_pheno_score_summary = suppressMessages(
  var_pheno_wide_filter_with_score %>%
  group_by(across(any_of(
   c(
     ID_COL,
     DATE_COL,
     GEOLOC_COL,
     MUT_TYPE_COL,
     VARIANT_CLASS_COL
   )
  ))) %>%
  summarise(
   ListMutationsSelected = str_c(aa_pattern[s_roi > 0 | s_moc > 0 | s_pm > 0],
                                 collapse = ","),
   nMutationsTotal = n(),
   across(starts_with("s_"), sum),
  ) %>%
  ungroup()
)

alert_colors = c(
  "red" = "#FF2400",     # alert 
  "orange" = "orange",   # weak alert
  "yellow" = "#FFEA00",   # accumulation alert
  "grey" = "slategrey"  # no alert
)
alert_level = factor(names(alert_colors), ordered = TRUE)

var_pheno_summary_wide = var_pheno_score_summary %>%
  pivot_wider(
    names_from = all_of(MUT_TYPE_COL),
    values_from = c(
      "ListMutationsSelected",
      "nMutationsTotal",
      "s_moc",
      "s_roi",
      "s_pm"
    ),
    values_fill = list(
      ListMutationsSelected = "",
      nMutationsTotal = 0 ,
      s_moc = 0 ,
      s_roi = 0,
      s_pm = 0
    )
  )

## This fonction returns the alert levels
## Version 1 with hard set thresholds
compute_alert_levels_v1 <- function(pheno_table_wide){
  ## WARNING: all variable names are hard coded in this function
  
  if(! "ListMutationsSelected_M" %in% colnames(pheno_table_wide)){pheno_table_wide$ListMutationsSelected_M = 0}
  if(! "ListMutationsSelected_D" %in% colnames(pheno_table_wide)){pheno_table_wide$ListMutationsSelected_D = 0}
  if(! "ListMutationsSelected_I" %in% colnames(pheno_table_wide)){pheno_table_wide$ListMutationsSelected_I = 0}
  
  if(! "nMutationsTotal_M" %in% colnames(pheno_table_wide)){pheno_table_wide$nMutationsTotal_M = 0}
  if(! "nMutationsTotal_D" %in% colnames(pheno_table_wide)){pheno_table_wide$nMutationsTotal_D = 0}
  if(! "nMutationsTotal_I" %in% colnames(pheno_table_wide)){pheno_table_wide$nMutationsTotal_I = 0}
  
  if(! "s_roi_M" %in% colnames(pheno_table_wide)){pheno_table_wide$s_roi_M = 0}
  if(! "s_roi_D" %in% colnames(pheno_table_wide)){pheno_table_wide$s_roi_D = 0}
  if(! "s_roi_I" %in% colnames(pheno_table_wide)){pheno_table_wide$s_roi_I = 0}
  
  if(! "s_pm_M" %in% colnames(pheno_table_wide)){pheno_table_wide$s_pm_M = 0}
  if(! "s_pm_D" %in% colnames(pheno_table_wide)){pheno_table_wide$s_pm_D = 0}
  if(! "s_pm_I" %in% colnames(pheno_table_wide)){pheno_table_wide$s_pm_I = 0}
  
  if(! "s_moc_M" %in% colnames(pheno_table_wide)){pheno_table_wide$s_moc_M = 0}
  if(! "s_moc_D" %in% colnames(pheno_table_wide)){pheno_table_wide$s_moc_D = 0}
  if(! "s_moc_I" %in% colnames(pheno_table_wide)){pheno_table_wide$s_moc_I = 0}
  
  pheno_table_wide_with_alert = pheno_table_wide %>%
    mutate(
      s_moc_roi_tot = (s_moc_M + s_moc_D + s_roi_M + s_roi_D + s_moc_I + s_roi_I),
      alert_level =
       case_when(
         (((s_moc_M + s_moc_D + s_moc_I) >= 1 & (s_moc_M + s_moc_D + s_moc_I + s_roi_M + s_roi_D + s_roi_I) >= 5)) ~ "red", # alert
         (((s_moc_M + s_moc_D + s_moc_I) >= 1 & (s_moc_M + s_moc_D + s_moc_I + s_roi_M + s_roi_D + s_roi_I) >= 3)) ~ "orange", # weak alert
         (((s_moc_M + s_moc_D + s_moc_I) == 0) & ((s_roi_M + s_roi_D + s_roi_I) >= 8)) ~ "yellow", # accumulation alert roi
         (((s_moc_M + s_moc_D + s_moc_I) == 0) & ((s_pm_M + s_pm_D + s_pm_I) >= 25)) ~ "yellow", # accumulation alert pm
         TRUE ~ "grey",
       ))
  return(pheno_table_wide_with_alert)
}

var_pheno_summary_wide_with_alert = compute_alert_levels_v1(var_pheno_summary_wide)
prediction_overview = var_pheno_summary_wide_with_alert %>% group_by(alert_level) %>% count()

# Write output with prediction overview
write.table(
  prediction_overview,
  file = file.path("prediction_overview.txt"),
  sep = "\t",
)

mutations_per_alert_level = suppressMessages(
  var_pheno_summary_wide_with_alert  %>%
  ungroup() %>%
  filter(alert_level != "grey") %>%
  select(ID, alert_level) %>%
  distinct() %>% ## WARNING THIS IS A HOTFIX
  inner_join(var_pheno_wide_filter %>% ungroup() %>% select(ID, aa_pattern)) %>%
  mutate(value = 1L) %>%
  group_by(alert_level) %>%
  nest() %>%
  mutate(data = map(
   data,
   ~ pivot_wider(
     .,
     names_from = aa_pattern,
     values_from = value,
     values_fill = 0
   )
  ))
)

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
  var_pheno_summary_wide_with_alert %>%
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
  inner_join(var_pheno_wide_filter) %>%
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
      s_moc_M, s_moc_D, s_moc_I,
      s_roi_M,	s_roi_D,	s_roi_I,
      s_pm_M, s_pm_D, s_pm_I,
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
