# Function to filter bacterial genera and normalize abundance
de_other_domain_genus <- function(df) {
  df_s <- df %>% 
    tibble::rownames_to_column("genus") %>% 
    left_join(taxa %>% dplyr::select(domain, genus) %>% unique(), by = "genus") %>% 
    filter(domain == "Bacteria") %>% 
    tibble::column_to_rownames("genus") %>% 
    dplyr::select(-domain)
  
  # Normalize abundance data
  df_normalized <- apply(df_s, 2, function(x) x/sum(x))
  return(df_normalized)
}

# Set data directories (configure these paths as needed)
data_dir <- "path/to/your/data"
output_dir <- "path/to/your/output"

# Load and preprocess data
df_genus <- read.csv(file.path(data_dir, "relative_abundance_genus.csv"), row.names = 1)
bigtable <- read.csv(file.path(data_dir, "models_data_clean.csv")) %>% 
  filter(!grepl("Wangf", sample))

# Data intersection handling
study_names <- bigtable$sample %>% intersect(colnames(df_genus))
df_phylum_s <- df_genus[, study_names]
bigtable_s <- bigtable %>% filter(sample %in% study_names)

# Create presence/absence matrix
df_phylum_ag <- df_phylum_s[, bigtable_s$sample] %>% 
  t() %>% 
  as.data.frame() %>% 
  aggregate(by = list(bigtable_s$FertiType), FUN = mean) %>% 
  tibble::column_to_rownames("Group.1") %>% 
  t() %>% 
  as.data.frame()

df_phylum_ag_n <- ifelse(df_phylum_ag > 0, 1, 0) %>% as.data.frame()

# Generate UpSet plot
pdf(file.path(output_dir, "p4_UpSet_genus.pdf"), family = "ArialMT", width = 6, height = 6)
UpSetR::upset(df_phylum_ag_n, 
              order.by = "freq", 
              point.size = 3, 
              line.size = 1,
              mainbar.y.label = "Count of Intersection",
              sets.x.label = "Counts by Type",
              number.angles = 0,
              set_size.show = FALSE,
              text.scale = c(1.5, 1.3, 1.3, 1.3, 1.5, 1.5))
dev.off()

# Load taxonomy data
taxa <- read.csv(file.path(data_dir, "taxa_split_unique.csv")) %>% 
  dplyr::select(phylum, genus) %>% 
  unique()

# Helper function for exporting phylum counts
export_phylum_counts <- function(condition, filename) {
  df_phylum_ag_n %>% 
    filter(!!condition) %>% 
    mutate(genus = rownames(.)) %>% 
    left_join(taxa, by = "genus") %>% 
    pull(phylum) %>% 
    table() %>% 
    as.data.frame() %>% 
    write.csv(file.path(output_dir, filename)), row.names = FALSE)
}

# Export unique phylum distributions
export_phylum_counts(expression(IF == 1 & OF == 0 & IFOF == 0 & CK == 0), "p4_UpSet_pie_IF.csv")
export_phylum_counts(expression(IF == 0 & OF == 1 & IFOF == 0 & CK == 0), "p4_UpSet_pie_OF.csv")
export_phylum_counts(expression(IF == 0 & OF == 0 & IFOF == 1 & CK == 0), "p4_UpSet_pie_IFOF.csv")

# Export detailed intersection stats
write_filtered_data <- function(condition, filename) {
  df_phylum_ag_n %>% 
    filter(!!condition) %>% 
    mutate(genus = rownames(.)) %>% 
    left_join(taxa, by = "genus") %>% 
    write.csv(file.path(output_dir, filename)), row.names = FALSE)
}

write_filtered_data(expression(IF == 1 & OF == 0 & IFOF == 0 & CK == 0), "p4_UpSet_pie_stats_IF.csv")
write_filtered_data(expression(IF == 0 & OF == 1 & IFOF == 0 & CK == 0), "p4_UpSet_pie_stats_OF.csv")
write_filtered_data(expression(IF == 0 & OF == 0 & IFOF == 1 & CK == 0), "p4_UpSet_pie_stats_IFOF.csv")
