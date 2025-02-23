# Configure directories (user should set these)
data_dir <- "path/to/your/data"
output_dir <- "path/to/your/output"
fig_subdir <- file.path(output_dir, "p3_stackbar_subplots")

# Create output directory if not exists
dir.create(fig_subdir, showWarnings = FALSE, recursive = TRUE)

# Load and preprocess phylum data
df_phylum <- read.csv(file.path(data_dir, "relative_abundance_phylum.csv"), row.names = 1)
bigtable <- read.csv(file.path(data_dir, "models_data_clean.csv")) 

# Data intersection handling
study_names <- bigtable$sample %>% intersect(colnames(df_phylum))
df_phylum_s <- df_phylum[, study_names]
bigtable_s <- bigtable %>% filter(sample %in% study_names)

# Aggregate top 10 phyla
df_phylum_s.ave <- apply(df_phylum_s, 1, mean) %>% sort(decreasing = TRUE)
names_top10 <- names(df_phylum_s.ave)[1:10]
class_final <- rbind(
  df_phylum_s[names_top10, ],
  Others = apply(df_phylum_s[names_top10, ], 2, function(x) 1 - sum(x))
)

# Visualization parameters
library(colortools)
values <- rev(wheel("skyblue3", num = nrow(class_final)))

# Main visualization function
f_stackbar_project <- function(ferti_name, project, 
                              show_legend = TRUE, 
                              show_label_y = TRUE, 
                              show_breaks = TRUE) {
  # Load and preprocess data
  bigtable <- read.csv(file.path(data_dir, "models_data_clean.csv")) %>% 
    filter(FertiType == ferti_name) %>% 
    filter(!grepl("Wangf", sample))
  
  # Factor level formatting
  if(project == "scale.duration") {
    bigtable <- bigtable %>% mutate(
      scale.duration = case_when(
        scale.duration == "very low" ~ "≤ 15",
        scale.duration == "low" ~ "15-25",
        scale.duration == "median" ~ "25-35",
        scale.duration == "high" ~ "> 35"
      ),
      scale.duration = factor(scale.duration, 
                             levels = c("≤ 15", "15-25", "25-35", "> 35"))
    )
  }
  
  # Similar formatting blocks for other scales...
  # Do what you want

  # Prepare visualization data
  group_data <- bigtable %>% 
    mutate(Type = get(project)) %>% 
    filter(sample %in% colnames(class_final))
  
  # Create stacked bar plot
  library(ggplot2)
  p <- ggplot(class.gg.group, aes(Type, 100 * Abundance, fill = Order)) +
    geom_bar(stat = "summary", fun = mean, position = 'stack') +
    scale_fill_manual(values = values) +
    labs(x = '', y = ifelse(show_label_y, 'Relative Abundance (%)', '')) +
    theme_minimal() +
    coord_flip()
  
  # Save visualization
  ggsave(p, filename = file.path(fig_subdir, 
                                paste0("p3_stackbar_fig_", ferti_name, "_", project, ".pdf")),
       width = 8, height = 6)
  
  # Export aggregated data
  write.csv(df_phylum_ag_sort, 
           file.path(fig_subdir, 
                    paste0("p3_stackbar_data_", ferti_name, "_", project, ".csv")))
}

# Batch processing
group_names <- c("scale.SOC", "scale.TN", "scale.pH", "scale.duration")
ferti_names <- c("IF", "IFOF", "OF")

# Generate all combinations
for (group_name in group_names) {
  for (ferti_name in ferti_names) {
    f_stackbar_project(ferti_name, group_name)
  }
}
