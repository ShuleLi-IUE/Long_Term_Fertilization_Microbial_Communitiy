# Note to users: 
# 1. Replace 'path/to/your/data' and 'path/to/your/output' with your actual directories
# 2. Ensure required packages are installed
# 3. Set appropriate parameters (ferti_name, topk) for your analysis

# Load required libraries
library(randomForest)
library(dplyr)

# Initialize parameters
ferti_name = "IF"

# Load data (replace paths with your own data directory)
data_dir <- "path/to/your/data"  # Set your data directory here
df_rf_genus <- read.csv(file.path(data_dir, "relative_abundance_genus.csv"), row.names = 1)
bigtable <- read.csv(file.path(data_dir, "models_data_clean.csv")) %>% 
  filter(FertiType == ferti_name | FertiType == "CK")

# Data preprocessing
study_names <- bigtable$sample %>% intersect(colnames(df_rf_genus))
df_rf_genus_s <- df_rf_genus[, study_names] %>% filter(rowMeans(.) > 0.0001)
bigtable_s <- bigtable %>% filter(sample %in% study_names)

# Prepare data for machine learning
otu_rf <- df_rf_genus_s %>% t %>% as.data.frame
otu_rf_g <- otu_rf %>% 
  tibble::rownames_to_column("sample") %>% 
  left_join(bigtable_s %>% dplyr::select(sample, FertiType), by = "sample")

# Export processed data
output_dir <- "path/to/your/output"  # Set your output directory here
otu_rf_g %>% janitor::clean_names() %>% 
  write.csv(file.path(output_dir, paste0("df_for_ml_", ferti_name, ".csv")), row.names = F)

# 10-fold cross-validation
set.seed(123)
otu.cv <- rfcv(otu_rf, factor(bigtable_s$FertiType), 
               scale = NULL, cv.fold = 10, step = 2)

# Visualization of cross-validation results
ferti_name = "IFOF"
otu_train.cv <- data.frame(
  otus = as.integer(names(otu.cv$error.cv)),
  value = as.numeric(otu.cv$error.cv)
)

topk = 112
p <- ggplot(otu_train.cv, aes(otus, value)) +
  geom_line(color = "black", alpha = 0.7, size = 0.6) +
  geom_point(color = "brown3", alpha = 0.5, size = 3) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text = element_text(size = 10, color = "black")) +
  labs(x = 'Number of variables', y = '10-fold cross-validation error') +
  geom_vline(xintercept = topk, linetype = "longdash", color = "brown3")

# Model training and feature importance
otu_rf_g <- read.csv(file.path(output_dir, paste0("df_for_ml_", ferti_name, ".csv")))
otu_rf_clean <- otu_rf_g %>% dplyr::select(-c("sample", "uncultured")) %>% janitor::clean_names()

set.seed(123)
model.otu_rf <- randomForest(as.factor(ferti_type)~., data = otu_rf_clean, 
                            importance = TRUE, ntree = 500)

# Save important features
importance_otu <- data.frame(model.otu_rf$importance)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
otu_select <- rownames(importance_otu)[1:topk]

# Permutation importance analysis
library(rfPermute)
otu_rf_topk <- otu_rf_clean[, c(otu_select, "ferti_type")]
otu_rfP_topk <- rfPermute(as.factor(ferti_type)~., data = otu_rf_topk, 
                         importance = TRUE, ntree = 500)

# Visualization of results
importance_otu.scale <- data.frame(importance(otu_rfP_topk, scale = TRUE), 
                                  check.names = FALSE)
# Process feature importance results
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$MeanDecreaseAccuracy, decreasing = TRUE), ]
write.csv(importance_otu.scale, 
          file = file.path(output_dir, paste0("p5_", ferti_name, "_rf_importance_selected.csv")))

# Prepare visualization data
importance_otu.scale$OTU_name <- rownames(importance_otu.scale)
importance_otu.scale$OTU_name <- factor(importance_otu.scale$OTU_name, 
                                       levels = importance_otu.scale$OTU_name)

# Calculate p-values from permutation test
importance_otu.scale.pval <- (otu_rfP_topk$pval)[, , 2]

# Prepare abundance data for visualization
library(reshape2)
topk_abun_m <- melt(otu_rf_topk, id.vars = "ferti_type", 
                   variable.name = "genus_name", value.name = "abundance") %>% 
  aggregate(. ~ ferti_type + genus_name, FUN = "mean") %>% 
  mutate(abundance = 100 * abundance)

# Sort data by feature importance
topk_abun_m_sorted <- topk_abun_m %>% 
  arrange(match(genus_name, rownames(importance_otu.scale)),
          match(ferti_type, c("CK", "IF")))

# Prepare annotation positions
text_location <- topk_abun_m_sorted %>% 
  group_by(genus_name) %>% 
  summarise(abundance = sum(abundance)) %>% 
  arrange(match(genus_name, rownames(importance_otu.scale)))

# Determine dominant treatment group
larger_group <- topk_abun_m_sorted %>%
  group_by(genus_name) %>%
  mutate(larger_ferti_type = ifelse(abundance == max(abundance), ferti_type, NA)) %>% 
  pull(larger_ferti_type) %>% 
  na.omit() %>% 
  as.vector()

# Create stacked bar plot
p1 <- ggplot(topk_abun_m_sorted, aes(x = genus_name, y = abundance)) +
  geom_bar(aes(fill = ferti_type), stat = "identity", position = "stack", width = 0.9) +
  theme_classic() +
  labs(x = "", y = "Relative abundance (%)") +
  scale_fill_manual(values = c("CK" = "#c9282d", "IFOF" = "#346b95")) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_fixed(ratio = 1.9)

# Save upper plot
ggsave(p1, filename = file.path(output_dir, paste0("p5_", ferti_name, "_rf_fig_upper.pdf")),
       height = 3)

# Prepare taxonomy data (assuming taxa_clean exists)
plot_df <- importance_otu.scale %>% 
  left_join(taxa_clean, by = c("OTU_name" = "genus")) %>% 
  mutate(phylum_s = ifelse(phylum %in% names(fill_colors), phylum, "Other"))

# Define color palette
fill_colors <- c(
  "Proteobacteria" = "#FFC300",
  "Acidobacteriota" = "#82E0AA",
  "Actinobacteriota" = "#9B59B6",
  # ... (keep other color definitions)
  "Other" = "#95A5A6"
)

# Create importance plot
p2 <- ggplot(plot_df) +
  geom_col(aes(x = OTU_name, y = -MeanDecreaseAccuracy, fill = phylum_s), 
           width = 0.8, color = NA) +
  labs(y = 'Mean Decrease Accuracy', fill = 'Phylum') +
  theme_classic() +
  scale_fill_manual(values = fill_colors) +
  coord_fixed(ratio = 0.4)

# Combine and save plots
library(patchwork)
combined_plot <- p1 / p2
ggsave(combined_plot, 
       filename = file.path(output_dir, paste0("p5_", ferti_name, "_rf_fig_ALL.pdf")))

