library(dplyr)
# Read relative abundance data, replace path for privacy
abun_family <- read.csv("path/to/relative_abundance_family.csv", row.names = 1)
abun_family_f <- abun_family

group_names = c("scale.SOC", "scale.TN", "scale.pH", "scale.duration")
ferti_names = c("IF", "IFOF", "OF")

# Loop through fertilizer and group names
for (ferti_name in ferti_names) {
  for (group_name in group_names) {
    # Read and filter model data, replace path for privacy
    bigtable <- read.csv("path/to/models_data_clean.csv") %>% filter(FertiType == ferti_name) 
    study_names = bigtable$sample %>% intersect(colnames(abun_family_f))
    ARG_beta <- abun_family_f[, study_names] %>% t %>% as.data.frame()
    bigtable_s <- bigtable %>% filter(sample %in% study_names)
    
    library(vegan)
    # Calculate Bray-Curtis dissimilarity
    df_dist <- vegdist(ARG_beta, method = "bray", binary = F)
    # Perform PCoA
    ARG_pcoa <- cmdscale(df_dist, k = 3, eig = T)
    ARG_pcoa_points <- as.data.frame(ARG_pcoa$points)
    sum_eig <- sum(ARG_pcoa$eig)
    eig_percent <- round(ARG_pcoa$eig / sum_eig * 100, 1)
    colnames(ARG_pcoa_points) <- paste0("PCoA", 1:3)
    ARG_pcoa_result <- ARG_pcoa_points %>% mutate(sample = rownames(.)) %>% inner_join(bigtable_s, by = "sample")
    
    library(ggplot2)
    library(ggthemes)
    ARG_pcoa_result_g <- ARG_pcoa_result %>% mutate(Group = factor(get(group_name)))
    
    # Perform PERMANOVA
    otu.adonis = adonis2(df_dist ~ Group, data = ARG_pcoa_result_g, distance = "bray", permutations = 9999)
    
    p <- ggplot(ARG_pcoa_result_g, aes(x = PCoA1, y = PCoA2, color = get(group_name))) +
      geom_point(size = 4) +
      stat_ellipse(level = 0.6) +
      labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
           y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
           color = group_name) +
      theme_base() +
      theme(legend.position = c(0.01, 0.01),
            legend.justification = c(0.01, 0.01),
            legend.box.background = element_rect(fill = "white",
                                                 colour = "black",
                                                 size = 0.5),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 10.5),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 11.5)) +
      scale_color_brewer(palette = "Set2")
    # Save plot, replace path for privacy
    ggsave(p, width = 6, height = 6, filename = paste0("path/to/save/nonote_p2_", ferti_name, "_", group_name, ".pdf"), family = "ArialMT")
  }
}

library(ape)
library(ggplot2)
library(grid)
library(dplyr)
library(multcomp)
library(patchwork)
library(agricolae)

group_name = "FertiType"
ARG_pcoa_result_g <- ARG_pcoa_result %>% mutate(Group = factor(get(group_name)))

# Calculate max values for annotations
yd1 <- ARG_pcoa_result_g %>% group_by(Group) %>% summarise(Max = max(PCoA1))
yd2 <- ARG_pcoa_result_g %>% group_by(Group) %>% summarise(Max = max(PCoA2))
yd1$Max <- yd1$Max + max(yd1$Max) * 0.1
yd2$Max <- yd2$Max + max(yd2$Max) * 0.1

# ANOVA and post - hoc tests for PCoA1
fit1 <- aov(PCoA1 ~ Group, data = ARG_pcoa_result_g)
tuk1 <- glht(fit1, linfct = mcp(Group = "Tukey"))
res1 <- cld(tuk1, alpha = 0.05, decreasing = T)
# ANOVA and post - hoc tests for PCoA2
fit2 <- aov(PCoA2 ~ Group, data = ARG_pcoa_result_g)
tuk2 <- glht(fit2, linfct = mcp(Group = "Tukey"))
res2 <- cld(tuk2, alpha = 0.05, decreasing = T)
lsd1 = LSD.test(fit1, "Group", p.adj = "BH")
lsd2 = LSD.test(fit2, "Group", p.adj = "BH")

test <- data.frame(PCoA1 = res1$mcletters$Letters, PCoA2 = res2$mcletters$Letters,
                   yd1 = yd1$Max, yd2 = yd2$Max, Group = yd1$Group) %>% 
  left_join(lsd2$groups %>% dplyr::select(PCoA_BH2 = groups) %>% mutate(Group = rownames(.)), by = "Group") %>% 
  left_join(lsd1$groups %>% dplyr::select(PCoA_BH1 = groups) %>% mutate(Group = rownames(.)), by = "Group")

otu.adonis = adonis2(df_dist ~ Group, data = ARG_pcoa_result_g, distance = "bray", permutations = 9999)

cbbPalette <- c("#B2182B", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC", "#66CC99", "#99999", "#ADD1E5")

point_size = 4
point_alpha = 0.9

# Create boxplot for PCoA1
p1 <- ggplot(ARG_pcoa_result_g, aes(Group, PCoA1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test, aes(x = Group, y = yd1, label = PCoA_BH1),
            size = 5, color = "black", fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.ticks.length = unit(0.4, "lines"), 
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_blank(),
        legend.position = "none")

# Create boxplot for PCoA2
p3 <- ggplot(ARG_pcoa_result_g, aes(Group, PCoA2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test, aes(x = Group, y = yd2, label = PCoA_BH2),
            size = 5, color = "black", fontface = "bold") +
  scale_fill_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.ticks.length = unit(0.4, "lines"), 
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 10, angle = 45,
                                 vjust = 1, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "none")

# Create scatter plot for PCoA1 and PCoA2
p2 <- ggplot(ARG_pcoa_result_g, aes(PCoA1, PCoA2)) +
  geom_point(aes(fill = Group), size = point_size, pch = 21, alpha = point_alpha, stroke = 0.5) +
  scale_fill_manual(values = cbbPalette, name = "Group") +
  labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
       y = paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text = element_text(size = 10)) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(), 
        axis.title = element_text(color = 'black', size = 12),
        axis.ticks.length = unit(0.4, "lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour = 'black', size = 12, vjust = 7, face = "bold"),
        axis.title.y = element_text(colour = 'black', size = 12, vjust = -2, face = "bold"),
        axis.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.position = c(0.63, 0.15),
        legend.background = element_rect(colour = "black")) +
  guides(fill = guide_legend(ncol = 1))

# Create text plot for PERMANOVA results
p4 <- ggplot(ARG_pcoa_result_g, aes(PCoA1, PCoA2)) +
  geom_text(aes(x = -0.5, y = 0.6, label = paste("PERMANOVA:",  "\nR2 = ", round(otu.adonis[[3]][1], 3),  "\np-value < ", otu.adonis[[5]][1], "\nF-value = ", round(otu.adonis[[4]][1], 1), sep = "")),
            size = 3) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# Combine plots
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
print(p5)

# Save combined plot and data, replace paths for privacy
ggsave(p5, width = 8, height = 8, filename = "path/to/save/figure.pdf", family = "ArialMT")
write.csv(ARG_pcoa_result, "path/to/save/table.csv")
