# Initialize fertilizer name and group name
# Set the initial fertilizer type to "IF"
ferti_name = "IF"
# Set the initial group type to "scale.duration"
group_name = "scale.duration"
# Commented out arrays of group names and fertilizer names
# group_names = c("scale.SOC", "scale.TN", "scale.pH", "scale.duration")
# ferti_names = c("IF", "IFOF", "OF")

# Load necessary libraries
# Load the ggplot2 library for data visualization
library(ggplot2)
# Load the ggpubr library for publication-ready plots
library(ggpubr)

# Define a function to generate a ggplot object
# This function takes a data table, a diversity type, and a test type as input
ggf <- function(table, diversity_type, test_type) {
  p <- ggplot(table, aes(Type, get(diversity_type), group = Type, colour = Type)) +
    # Add jittered points to the plot with transparency and size settings
    geom_jitter(alpha = 0.3, size = 2.5) +
    # Add a boxplot to the plot with transparency, size, width, and fill color settings
    geom_boxplot(alpha = 0.5, size = 1.5, width = 0.7, fill = color) + 
    # Set the color scale manually using limits and values
    scale_color_manual(limits = limits_type,
                       values = value) + 
    # Add a point representing the mean value of each group
    stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, colour = "black") +
    # Apply a black and white theme to the plot
    theme_bw() +
    # Set the coordinate system
    coord_cartesian() +
    # Add text annotations to the plot
    geom_text(data = bigtable_annotated,aes(x = Type,y = yd1,label = PCoA_BH1),
              size = 12,color = "black",fontface = "bold", size.unit = "pt") +
    # Apply a base theme with optional settings
    theme_base(
      # base_line_size = 1
    ) + 
    # Set the plot title, x-axis label, and y-axis label
    labs(title = NULL, x = "Duration (yrs)", y = diversity_type) + 
    # Customize the theme settings
    theme(legend.position = "none", # Whether to add a legend
          plot.title = element_text(size = 15,
                                    hjust = 0.5,
                                    colour = "black",
                                    face = "bold"),
          axis.title.y = element_text(size = 12, color = "black", vjust = 1.9,hjust = 0.5, angle = 90),
          axis.title.x = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 12,
                                      color = "black"),
          # face = "bold"),
          legend.text = element_text(size = 10,
                                     color = "black",
                                     # face = "bold"
          ),
          axis.text.x = element_text(size = 10,
                                     color = "black",
                                     # face = "bold",
                                     vjust = 0.5,
                                     hjust = 0.5,
                                     angle = 0),
          axis.text.y = element_text(size = 10,
                                     color = "black",
                                     # face = "bold",
                                     vjust = 0.5,
                                     hjust = 0.5,
                                     angle = 0)
    )
  return(p)
}

# Define arrays of group names and fertilizer names
group_names = c("scale.SOC", "scale.TN", "scale.pH", "scale.duration")
ferti_names = c("IF", "IFOF", "OF")
# Commented out array of index names
# index_names = c("Shannon", "Richness", "InvSimpson", "Chao1")
index_names = c("Pielou")

# Loop through each index name, fertilizer name, and group name
for (index_name in index_names) {
  for (ferti_name in ferti_names) {
    for (group_name in group_names) {
      # Read the data from a CSV file and filter it based on fertilizer type and sample name
      # Replace the file path with a placeholder for privacy
      bigtable <- read.csv("path/to/your/data.csv") %>% filter(FertiType == ferti_name) %>% filter(!grepl("Wangf", sample))
      
      # Modify the scale.duration variable if the group name is "scale.duration"
      if (group_name == "scale.duration") {
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "low", "median", bigtable$scale.duration)
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "very low", "low", bigtable$scale.duration)
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "low", "≤ 15", bigtable$scale.duration)
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "median", "15-25", bigtable$scale.duration)
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "high", "25-35", bigtable$scale.duration)
        bigtable$scale.duration <- ifelse(bigtable$scale.duration == "very high", "> 35", bigtable$scale.duration)
        bigtable$scale.duration <- factor(bigtable$scale.duration, levels = c("≤ 15", "15-25", "25-35", "> 35"))
      }
      
      # Modify the scale.TN variable if the group name is "scale.TN"
      if (group_name == "scale.TN") {
        bigtable$scale.TN <- ifelse(bigtable$scale.TN == "low", "≤ 1.5", bigtable$scale.TN)
        bigtable$scale.TN <- ifelse(bigtable$scale.TN == "median", "1.5-3", bigtable$scale.TN)
        bigtable$scale.TN <- ifelse(bigtable$scale.TN == "high", "> 3", bigtable$scale.TN)
        bigtable$scale.TN <- factor(bigtable$scale.TN, levels = c("≤ 1.5", "1.5-3", "> 3"))
      }
      
      # Modify the scale.pH variable if the group name is "scale.pH"
      if (group_name == "scale.pH") {
        bigtable$scale.pH <- ifelse(bigtable$scale.pH == "acid", "≤ 6.5", bigtable$scale.pH)
        bigtable$scale.pH <- ifelse(bigtable$scale.pH == "neutral", "6.5-7.5", bigtable$scale.pH)
        bigtable$scale.pH <- ifelse(bigtable$scale.pH == "base", "> 7.5", bigtable$scale.pH)
        bigtable$scale.pH <- factor(bigtable$scale.pH, levels = c("≤ 6.5", "6.5-7.5", "> 7.5"))
      }
      
      # Modify the scale.SOC variable if the group name is "scale.SOC"
      if (group_name == "scale.SOC") {
        bigtable$scale.SOC <- ifelse(bigtable$scale.SOC == "low", "≤ 10", bigtable$scale.SOC)
        bigtable$scale.SOC <- ifelse(bigtable$scale.SOC == "median", "10-20", bigtable$scale.SOC)
        bigtable$scale.SOC <- ifelse(bigtable$scale.SOC == "high", "20-30", bigtable$scale.SOC)
        bigtable$scale.SOC <- ifelse(bigtable$scale.SOC == "very high", "> 30", bigtable$scale.SOC)
        bigtable$scale.SOC <- factor(bigtable$scale.SOC, levels = c("≤ 10", "10-20", "20-30", "> 30"))
      }
      
      # Intersect the sample names with the column names of another data frame
      study_names = bigtable$sample %>% intersect(colnames(abun_family_f))
      # Filter the data based on the study names and add a new variable "Type"
      bigtable_s <- bigtable %>% filter(sample %in% study_names) %>% mutate(Type = get(group_name))
      
      # # Commented out code for data transformation
      # df_g <- bigtable_s %>% mutate(Group = get(group_name),
      #                               Value = Shannon) %>% dplyr::select(Group, Value)
      
      # Summarize the data by the "Type" variable
      Type_summary <- bigtable_s %>% group_by(Type) %>% count
      print(Type_summary)
      # Extract the unique "Type" values
      limits_type <- Type_summary %>% dplyr::select(Type) %>% pull
      
      # Calculate the maximum value of the index name for each "Type" and add a small offset
      yd1 <- bigtable_s %>% group_by(Type) %>% summarise(Max = max(get(index_name)))
      yd1$Max <- yd1$Max + max(yd1$Max)*0.03
      
      # Perform an ANOVA test
      fit1 <- aov(get(index_name)~Type,data = bigtable_s)
      # Perform a Tukey's HSD test
      tuk1<-glht(fit1, linfct=mcp(Type="Tukey"))
      # Calculate the compact letter display
      res1 <- cld(tuk1,alpha=0.05, decreasing = T)
      # Perform an LSD test with BH adjustment
      lsd1 = LSD.test(fit1, "Type", p.adj="BH")
      
      # Create a data frame with annotations
      bigtable_annotated <- data.frame(yd1 = yd1$Max,Type = yd1$Type) %>% 
        left_join(lsd1$groups %>% dplyr::select(PCoA_BH1 = groups) %>% mutate(Type=rownames(.)), by="Type")
      
      # Get color palettes
      color <- get_lsl_palette(n = nrow(Type_summary), type = "normal") # Original color
      value <- get_lsl_palette(n = nrow(Type_summary), type = "desat") # Desaturated color
      
      # Generate a plot using the ggf function
      p0 <- ggf(bigtable_s , index_name)
      print(p0)
      
      # Save the plot as a PDF file
      # Replace the file path with a placeholder for privacy
      ggsave(p0, width = 5, height = 4.5,
             filename = paste0("path/to/save/figures/p3_", index_name, "_", ferti_name, "_", group_name, ".pdf"),
             family = "ArialMT")
    }
  }
}
