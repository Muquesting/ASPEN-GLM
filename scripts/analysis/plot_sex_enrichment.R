library(ggplot2)
library(dplyr)

# Load results
results <- read.csv("results/analysis/sex_enrichment_only_vs_overlap.csv")

# Create visualization
results$Group <- factor(results$only_group, 
                       levels = c("scDALI_only", "ASPEN_only", "glmmTMB_only", "GAMLSS_only"))

# Add significance markers
results$significance <- ifelse(results$p_value < 0.001, "***",
                              ifelse(results$p_value < 0.01, "**",
                                    ifelse(results$p_value < 0.05, "*", "ns")))

# Plot
p <- ggplot(results, aes(x = Group, y = prop_only_sexsig * 100, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = results$prop_overlap_sexsig[1] * 100, 
             linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = significance), vjust = -0.5, size = 6) +
  geom_text(aes(label = paste0(round(prop_only_sexsig * 100, 1), "%")), 
            vjust = 1.5, color = "white", fontface = "bold") +
  annotate("text", x = 4, y = results$prop_overlap_sexsig[1] * 100 + 2, 
           label = "Overlap baseline", color = "red", hjust = 1) +
  labs(title = "Sex-Significant Gene Enrichment: 'Only' Groups vs 'Overlap'",
       subtitle = "Fisher's Exact Test (* p<0.05, ** p<0.01, *** p<0.001)",
       x = "", y = "% Sex-Significant Genes (p_sex < 0.05)",
       caption = "Red dashed line = proportion in 'Both/Overlap' group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = c("scDALI_only" = "#E74C3C", 
                                "ASPEN_only" = "#3498DB",
                                "glmmTMB_only" = "#2ECC71",
                                "GAMLSS_only" = "#F39C12"))

ggsave("results/analysis/sex_enrichment_barplot.png", p, width = 10, height = 6)
message("Saved: results/analysis/sex_enrichment_barplot.png")
