#How to test significance of relationships 
NoC_WF_Score_Alpha$Score_Group <- ifelse(NoC_WF_Score_Alpha$Score <= 1, "Low (≤1)", 
                                         ifelse(NoC_WF_Score_Alpha$Score >= 3, "Severe (≥3)", "Moderate (1-3)"))
  
  # Linear model (Observed as the dependent variable, Score as the independent variable)
  model <- lm(Score ~ Shannon, data = NoC_WF_Score_Alpha)
  summary(model)
  
  # Quadratic model (Observed as the dependent variable, Score as the independent variable)
  quadratic_model <- lm(Score~ poly(Shannon, 2), data = NoC_WF_Score_Alpha)
  
  # Cubic model (Observed as the dependent variable, Score as the independent variable)
  cubic_model <- lm(Score ~ poly(Shannon, 3), data = NoC_WF_Score_Alpha)
  summary(cubic_model)
  # Compare models using ANOVA
  model_comparison <- anova(model, quadratic_model)
  
  # Summarize each model
  model_summaries <- list(
    linear_summary = summary(model),
    quadratic_summary = summary(quadratic_model)
  )
  
  # Return results
  return(list(
    model_comparison = model_comparison,
    model_summaries = model_summaries
  ))
}

print(model_comparison)
print(model_summaries)


result <- run_model_analysis("ExtrudedFG", Shannon)

# Access the model comparison (ANOVA results)
print(result)
print(result$model_comparison)
print(result$model_summaries$linear_summary)
print(result$model_summaries$quadratic_summary)
print(result$model_summaries$cubic_summary)


Shannon$Substrate <- str_replace(Shannon$Substrate, "Extruded FG", "Extruded")


D <- ggplot(NoC_WF_Score_Alpha, aes(x = Shannon, y = Score, color = Score_Group)) +
  geom_point(aes(shape = Inoc), size = 3, alpha = 0.7) +
  geom_smooth(data = subset(NoC_WF_Score_Alpha, Substrate == "Peatlite"), 
              method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = c("Low (≤1)" = "#009E73", 
                                "Moderate (1-3)" = "#E69F00", 
                                "Severe (≥3)" = "#D55E00")) +
  labs(
    x = "Shannon Index Score",
    y = "Disease Score"
  ) +  
  theme_minimal(base_size = 14) +
  facet_wrap(~ Substrate) +
  theme(legend.position = "none")


  "Disc-refined GF" = "Disc-refined ForestGold®",
  "Extruded FG" = "Extruded Green-Fibre®",
  "Hammermilled-PTS" = "Hammer-milled Pine Tree Substrate (PTS)",
  "Peatlite" = "Peatlite"
)


# Fit the linear model to the subset of data for Peatlite
lm_model <- lm(Score ~ Shannon, data = subset(NoC_WF_Score_Alpha_Inoculated, Substrate == "Peatlite"))

# Perform ANOVA on the model
anova_result <- anova(lm_model)
print(anova_result)

# Extract R-squared value
r_squared <- summary(lm_model)$r.squared
print(r_squared)
# Create the plot


H<-ggplot(NoC_WF_Score_Alpha, aes(x = Simpson, y = Score, color = Score_Group)) +
  geom_point(aes(shape = Inoc), size = 3, alpha = 0.7) + 
  scale_color_manual(values = c("Low (≤1)" = "#009E73", 
                                "Moderate (1-3)" = "#E69F00", 
                                "Severe (≥3)" = "#D55E00")) +
  labs(
       x = "Simpson DIversity Index (1-D)",
       y = "Disease Score", 
       color = "Disease Score",
       shape = "Inoculation Status") +  # Fixed the syntax issue here
  theme_minimal(base_size = 14) +
  facet_wrap(~ Substrate) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )




ibrary(gridExtra)
library(grid)

grid.arrange(D, H, ncol = 2, top = textGrob("Relationship Between Fungal Community Alpha Diversity Metrics and Disease Score was Similar Between Substrates", 
                                               gp = gpar(fontface = "bold", fontsize = 14)))




library(dplyr)

# Add a column to identify each dataset
NoC_WF_Score_Alpha_Inoculated <- NoC_WF_Score_Alpha_Inoculated %>% 
  mutate(Type = "Fungi")

Bacteria_Inoc <- Bacteria_Inoc %>% 
  mutate(Type = "Bacteria")

Bacteria_Inoc$Inoc <- as.factor(Bacteria_Inoc$Inoc)
NoC_WF_Score_Alpha_Inoculated$Inoc <- as.factor(NoC_WF_Score_Alpha_Inoculated$Inoc)
# Combine datasets
combined_df <- bind_rows(NoC_WF_Score_Alpha_Inoculated, Bacteria_Inoc)

# Plot
combined_plot <- ggplot(combined_df, aes(x = Shannon, y = Score, fill = Substrate, shape = Type)) +
  geom_point(size = 4, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("Fungi" = 24, "Bacteria" = 21)) +  # triangle for fungi, circle for bacteria
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("")

combined_plot
ggsave("Dis_Score_Fungi_Bac.png", plot = combined_plot, width =9, height = 6, dpi = 800)