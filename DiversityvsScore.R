# Step 1: Create a mapping from Experiment to Trial FUNGI
merged.df.fungi$Trial_numeric <- as.integer(sub("T", "", merged.df.fungi$Trial))  # Remove 'T' and convert to integer

# Step 2: Merge based on matching Num (ID) and the mapped Trial_numeric (Experiment)
merged_with_score.trial <- merge(merged.df.fungi, COMBO_combined_data, by.x = c("Num", "Trial_numeric"), by.y = c("ID", "Experiment"), all.x = TRUE)


merged_with_score.trial$sample.id <- NULL


working <- merged_with_score.trial

rm(metadataF, merged_with_score, merged.df.fungi, merged_with_score.trial)
rm(Next.F, iNextEst_Fungi_data)

working$inoc.x<- NULL
working$LCL<- NULL
working$UCL<- NULL
working$is.control<- NULL
working$Block <- NULL
working$Rep <- NULL
working$plate <- NULL
working$col <- NULL
working$con <- NULL
working$Treatment <- NULL
working$T.T <- NULL
working$Trial <- NULL
working$Merge.ID <- NULL
working$WF <- NULL
working$Inoc.y <- NULL
working$row <- NULL
working$sam <- NULL
working$Type <- NULL

#Richnesss Correlation 
Richness <- working[working$Diversity == "Species richness", ]
Richness$is.pythium <- as.factor(Richness$is.pythium)
ggplot(Richness, aes(x = Score, y = Observed, color = is.pythium)) +  
  geom_point() +  
  geom_smooth(method = "lm", aes(group = is.pythium), se = FALSE) +  
  labs(
    title = "Linear Relationship between Observed Species Richness and Disease Severity Score",
    x = "Disease Score",
    y = "Observed Richness, q = 0"
  ) +
  theme_minimal() +
  facet_wrap(~Substrate)

#Testing Linear, Qud, and Cubic
ggplot(Shannon, aes(x = Score, y = Observed)) +  
  geom_point() +  # Scatter plot points
  #geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = TRUE) +  # Linear fit
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", se = TRUE) +  # Quadratic fit
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "green", se = TRUE) +  # Cubic fit
  labs(
    title = "Relationship between Shannon Entropy and Disease Score Varied by Substrate",
    x = "Disease Score",
    y = "Shanonn Entropy, q = 1"
  ) +
  theme_minimal() +
  facet_wrap(~Substrate)


ggplot(Richness, aes(x = Observed, y = Score, color = is.pythium)) +  
  geom_point() +  
  geom_smooth(
    data = Richness[Richness$Substrate == "Peatlite", ],  # Apply regression only to Peatlite
    method = "lm",
    aes(group = is.pythium, color = is.pythium),
    se = FALSE
  ) +  
  labs(
    title = "Linear Relationship between Disease Severity Score and Observed Species Richness",
    x = "Observed Richness, q = 0",  # Now the predictor (independent variable)
    y = "Disease Score"  # Now the dependent variable
  ) +
  theme_minimal() +
  facet_wrap(~Substrate)

#How to test significance of relationships 

# Define the function to run the analysis
run_model_analysis <- function(substrate_level, data) {
  substrate_data <- subset(data, Substrate == substrate_level)
  
  # Linear model (Observed as the dependent variable, Score as the independent variable)
  linear_model <- lm(Observed ~ Score, data = substrate_data)
  
  # Quadratic model (Observed as the dependent variable, Score as the independent variable)
  quadratic_model <- lm(Observed ~ poly(Score, 2), data = substrate_data)
  
  # Cubic model (Observed as the dependent variable, Score as the independent variable)
  cubic_model <- lm(Observed ~ poly(Score, 3), data = substrate_data)
  
  # Compare models using ANOVA
  model_comparison <- anova(linear_model, quadratic_model, cubic_model)
  
  # Summarize each model
  model_summaries <- list(
    linear_summary = summary(linear_model),
    quadratic_summary = summary(quadratic_model),
    cubic_summary = summary(cubic_model)
  )
  
  # Return results
  return(list(
    model_comparison = model_comparison,
    model_summaries = model_summaries
  ))
}


result <- run_model_analysis("ExtrudedFG", Shannon)

# Access the model comparison (ANOVA results)
print(result)
print(result$model_comparison)
print(result$model_summaries$linear_summary)
print(result$model_summaries$quadratic_summary)
print(result$model_summaries$cubic_summary)


Shannon$Substrate <- str_replace(Shannon$Substrate, "Extruded FG", "Extruded")











#Shannon 
Shannon <- working[working$Diversity == "Shannon diversity", ]

ggplot(Shannon, aes(x = Score, y = Observed, color = Substrate)) +  
  geom_point() +  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +  
  labs(
    title = "Shannon Entropy vs. Disease Score by Substrate",
    x = "Disease Score",
    y = "Shanonn Entropy, q = 1"
  ) +
  theme_minimal() +
  facet_wrap(~Substrate)


#different lm


#Simpson 
Simpson <- working[working$Diversity == "Simpson diversity", ]

ggplot(Richness, aes(x = Score, y = Observed, color = Substrate)) +  
  geom_point() +  
  geom_smooth(method = "lm", aes(group = Substrate), se = TRUE) +  
  labs(
    title = "Simpson Diversity vs. Disease Score by Substrate",
    x = "Disease Score",
    y = "Simpson diversity, q = 2"
  ) +
  theme_minimal() +
  facet_wrap(~Substrate)