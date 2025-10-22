#How to make the ps.nobins

#Bacterital WF Alpha Diversity 
sample_sums <- sample_sums(psB)
low_read_samples <- sample_sums[sample_sums < 450]
low_read_samples

psB.450 <- prune_samples(sample_sums(psB) >= 450, psB)
psB.450 <- prune_taxa(taxa_sums(psB.450) > 0, psB.450)

psB.450 

#Script works with psB.450 to build ps without bin samples 

#Pull sample metadata 
B450_ALLDATA <- sample_data(psB.450)
View(B450_ALLDATA)

#Prune Bin Samples, bin = samples from substrate bin stocks 
ps.nobin <- prune_samples(!B450_ALLDATA$is.control,psB.450)

Bacteria_Full_Data <- data.frame(Bacteria_Full_Data)
rownames(Bacteria_Full_Data) <- Bacteria_Full_Data$Merge.ID
all(sample_names(ps.nobin) %in% rownames(Bacteria_Full_Data))
all(rownames(Bacteria_Full_Data) %in% sample_names(ps.nobin))

sample_data(ps.nobin) <- Bacteria_Full_Data
head(sample_names(ps.nobin))
head(rownames(Bacteria_Full_Data))

NOBIN_DATA<- sample_data(ps.nobin) # n = 119
NOBIN_DATA$Substrate <- as.factor(NOBIN_DATA$Substrate)
NOBIN_DATA$Inoc <- as.factor(NOBIN_DATA$Inoc)
NOBIN_DATA$Treatment <- as.factor(NOBIN_DATA$Treatment)
levels(NOBIN_DATA$Substrate) <- c(
  "Disc-refined GF" = "Disc-refined ForestGold®",
  "Extruded FG" = "Extruded Green-Fibre®",
  "Hammermilled-PTS" = "Hammer-milled Pine Tree Substrate (PTS)",
  "Peatlite" = "Peatlite"
)


NOBIN_DATA$Inoc <- factor(NOBIN_DATA$Inoc, levels = c("0", "1"), labels = c("Water Control", "Inoculated"))
View(NOBIN_DATA)

sample_data(ps.nobin) <- NOBIN_DATA
OTU_nobin <- otu_table(ps.nobin)
which(colSums(OTU_nobin)==0)

#Only removes sample data, need to remove ASVS that are now 0 total counts 
ps.nobin<-prune_taxa(taxa_sums(ps.nobin) > 0, ps.nobin)
OTU_nobin <- otu_table(ps.nobin)
which(colSums(OTU_nobin)==0)
saveRDS(psB, "~/Desktop/psB.RDS")
saveRDS(ps.nobin,"~/Desktop/ps.nobinB.RDS")


#Creating the data frame was done in two parts: 

#Getting the alpha-diversity measurement data
diversity_metrics <- estimate_richness(ps.nobin, measures = c("Shannon", "Simpson"))
NOBIN_DATA$Shannon <- diversity_metrics$Shannon
NOBIN_DATA$Simpson <- diversity_metrics$Simpson
NOBIN_DATA<-data.frame(sample_data(NOBIN_DATA))

###Getting block column and disease score column added to this data frame was done by sourcing them
#from meta data for ALL samples (including bin and fungal samples)

#Example: 
#NOBIN_DATA <- NOBIN_DATA %>%
#mutate(Block = NoC_WF_Score_Alpha$Block[match(sam, NoC_WF_Score_Alpha$samplenumber)])



#Load in Bacteria_Full_Data and save it as NOBIN_DATA to continue the analysis, 
#this sheet now has the alpha diversity metrics and the block and disease score

NOBIN_DATA <- Bacteria_Full_Data #This is the complete metadata for bacterial samples,
#no bin samples, filtering complete, all columns present 



#Subset by Inoc factor to get rid of interaction, create dis object (matrix)
Inoc_ps_tr_clr <- phyloseq::subset_samples(ps_tr_clr, Inoc =="Inoculated")
DATA.Inoc <- data.frame(sample_data(Inoc_ps_tr_clr))
Inoc.matrix <- phyloseq::distance(Inoc_ps_tr_clr, method ="euclidean")


W_ps_tr_clr <- phyloseq::subset_samples(ps_tr_clr, Inoc =="Water Control")
DATA.W <- data.frame(sample_data(W_ps_tr_clr))
W.matrix <- phyloseq::distance(W_ps_tr_clr, method ="euclidean")



#Unused
# Extract treatment labels from Meta
treatment_labels <- NOBIN_DATA$Treatment
names(treatment_labels) <- rownames(NOBIN_DATA)  # match sample IDs

# Replace row and column names in the matrix
rownames(matrix_nobin) <- treatment_labels[rownames(matrix_nobin)]
colnames(matrix_nobin) <- treatment_labels[colnames(matrix_nobin)]

# Melt the matrix for ggplot
df_melt <- melt(matrix_nobin)
# Optional: keep treatment label order consistent
df_melt$Var1 <- factor(df_melt$Var1, levels = unique(treatment_labels))
df_melt$Var2 <- factor(df_melt$Var2, levels = unique(treatment_labels))

# Plot
ggplot(df_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "#CDE7FF", high = "#D50000" ) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####See if factors associate with  cluster of dendrogram

#create hierarchal clustering object
hc <- hclust(ps_dist_matrix_clr_euc, method = "average")
clusters <- cutree(hc, k = 4)
table(clusters)  # check cluster sizes

metadata <- data.frame(sample_data(ps_tr_clr))  # This is the correct source
metadata$Cluster <- factor(clusters[rownames(metadata)])

# For Substrate
table_substrate <- table(metadata$Cluster, metadata$Substrate)
chisq.test(table_substrate)

# For Inoc
table_inoc <- table(metadata$Cluster, metadata$Inoc)
chisq.test(table_inoc)

# For Treatment
table_treatment <- table(metadata$Cluster, metadata$Treatment)
chisq.test(table_treatment)