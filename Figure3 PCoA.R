# PCoA
# Read OTU abundance table
bac.otu<- read.delim("fungi.txt",row.names = 1)
# Read environmental factors/group information table
env.df<-read.csv("Samples group category.csv", row.names = 1)
# Filter samples in the OTU table to match those in the group information table
bac.otu<-bac.otu[,colnames(bac.otu)%in%env.df$Plots]
# Transpose OTU table (samples as rows, OTUs as columns)
bac.otu<-t(bac.otu)
# Sort samples to ensure consistency with the group information table order
bac.otu<-bac.otu[order(row.names(bac.otu)),]
# Sort group information table
env.df<-env.df[order(env.df$Plots),]
# Add group information 'site' to the OTU table
bac.otu1<-data.frame(bac.otu,site=env.df$site)
# Aggregate and calculate mean abundance by 'site'
bac.otu2<-aggregate(.~site,bac.otu1,mean)
# Set row names to 'site'
row.names(bac.otu2)<-bac.otu2$site
# Delete unnecessary column from env.df (assuming it's the Plots column)
env.df<-env.df[,-2]
# Ensure uniqueness of 'site' in env.df
env.df<-unique(env.df)
# Merge the mean abundance OTU table and the group information table
bac.otu3<-merge(bac.otu2,env.df,by="site")
# Set row names to 'site'
row.names(bac.otu3)<-bac.otu3$site
# Delete redundant column ('site' column)
bac.otu3<-bac.otu3[,-1]
# Delete extra column (assuming it's an OTU data column)
bac.otu3<-bac.otu3[,-12722]
# Load vegan package
library(vegan)
# Calculate Bray-Curtis distance matrix
dist<-vegdist(bac.otu3,method = "bray")
# for function use method = "euclidean"
### Distance to centroid difference (Beta diversity variance analysis/within-group dispersion analysis)
# Calculate distance to group centroid
mod<-betadisper(dist,env.df$Years)
pcoa <- capscale(dist ~ 1)
eig <- pcoa$CA$eig
eig_perc <- eig / sum(eig) * 100
eig_perc[1:2]  # PC1 and PC2 percentage

# ANOVA for distances
anova(mod)
# Tukey HSD post hoc multiple comparison for distances
TukeyHSD(mod)
# Permutation test for distances
p_dis<-permutest(mod, pairwise = TRUE, permutations = 999)
### Community difference (PERMANOVA Permutational Multivariate Analysis of Variance)
# Load pairwiseAdonis package
# Install devtools package
library(devtools)
# Install pairwiseAdonis from GitHub
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# Perform overall PERMANOVA analysis
adonis2(bac.otu3~ Years,data=env.df)
# Perform pairwise PERMANOVA comparison
p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "bray", p.adjust.m= "BH")
# for function use : p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "euclidean", p.adjust.m= "BH")
# Reorder factor levels
mod$group<-factor(mod$group,levels = c('DG','RG1','RG2','RG3','CK'))
# Reorder factor levels for centroid names
row.names(mod[["centroids"]])<-factor(row.names(mod[["centroids"]]),levels = c('DG','RG1','RG2','RG3','CK'))
# Plot Boxplot
boxplot(mod)
# Extract PCoA coordinates and group information
dt<-data.frame(mod[["vectors"]],group=env.df$Years)
# Ensure to use base R's mean function
base::mean
# Calculate mean of PCoA1 and PCoA2 for each group (requires doBy package)
library(doBy)
mean=summaryBy(PCoA1+PCoA2~ group, dt, FUN = base::mean)
# Merge means back to the original PCoA coordinate dataframe
dt1<- merge(dt, mean, by = 'group')
# Reorder factor levels
dt1$group<-factor(dt1$group,levels = c('DG','RG1','RG2','RG3','CK'))
# Rename mean columns
colnames(dt1)[57:58]<-c("PCoA1.mean","PCoA2.mean")


# Load plotting and helper packages
library(ggplot2)
library(patchwork)
library(doBy)
library(dplyr)
library(colorspace)   # Used to darken colors

# Define theme (enlarge font + padding around edges)
mytheme <- theme_bw(base_size = 28) + 
  theme(
    legend.title = element_blank(),
    legend.position = 'right',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 22, color = "black"),
    plot.margin = margin(20,20,20,20) # Padding around edges
  )

# Replace group labels
env.df$Years <- factor(env.df$Years,
                       levels = c("DG","RG1","RG2","RG3","CK"),
                       labels = c("Degraded","Early","Middle","Late","Intact"))

# Define color palette (original color + darker stroke)
cols <- c("#be564a","#f9b483","#576289","#6c958f","#3e9fd1")
cols_darker <- darken(cols, amount = 0.4)  # Darken by 40%

# PCoA Plot ---------------------------------------------------------
p1 <- ggplot(data = dt1, aes(x = PCoA1, y = PCoA2, fill = group)) +
  geom_point(shape = 21, size = 5, stroke = 0.8,
             color = cols_darker[as.numeric(factor(dt1$group))]) +  # Darker stroke
  geom_polygon(data = dt1 %>% group_by(group) %>% 
                 slice(chull(PCoA1, PCoA2)),
               aes(x = PCoA1, y = PCoA2, group = group, colour = group),
               fill = NA, size = 1, linetype = "solid") +
  geom_point(aes(x = PCoA1.mean, y = PCoA2.mean, fill = group),
             shape = 24, size = 6, stroke = 0.9,
             color = cols_darker[as.numeric(factor(dt1$group))]) + # Darker stroke for mean point
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = "PC1 (11.08%)", y = "PC2 (7.78%)") +
  mytheme

# Distance to centroid -------------------------------------------------
distence <- data.frame(dis = mod[["distances"]], group = env.df$Years)

p2 <- ggplot(distence, aes(x = group, y = dis, fill = group)) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.2),
               width = 0, size = 0.6) +
  geom_boxplot(color = "black", size = 0.6, outlier.shape = NA) +
  geom_jitter(aes(fill = group),
              shape = 21, size = 3.8, stroke = 0.6,
              color = cols_darker[as.numeric(distence$group)], # Darker stroke
              width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = cols) +
  labs(y = "Distance to centroid", x = NULL) +
  mytheme +
  theme(legend.position = "none")

# Arrange plots side by side using patchwork
combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1))

# Save to PDF
ggsave("PCoA_and_Distance_combined fungi.pdf", combined, width = 16, height = 6)