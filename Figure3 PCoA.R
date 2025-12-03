############################## Fungi PCOA #########################################################

# Read OTU abundance table
bac.otu<- read.delim("fungi.txt",row.names = 1)

# Read environmental factors / grouping information table
env.df<-read.csv("Samples group category.csv", row.names = 1)

# Filter samples in the OTU table to match samples in the grouping table
bac.otu<-bac.otu[,colnames(bac.otu)%in%env.df$Plots]

# Transpose OTU table (samples as rows, OTUs as columns)
bac.otu<-t(bac.otu)

# Sort samples to ensure the same order as in the grouping table
bac.otu<-bac.otu[order(row.names(bac.otu)),]

# Sort grouping information table
env.df<-env.df[order(env.df$Plots),]

# Add grouping information 'site' to the OTU table
bac.otu1<-data.frame(bac.otu,site=env.df$site)

# Aggregate and calculate mean abundance by 'site'
bac.otu2<-aggregate(.~site,bac.otu1,mean)

# Set row names to 'site'
row.names(bac.otu2)<-bac.otu2$site

# Remove redundant column from env.df (assumed to be the Plots column)
env.df<-env.df[,-2]

# Ensure uniqueness of 'site' in env.df
env.df<-unique(env.df)

# Merge averaged OTU table and grouping information table
bac.otu3<-merge(bac.otu2,env.df,by="site")

# Set row names to 'site'
row.names(bac.otu3)<-bac.otu3$site

# Remove redundant column ('site' column)
bac.otu3<-bac.otu3[,-1]

# Remove additional column (assumed to be one OTU column)
bac.otu3<-bac.otu3[,-12722]

# Load vegan package
library(vegan)

# Calculate Bray–Curtis distance matrix
dist<-vegdist(bac.otu3,method = "bray")

# for the function use method = "euclidean"

### Distance to centroid differences (beta diversity variance / within-group dispersion analysis)

# Calculate distances to group centroids
mod<-betadisper(dist,env.df$Years)

pcoa <- capscale(dist ~ 1)

eig <- pcoa$CA$eig

eig_perc <- eig / sum(eig) * 100

eig_perc[1:2]  # Percentages of PC1 and PC2



# Perform ANOVA on distances
anova(mod)

# Perform Tukey HSD post-hoc multiple comparisons on distances
TukeyHSD(mod)

# Perform permutation test on distances
p_dis<-permutest(mod, pairwise = TRUE, permutations = 999)

### Community differences (PERMANOVA permutation-based multivariate ANOVA)

# Load pairwiseAdonis package

# Install devtools package
library(devtools)

# Install pairwiseAdonis from GitHub
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

# Perform overall PERMANOVA analysis
adonis2(bac.otu3~ Years,data=env.df)

# Perform pairwise PERMANOVA comparisons
p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "bray", p.adjust.m= "BH")

# for the function use: p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "euclidean", p.adjust.m= "BH")

# Reorder factor levels
mod$group<-factor(mod$group,levels = c('DG','RG1','RG2','RG3','CK'))

# Reorder factor levels for centroid names
row.names(mod[["centroids"]])<-factor(row.names(mod[["centroids"]]),levels = c('DG','RG1','RG2','RG3','CK'))

# Plot boxplot
boxplot(mod)

# Extract PCoA coordinates and grouping information
dt<-data.frame(mod[["vectors"]],group=env.df$Years)

# Ensure using the base R mean function
base::mean

# Calculate mean values of PCoA1 and PCoA2 for each group (requires doBy package)
library(doBy)

mean=summaryBy(PCoA1+PCoA2~ group, dt, FUN = base::mean)

# Merge mean values back to the original PCoA coordinate data frame
dt1<- merge(dt, mean, by = 'group')

# Reorder factor levels
dt1$group<-factor(dt1$group,levels = c('DG','RG1','RG2','RG3','CK'))

# Rename mean value columns
colnames(dt1)[57:58]<-c("PCoA1.mean","PCoA2.mean")

# Load plotting and helper packages
library(ggplot2)

library(patchwork)

library(doBy)

library(dplyr)

library(colorspace)   # Used to darken colors

# Define theme (larger fonts + outer margins)
mytheme <- theme_bw(base_size = 28) + 

  theme(

    legend.title = element_blank(),

    legend.position = 'right',

    panel.grid.major = element_blank(),

    panel.grid.minor = element_blank(),

    panel.background = element_blank(),

    axis.title = element_text(size = 22, face = "bold"),

    axis.text  = element_text(size = 22, color = "black"),

    plot.margin = margin(20,20,20,20) # Outer margins on all sides

  )


# Replace group labels
env.df$Years <- factor(env.df$Years,

                       levels = c("DG","RG1","RG2","RG3","CK"),

                       labels = c("Degraded","Early","Middle","Late","Intact"))



# Define colors (original fill colors + darker outlines)
cols <- c("#be564a","#f9b483","#576289","#6c958f","#3e9fd1")

cols_darker <- darken(cols, amount = 0.4)  # Darken by 40%

# PCoA plot ---------------------------------------------------------

p1 <- ggplot(data = dt1, aes(x = PCoA1, y = PCoA2, fill = group)) +

  geom_point(shape = 21, size = 5, stroke = 0.8,

             color = cols_darker[as.numeric(factor(dt1$group))]) +  # Dark outline for points

  

  # ✅ Modified part: convex hull border color is transparent, fill color is kept

  geom_polygon(

    data = dt1 %>% group_by(group) %>% slice(chull(PCoA1, PCoA2)),

    aes(x = PCoA1, y = PCoA2, group = group, fill = group),

    color = NA,         # ⬅️ Transparent hull border

    alpha = 0.2,        # ⬅️ Keep semi-transparent filled hull

    size = 1, 

    linetype = "solid"

  ) +

  
  geom_point(aes(x = PCoA1.mean, y = PCoA2.mean, fill = group),

             shape = 24, size = 6, stroke = 0.9,

             color = cols_darker[as.numeric(factor(dt1$group))]) + # Dark outline for group mean points

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

              color = cols_darker[as.numeric(distence$group)], # Dark outline for points

              width = 0.2, alpha = 0.8) +

  scale_fill_manual(values = cols) +

  labs(y = "Distance to centroid", x = NULL) +

  mytheme +

  theme(legend.position = "none")

# Arrange plots side by side using patchwork
combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1))

# Save to PDF
ggsave("PCoA_and_Distance_combined_fungi.pdf", combined, width = 16, height = 6)






########################################### Bacteria PCOA #################################################
# Read OTU abundance table
bac.otu<- read.delim("bacteria.txt",row.names = 1)
# Read environmental factor / grouping information table
env.df<-read.csv("Samples group category.csv", row.names = 1)
# Filter samples in OTU table to match samples in grouping table
bac.otu<-bac.otu[,colnames(bac.otu)%in%env.df$Plots]
# Transpose OTU table (samples as rows, OTUs as columns)
bac.otu<-t(bac.otu)
# Sort samples to ensure the same order as in the grouping table
bac.otu<-bac.otu[order(row.names(bac.otu)),]
# Sort grouping information table
env.df<-env.df[order(env.df$Plots),]
# Add grouping information 'site' to the OTU table
bac.otu1<-data.frame(bac.otu,site=env.df$site)
# Aggregate by 'site' and calculate mean abundance
bac.otu2<-aggregate(.~site,bac.otu1,mean)
# Set row names to 'site'
row.names(bac.otu2)<-bac.otu2$site
# Remove redundant column in env.df (assumed to be the Plots column)
env.df<-env.df[,-2]
# Ensure uniqueness of 'site' in env.df
env.df<-unique(env.df)
# Merge mean OTU table with grouping information table
bac.otu3<-merge(bac.otu2,env.df,by="site")
# Set row names to 'site'
row.names(bac.otu3)<-bac.otu3$site
# Remove redundant column ('site' column)
bac.otu3<-bac.otu3[,-1]
# Remove redundant column ('site' column)
bac.otu3<-bac.otu3[,-1]
# Remove additional column (assumed to be one OTU column)
bac.otu3<-bac.otu3[,-43918]
# Load vegan package
library(vegan)
# Calculate Bray–Curtis distance matrix
dist<-vegdist(bac.otu3,method = "bray")
# for function use method = "euclidean"
### Differences in distance to centroids (beta diversity variance / within-group dispersion analysis)
# Calculate distances to group centroids
mod<-betadisper(dist,env.df$Years)
pcoa <- capscale(dist ~ 1)
eig <- pcoa$CA$eig
eig_perc <- eig / sum(eig) * 100
eig_perc[1:2]  # Percentages of PC1 and PC2


# Perform ANOVA on distances
anova(mod)
# Perform Tukey HSD post-hoc multiple comparisons on distances
TukeyHSD(mod)
# Perform permutation test on distances
p_dis<-permutest(mod, pairwise = TRUE, permutations = 999)
### Community differences (PERMANOVA permutation-based multivariate ANOVA)
# Load pairwiseAdonis package
# Install devtools package
library(devtools)
# Install pairwiseAdonis from GitHub
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# Perform overall PERMANOVA analysis
adonis2(bac.otu3~ Years,data=env.df)
# Perform pairwise PERMANOVA
p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "bray", p.adjust.m= "BH")
# for function use : p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "euclidean", p.adjust.m= "BH")
# Reorder factor levels
mod$group<-factor(mod$group,levels = c('DG','RG1','RG2','RG3','CK'))
# Reorder factor levels for centroid names
row.names(mod[["centroids"]])<-factor(row.names(mod[["centroids"]]),levels = c('DG','RG1','RG2','RG3','CK'))
# Draw boxplot
boxplot(mod)
# Extract PCoA coordinates and grouping information
dt<-data.frame(mod[["vectors"]],group=env.df$Years)
# Ensure using base R mean function
base::mean
# Calculate mean PCoA1 and PCoA2 for each group (requires doBy package)
library(doBy)
mean=summaryBy(PCoA1+PCoA2~ group, dt, FUN = base::mean)
# Merge mean values back into the original PCoA coordinate data frame
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
library(colorspace)   # Used to darken colors

# Define theme (larger fonts + margins on all sides)
mytheme <- theme_bw(base_size = 28) + 
  theme(
    legend.title = element_blank(),
    legend.position = 'right',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 22, color = "black"),
    plot.margin = margin(20,20,20,20) # Margins on all sides
  )

# Replace group labels
env.df$Years <- factor(env.df$Years,
                       levels = c("DG","RG1","RG2","RG3","CK"),
                       labels = c("Degraded","Early","Middle","Late","Intact"))

# Define palette (original colors + darker outlines)
cols <- c("#be564a","#f9b483","#576289","#6c958f","#3e9fd1")
cols_darker <- darken(cols, amount = 0.4)  # Darken by 40%

# PCoA plot ---------------------------------------------------------
p1 <- ggplot(data = dt1, aes(x = PCoA1, y = PCoA2, fill = group)) +
  geom_point(shape = 21, size = 5, stroke = 0.8,
             color = cols_darker[as.numeric(factor(dt1$group))]) +  # Dark outline
  
  # ✅ Modified part: convex hull border color is transparent, fill color retained
  geom_polygon(
    data = dt1 %>% group_by(group) %>% slice(chull(PCoA1, PCoA2)),
    aes(x = PCoA1, y = PCoA2, group = group, fill = group),
    color = NA,         # ⬅️ Transparent hull border
    alpha = 0.2,        # ⬅️ Keep semi-transparent fill
    size = 1, 
    linetype = "solid"
  ) +
  
  geom_point(aes(x = PCoA1.mean, y = PCoA2.mean, fill = group),
             shape = 24, size = 6, stroke = 0.9,
             color = cols_darker[as.numeric(factor(dt1$group))]) + # Dark outline for mean points
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
              color = cols_darker[as.numeric(distence$group)], # Dark outline
              width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = cols) +
  labs(y = "Distance to centroid", x = NULL) +
  mytheme +
  theme(legend.position = "none")

# Arrange left and right using patchwork
combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1))

# Save to PDF
ggsave("PCoA_and_Distance_combined_bacteria.pdf", combined, width = 16, height = 6)




########################################################## Plant PCOA ##################################################
# Read OTU abundance table
bac.otu<- read.delim("plant.txt",row.names = 1)
# Read environmental factor / grouping information table
env.df<-read.csv("Samples group category.csv", row.names = 1)
# Filter samples in the OTU table to match samples in the grouping table
bac.otu<-bac.otu[,colnames(bac.otu)%in%env.df$Plots]
# Transpose OTU table (samples as rows, OTUs as columns)
bac.otu<-t(bac.otu)
# Sort samples to ensure the same order as in the grouping table
bac.otu<-bac.otu[order(row.names(bac.otu)),]
# Sort grouping information table
env.df<-env.df[order(env.df$Plots),]
# Add grouping information 'site' to the OTU table
bac.otu1<-data.frame(bac.otu,site=env.df$site)
# Aggregate by 'site' and calculate mean abundance
bac.otu2<-aggregate(.~site,bac.otu1,mean)
# Set row names to 'site'
row.names(bac.otu2)<-bac.otu2$site
# Remove redundant column in env.df (assumed to be the Plots column)
env.df<-env.df[,-2]
# Ensure uniqueness of 'site' in env.df
env.df<-unique(env.df)
# Merge mean OTU table with the grouping information table
bac.otu3<-merge(bac.otu2,env.df,by="site")
# Set row names to 'site'
row.names(bac.otu3)<-bac.otu3$site
# Remove redundant column ('site' column)
bac.otu3<-bac.otu3[,-1]
# Remove additional column (assumed to be one OTU column)
bac.otu3<-bac.otu3[,-102]
# Load vegan package
library(vegan)
# Calculate Bray–Curtis distance matrix
dist<-vegdist(bac.otu3,method = "bray")
# for function use method = "euclidean"
### Differences in distance to centroids (beta diversity variance analysis / within-group dispersion analysis)
# Calculate distances to group centroids
mod<-betadisper(dist,env.df$Years)
pcoa <- capscale(dist ~ 1)
eig <- pcoa$CA$eig
eig_perc <- eig / sum(eig) * 100
eig_perc[1:2]  # Percentages of PC1 and PC2


# Perform ANOVA on distances
anova(mod)
# Perform Tukey HSD post-hoc multiple comparisons on distances
TukeyHSD(mod)
# Perform permutation test on distances
p_dis<-permutest(mod, pairwise = TRUE, permutations = 999)
### Community differences (PERMANOVA permutation-based multivariate ANOVA)
# Load pairwiseAdonis package
# Install devtools package
library(devtools)
# Install pairwiseAdonis from GitHub
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# Perform overall PERMANOVA analysis
adonis2(bac.otu3~ Years,data=env.df)
# Perform pairwise PERMANOVA
p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "bray", p.adjust.m= "BH")
# for function use : p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "euclidean", p.adjust.m= "BH")
# Reorder factor levels
mod$group<-factor(mod$group,levels = c('DG','RG1','RG2','RG3','CK'))
# Reorder factor levels for centroid names
row.names(mod[["centroids"]])<-factor(row.names(mod[["centroids"]]),levels = c('DG','RG1','RG2','RG3','CK'))
# Draw boxplot
boxplot(mod)
# Extract PCoA coordinates and grouping information
dt<-data.frame(mod[["vectors"]],group=env.df$Years)
# Ensure using base R mean function
base::mean
# Calculate mean PCoA1 and PCoA2 for each group (requires doBy package)
library(doBy)
mean=summaryBy(PCoA1+PCoA2~ group, dt, FUN = base::mean)
# Merge mean values back into the original PCoA coordinate data frame
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
library(colorspace)   # Used to darken colors

# Define theme (larger fonts + margins on all sides)
mytheme <- theme_bw(base_size = 28) + 
  theme(
    legend.title = element_blank(),
    legend.position = 'right',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 22, color = "black"),
    plot.margin = margin(20,20,20,20) # Margins on all sides
  )

# Replace group labels
env.df$Years <- factor(env.df$Years,
                       levels = c("DG","RG1","RG2","RG3","CK"),
                       labels = c("Degraded","Early","Middle","Late","Intact"))

# Define color palette (original colors + darker outlines)
cols <- c("#be564a","#f9b483","#576289","#6c958f","#3e9fd1")
cols_darker <- darken(cols, amount = 0.4)  # Darken by 40%

# PCoA plot ---------------------------------------------------------
p1 <- ggplot(data = dt1, aes(x = PCoA1, y = PCoA2, fill = group)) +
  geom_point(shape = 21, size = 5, stroke = 0.8,
             color = cols_darker[as.numeric(factor(dt1$group))]) +  # Dark outline
  
  # ✅ Modified section: convex hull border color transparent, fill color retained
  geom_polygon(
    data = dt1 %>% group_by(group) %>% slice(chull(PCoA1, PCoA2)),
    aes(x = PCoA1, y = PCoA2, group = group, fill = group),
    color = NA,         # ⬅️ Transparent hull border
    alpha = 0.2,        # ⬅️ Keep semi-transparent fill
    size = 1, 
    linetype = "solid"
  ) +
  
  geom_point(aes(x = PCoA1.mean, y = PCoA2.mean, fill = group),
             shape = 24, size = 6, stroke = 0.9,
             color = cols_darker[as.numeric(factor(dt1$group))]) + # Dark outline for mean points
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
              color = cols_darker[as.numeric(distence$group)], # Dark outline
              width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = cols) +
  labs(y = "Distance to centroid", x = NULL) +
  mytheme +
  theme(legend.position = "none")

# Arrange plots left and right using patchwork
combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1))

# Save to PDF
ggsave("PCoA_and_Distance_combined_plant.pdf", combined, width = 16, height = 6)




####################################################### Functions PCOA ################################################# 
# Read OTU abundance table
bac.otu<- read.delim("function.txt",row.names = 1)
# Read environmental factors / grouping information table
env.df<-read.csv("Samples group category.csv", row.names = 1)
# Filter samples in OTU table to match those in grouping information table
bac.otu<-bac.otu[,colnames(bac.otu)%in%env.df$Plots]
# Transpose OTU table (samples as rows, OTU as columns)
bac.otu<-t(bac.otu)
# Sort samples to ensure consistency with grouping information table
bac.otu<-bac.otu[order(row.names(bac.otu)),]
# Sort grouping information table
env.df<-env.df[order(env.df$Plots),]

# Add grouping information 'site' into OTU table
bac.otu1<-data.frame(bac.otu,site=env.df$site)
# Aggregate by 'site' to calculate mean abundance
bac.otu2<-aggregate(.~site,bac.otu1,mean)
# Set row names as 'site'
row.names(bac.otu2)<-bac.otu2$site
# Remove extra column from env.df (assumed to be Plots column)
env.df<-env.df[,-2]
# Ensure uniqueness of 'site' in env.df
env.df<-unique(env.df)
# Merge averaged OTU table and grouping information table
bac.otu3<-merge(bac.otu2,env.df,by="site")
# Set row names as 'site'
row.names(bac.otu3)<-bac.otu3$site
# Remove extra column ('site' column)
bac.otu3<-bac.otu3[,-1]
# Remove extra column (assumed to be a certain OTU data column)
bac.otu3<-bac.otu3[,-9]
# Load vegan package
library(vegan)
# Calculate Bray-Curtis distance matrix
dist<-vegdist(bac.otu3,method = "bray")
#for function use method = "euclidean"
### Distance to centroid differences (Beta diversity variance analysis / within-group dispersion analysis)
# Calculate distance to group centroids
mod<-betadisper(dist,env.df$Years)
pcoa <- capscale(dist ~ 1)
eig <- pcoa$CA$eig
eig_perc <- eig / sum(eig) * 100
eig_perc[1:2]  # Percentages of PC1 and PC2

# ANOVA on distances
anova(mod)
# Tukey HSD post-hoc test on distances
TukeyHSD(mod)
# Permutation test on distances
p_dis<-permutest(mod, pairwise = TRUE, permutations = 999)
### Community difference (PERMANOVA permutation multivariate analysis of variance)
# Load pairwiseAdonis package
# Install devtools package
library(devtools)
# Install pairwiseAdonis from GitHub
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# Perform overall PERMANOVA analysis
adonis2(bac.otu3~ Years,data=env.df)
# Pairwise PERMANOVA comparisons
p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "bray", p.adjust.m= "BH")
# for function use : p_perma<-pairwiseAdonis::pairwise.adonis(bac.otu3, env.df$Years,sim.method = "euclidean", p.adjust.m= "BH")
# Reorder factor levels
mod$group<-factor(mod$group,levels = c('DG','RG1','RG2','RG3','CK'))
# Reorder centroid factor levels
row.names(mod[["centroids"]])<-factor(row.names(mod[["centroids"]]),levels = c('DG','RG1','RG2','RG3','CK'))
# Plot Boxplot
boxplot(mod)
# Extract PCoA coordinates and grouping information
dt<-data.frame(mod[["vectors"]],group=env.df$Years)
# Ensure using base R mean function
base::mean
# Calculate group means of PCoA1 and PCoA2 (requires doBy package)
library(doBy)
mean=summaryBy(PCoA1+PCoA2~ group, dt, FUN = base::mean)
# Merge means back into original PCoA coordinate dataframe
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
library(colorspace)   # For darkening colors

# Define theme (larger font + plot margins)
mytheme <- theme_bw(base_size = 28) + 
  theme(
    legend.title = element_blank(),
    legend.position = 'right',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 22, color = "black"),
    plot.margin = margin(20,20,20,20) # Margins around plot
  )

# Replace grouping labels
env.df$Years <- factor(env.df$Years,
                       levels = c("DG","RG1","RG2","RG3","CK"),
                       labels = c("Degraded","Early","Middle","Late","Intact"))

# Define colors (original colors + dark outline)
cols <- c("#be564a","#f9b483","#576289","#6c958f","#3e9fd1")
cols_darker <- darken(cols, amount = 0.4)  # Darken by 40%

# PCoA plot ---------------------------------------------------------
p1 <- ggplot(data = dt1, aes(x = PCoA1, y = PCoA2, fill = group)) +
  geom_point(shape = 21, size = 5, stroke = 0.8,
             color = cols_darker[as.numeric(factor(dt1$group))]) +  # Dark outline
  
  # Modified: hull border transparent, fill retained
  geom_polygon(
    data = dt1 %>% group_by(group) %>% slice(chull(PCoA1, PCoA2)),
    aes(x = PCoA1, y = PCoA2, group = group, fill = group),
    color = NA,         # Transparent border
    alpha = 0.2,        # Semi-transparent fill
    size = 1, 
    linetype = "solid"
  ) +
  
  geom_point(aes(x = PCoA1.mean, y = PCoA2.mean, fill = group),
             shape = 24, size = 6, stroke = 0.9,
             color = cols_darker[as.numeric(factor(dt1$group))]) + # Mean point outline
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
              color = cols_darker[as.numeric(distence$group)], # Dark outline
              width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = cols) +
  labs(y = "Distance to centroid", x = NULL) +
  mytheme +
  theme(legend.position = "none")

# Arrange with patchwork
combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1))

# Save to PDF
ggsave("PCoA_and_Distance_combined_function.pdf", combined, width = 16, height = 6)



