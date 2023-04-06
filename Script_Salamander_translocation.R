#################################################################################
##      Microbiota dynamics during translocation of A. maculatum larvae,       ##
##                      Fieschi-Meric et al., 2022                             ##
#################################################################################

# Code for preprocessing the sequences is currently available upon request 
# and will be published in Open-Source after our other articles using the same
# dataset are published.

load("Larval_microbiota.RData")
load("Larval_microbiota_unrarefied.RData")
library(ggplot2)
library(phyloseq)
library(dplyr)
library(ggpubr)
library(btools)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(gratia)
library(mgcv)
library(magrittr)
library(data.table)
library(ape)
library(rbiom)
library(dendextend)
library(vegan)
library(ade4)
library(forcats)
library(dunn.test)
library(pairwiseAdonis)
library(MiscMetabar)
library(ggVennDiagram)
library(MicrobiotaProcess)
library(microbiomeMarker)
library(speedyseq)
library(multcomp)
library(DESeq2)


### ------------------------------- Step 1: Datasets prep -----------------------

# Initial samples -> variation bt sites at D0
initial_ps = subset_samples(larvae_ps, (Lake_transfer == "INI"))
initial_ps = filter_taxa(initial_ps, function(x) sum(x) > 0, TRUE)
initial_inhib_ps = subset_taxa(initial_ps, Bd_inhibition == "inhibitory")

# Ctrl + Transferred samples -> variation at D15
transferred_ps = subset_samples(larvae_ps, (Date == "D15"))
transferred_ps = filter_taxa(transferred_ps, function(x) sum(x) > 0, TRUE)
transferred_inhib_ps = subset_taxa(transferred_ps, Bd_inhibition == "inhibitory")

# All samples
final_larvae_inhib_ps = subset_taxa(larvae_ps, Bd_inhibition == "inhibitory")


### ------------------------------- Step 2: Venn diagramms ----------------------

# Shared phylotypes between sites at D0
vennlist <- get_vennlist(initial_ps,
                         factorNames="LakeOrigin")
ggVennDiagram(vennlist, label_alpha = 0, label_size = 5, label = "count", edge_size = 2) 

# Shared Bd-inhibitory phylotypes between sites at D0
vennlist <- get_vennlist(initial_inhib_ps,
                         factorNames="LakeOrigin")
ggVennDiagram(vennlist, label_alpha = 0, label_size = 5, label = "count", edge_size = 2)


### ------------------------------- Step 3: Relative abundance plots ------------

# All phylotypes - at phylum level - FIGURE 2
ps_phylum = tax_glom(larvae_ps,taxrank = "Phylum")
ps_phylum = ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
ps_phylum = psmelt(ps_phylum)
ps_barplot = ps_phylum %>%
  group_by(Phylum, Origin, Final_lake) %>%
  summarise(nb = sum(Abundance))
ps_barplot$Phylum = as.character(ps_barplot$Phylum)
ps_barplot$Phylum[ps_barplot$nb < 0.005] = "Phylum < 0.5% abund."
ggplot(ps_barplot %>% mutate(Origin = factor(Origin, levels = c('INI', 'CTRL', 'BL', 'SP', 'LL'))),
               aes(x = Origin, y = nb, fill = Phylum, na.rm = TRUE)) +
  geom_bar(stat="identity",na.rm = TRUE, position="fill") +
  facet_grid(~ Final_lake, scales="free") +
  ylab("Relative Abundance \n") + xlab("Origin") 

# Bd-inhibitory phylotypes - phylum level - FIGURE 3
taxo <- data.frame(as(tax_table(larvae_ps), "matrix"))
temp_larvae_ps <- larvae_ps %>% transmute_tax_table(Phylum_activity = Phylum_activity, 
                                                          Bd_phylum = Bd_phylum)
ps_phylum = tax_glom(temp_larvae_ps,taxrank = "Bd_phylum")
ps_phylum = ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
ps_phylum = psmelt(ps_phylum)
ps_barplot = ps_phylum %>%
  group_by(Bd_phylum, Origin, Final_lake) %>%
  summarise(nb = sum(Abundance))
ps_barplot$Bd_phylum = factor(ps_barplot$Bd_phylum, levels = c('not inhibitory', 'unknown','Bacteroidetes', 'Proteobacteria'))
ps_barplot <- ps_barplot %>%
  group_by(Origin, Final_lake) %>%
  mutate(tot=sum(nb)) %>%
  ungroup() %>%
  mutate(freq = (nb/tot)*100) %>%
  filter(Bd_phylum == "Proteobacteria" | Bd_phylum == "Bacteroidetes" )
ggplot(ps_barplot %>% mutate(Origin = factor(Origin, levels = c('INI', 'CTRL', 'BL', 'SP', 'LL'))),
       aes(x = Origin, y = freq, fill = Bd_phylum)) +
  geom_bar(stat="identity") +
  facet_grid(~ Final_lake, scales="free") +
  ylab("Relative Abundance (%)") + xlab("Origin")


### ------------------------------- Step 4: Models for alpha diversity ----------
#Idem for Bd-inhibitory ASVs (replace with final_larvae_inhib_ps or use the .csv files)

tab = estimate_richness(larvae_ps, split = TRUE, measures = NULL)
tab = rownames_to_column(tab, var = "X.SampleID") 
tab$X.SampleID <- gsub("X","",as.character(tab$X.SampleID))
metadata <- read.csv("metadata.csv", sep=";")
metadata <- metadata %>% mutate(X.SampleID = as.character(X.SampleID))
alpha = tab %>% inner_join(metadata, by = "X.SampleID")

# Correct the class of the covariates
lapply(alpha, class)
alpha <- alpha %>% 
  mutate(Individual_ID = factor(Individual_ID),
         LakeOrigin = factor(LakeOrigin),
         Lake_transfer2 = factor(Lake_transfer2),
         Status = factor(Status),
         Final_lake = factor(Final_lake),
         Date = factor(Date),
         Transfer = factor(Transfer))

#____Initial differences - Chao1
mod1 = aov(Chao1 ~ LakeOrigin + Total_lenght, data = alpha %>% filter(Status == "INI")) 
shapiro.test(residuals(object = mod1))
kruskal.test(Chao1 ~ LakeOrigin, data = alpha %>% filter(Status == "INI"))
alpha_ini <- alpha %>% filter(Status == "INI")
dunn.test(alpha_ini$Chao1, alpha_ini$LakeOrigin, method="bonferroni")
alphaini <- alpha %>% filter(Status == "INI")
cor.test(alphaini$Total_lenght, alphaini$Chao1, method = "pearson")

#____Initial differences - Shannon
mod2 = aov(Shannon ~ LakeOrigin + Total_lenght, data = alpha %>% filter(Status == "INI"))
shapiro.test(residuals(object = mod2))
bartlett.test(Shannon ~ LakeOrigin, data = alpha %>% filter(Status == "INI"))
summary(mod2)

#____Ontogenic differences - Chao1
mod3 = lmer(Chao1 ~ Transfer + Date + (1|LakeOrigin) + (1|Individual_ID), data = alpha)
shapiro.test(residuals(object = mod3))   
mod3 = lmer(log(Chao1) ~ Transfer + Date + (1|LakeOrigin) + (1|Individual_ID), data = alpha)
shapiro.test(residuals(object = mod3))
anova(mod3)

#____Ontogenic differences - Shannon
mod4 = lmer(Shannon ~ Transfer + Date + (1|LakeOrigin) + (1|Individual_ID), data = alpha)
shapiro.test(residuals(object = mod4))  
anova(mod4)

#____Translocation status differences - Chao1
alphawide <- read.csv2("alpha-wide.csv", sep=";")  # Must use wide format
alphawide <- alphawide %>%
  mutate(LakeOrigin = factor(LakeOrigin),
         Final_lake = factor(Final_lake),
         Final_status = factor(Final_status))
mod5 = lmer(Chao1.old ~ Chao1.ini + Final_status + (1|LakeOrigin), data = alphawide) 
shapiro.test(residuals(object = mod5))                             # not normal... p<0.05
mod5 = lmer(log(Chao1.old) ~ Chao1.ini + Final_status + (1|LakeOrigin), data = alphawide) 
shapiro.test(residuals(object = mod5))                             # normal!
anova(mod5)

#____Translocation status differences in Shannon
mod6 = lmer(Shannon.old ~ Shannon.ini + Final_status + (1|LakeOrigin), data = alphawide) 
shapiro.test(residuals(object = mod6))                             # normal !
anova(mod6)

#____Plots (replace Chao1/Shannon)
ggplot((alpha %>% mutate(Lake_transfer2 = factor(Lake_transfer2, levels=c("INI", "BL", "LL", "SP")))),
                aes(x = Lake_transfer2, y = Chao1, fill = LakeOrigin)) + geom_boxplot()


### ----------------------------- Step 5: Beta diversity ------------------------
#Idem for Bd-inhibitory ASVs (replace with final_larvae_inhib_ps)

pick_new_outgroup = function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") 
  treeDT =
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  new.outgroup = treeDT[which.max(length)]$id
  return(new.outgroup) }
my.tree = phy_tree(larvae_ps)
out.group = pick_new_outgroup(my.tree)
out.group ## [1] "758d7c8c941a9519570882956fa2af98"
new.tree = ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
phy_tree(larvae_ps) = new.tree

#___Initial differences between sites
rbiom_weighted_A = rbiom::unifrac(otu_table(initial_ps),
                                  weighted=TRUE, tree=new.tree)
ini_var_meta <- data.frame(sample_data(initial_ps))
adonis2(rbiom_weighted_A ~ LakeOrigin,
        data = ini_var_meta, permutations = 9999)
pairwise.adonis2(rbiom_weighted_A ~ LakeOrigin, data = ini_var_meta, permutations = 9999)
disp_A = betadisper(rbiom_weighted_A, ini_var_meta$LakeOrigin)
permutest(disp_A, pairwise = TRUE, permutations = 9999)   

#___Ontogenic vs. Transloc differences 
rbiom_weighted_B = rbiom::unifrac(otu_table(larvae_ps),
                                  weighted=TRUE, tree=new.tree)
onto_var_meta <- data.frame(sample_data(larvae_ps))
adonis2(rbiom_weighted_B ~ Date + Transfer, Blocks(LakeOrigin),
        data = onto_var_meta, permutations = 9999)
disp_B = betadisper(rbiom_weighted_B, onto_var_meta$Date)
permutest(disp_B, pairwise = TRUE, permutations = 9999)   
disp_B = betadisper(rbiom_weighted_B, onto_var_meta$Transfer)
permutest(disp_B, pairwise = TRUE, permutations = 9999)   
ord_B = ordinate(larvae_ps, method = "PCoA", distance = rbiom_weighted_B)
plot_ordination(larvae_ps, ord_B, type='samples',
                color="Lake_transfer", justDF = F) +
  stat_ellipse() +
  facet_grid(cols = vars(LakeOrigin)) +
  labs(color = "Final location (D15)")

#___Translocation effect 
rbiom_weighted_C = rbiom::unifrac(otu_table(transferred_ps),
                                  weighted=TRUE, tree=new.tree)
final_var_meta <- data.frame(sample_data(transferred_ps))
adonis2(rbiom_weighted_C ~ Status, Blocks(LakeOrigin),
        data = final_var_meta, permutations = 9999)
disp_C = betadisper(rbiom_weighted_C, final_var_meta$Status)
permutest(disp_C, pairwise = TRUE, permutations = 9999)   
ord_C = ordinate(transferred_ps, method = "PCoA", distance = rbiom_weighted_C)
plot_ordination(transferred_ps, ord_C, type='samples',
                color="Lake_transfer2", shape="Status", justDF = F) +
  stat_ellipse() +
  facet_grid(cols = vars(LakeOrigin)) +
  labs(color = "Final location (D15)")



### ----------------------------- Step 6: Differential abundance ----------------
library(Rmisc)

#___Ontogenic effect: Idem for effect of translocation, for Bd-inhibitory ASVs, for dead vs. surviving larvae (replace with other ps objects or use .csv files)

# Subset groups of interest within unrarefied data
natural_unrar_ps = subset_samples(unrar_final_larvae_ps, (Status == "INI" | Status == "CTRL"))
natural_unrar_ps = filter_taxa(natural_unrar_ps, function(x) sum(x) > 0, TRUE) # ABSOLUTELY KEY!!!!

# DESeQ
deseq_env <- phyloseq_to_deseq2(natural_unrar_ps, ~ Status)
deseq_env2 = DESeq(deseq_env, test= "Wald", fitType="parametric",sfType="poscounts")
resultsNames(deseq_env2)
res_onto <- results(deseq_env2, name = "Status_INI_vs_CTRL", alpha = 0.01, pAdjustMethod="BH")
summary(res_onto)
res_onto = res_onto[order(res_onto$padj, na.last = NA), ]    
sigtab01 = res_onto[(res_onto$padj < 0.05), ]                 
summary(sigtab01)                                             
sigtab01 = cbind(as(sigtab01, "data.frame"), as(tax_table(natural_unrar_ps)[rownames(sigtab01), ], "matrix"))
sigtab01 = sigtab01[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class","Order", "Family", "Genus", "Species")]
possigtab01 = sigtab01[(sigtab01$log2FoldChange >= 0) & (sigtab01$padj <= 0.05),]
negsigtab01 = sigtab01[(sigtab01$log2FoldChange <= 0) & (sigtab01$padj <= 0.05),]
possigtab01[["change"]] <- (rep(c("Increase"),dim(possigtab01)[1]))
negsigtab01[["change"]] <- (rep(c("Decrease"),dim(negsigtab01)[1]))
contrasts_onto = rbind(possigtab01,negsigtab01)
contrasts_onto <- read.csv("Deseq2_contrasts_ontogeny.csv") # can replace with other .csv files here 
contrasts_onto <- contrasts_onto %>% dplyr::rename("ASV" = "X")
x = tapply(contrasts_onto$log2FoldChange, contrasts_onto$Phylum, function(x) max(x))
x = sort(x, TRUE)
contrasts_onto$Phylum = factor(as.character(contrasts_onto$Phylum), levels = names(x))
x = tapply(contrasts_onto$log2FoldChange, contrasts_onto$Genus, function(x) max(x))
x = sort(x, TRUE)
contrasts_onto$Genus = factor(as.character(contrasts_onto$Genus), levels = names(x))
ggplot(contrasts_onto, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -90))


### ----------------------------- Step 7: Potential causes of mortality ---------

df_temp <- read.csv2("Larval_survival.csv", sep=";")
alpha <- read.csv2(file = "alpha-long.csv", sep=";")
alpha <- alpha %>% filter(Date == "D0") %>%
  dplyr::select(SampleID, Chao1, Shannon, Total_lenght)
df_survival <- inner_join(df_temp, alpha, by="SampleID")
df_survival <- df_survival %>% mutate(Status = factor(Status, ordered = FALSE),
                                      LakeOrigin = factor(LakeOrigin, ordered = FALSE),
                                      FinalLake = factor(FinalLake, ordered = FALSE),
                                      Individual_ID = factor(Individual_ID, ordered = FALSE))

# Logistic regression
mod.1 <- glm(Survived ~ Total_lenght + Shannon + Status + LakeOrigin + FinalLake, data = df_survival, family = binomial) 
summary(mod.1)
exp(-2.67714)
exp(-2.67714 +1.96*1.08174)
exp(-2.67714 -1.96*1.08174)
1/(exp(-2.67714))


### ----------------------------- Step 8: Environmental controls -----------------

load("Water_microbiota.RData")

# Shared phylotypes between sites
vennlist <- get_vennlist(water_ps,
                         factorNames="LakeOrigin")
ggVennDiagram(vennlist, label_alpha = 0, label_size = 5, label = "count", edge_size = 2) 

# Relative abundance plot - at phylum level
ps_phylum = tax_glom(water_ps,taxrank = "Phylum")
ps_phylum = ps_phylum %>% transform_sample_counts(function(x) {x/sum(x)} )
ps_phylum = psmelt(ps_phylum)
ps_barplot = ps_phylum %>%
  group_by(Phylum, LakeOrigin) %>%
  dplyr::summarise(nb = sum(Abundance))
ps_barplot$Phylum = as.character(ps_barplot$Phylum)
ps_barplot$Phylum[ps_barplot$nb < 0.005] = "Phylum < 0.5% abund."
ggplot(ps_barplot %>% mutate(LakeOrigin = factor(LakeOrigin, levels = c('Bat', 'LostRay', 'SpeckledTrout'))),
                aes(x = LakeOrigin, y = nb, fill = Phylum, na.rm = TRUE)) +
  geom_bar(stat="identity",na.rm = TRUE, position="fill") +
  ylab("Relative Abundance \n") + xlab("Lake") 

# Alpha diversity
tab = estimate_richness(water_ps, split = TRUE, measures = NULL)
tab = rownames_to_column(tab, var = "X.SampleID") 
tab$X.SampleID <- gsub("X","",as.character(tab$X.SampleID))     
metadata_water <- read.csv("metadata_water.csv", sep=";")
metadata_water <- metadata_water %>% mutate(X.SampleID = as.character(X.SampleID))
alpha = tab %>% inner_join(metadata_water, by = "X.SampleID")
lapply(alpha, class)
alpha <- alpha %>% mutate(LakeOrigin = factor(LakeOrigin))
# Statistical test (replace with Shannon)
mod1 = aov(Chao1 ~ LakeOrigin, data = alpha) 
shapiro.test(residuals(object = mod1))                              
summary(mod1)
# Plots (replace with Shannon)
ggplot((alpha %>% mutate(LakeOrigin = factor(LakeOrigin, levels=c("Bat", "LostRay", "SpeckledTrout")))),
                  aes(x = LakeOrigin, y = Chao1, fill = LakeOrigin)) +
  geom_boxplot() +
  labs(y = 'Chao1 index',
       x = 'Lake')

# Beta diversity
my.tree = phy_tree(water_ps)
out.group = pick_new_outgroup(my.tree)
out.group ## [1] "3a7bcfd3dd1c23c2efae4b6b3fcd36a3"
new.tree = ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
phy_tree(water_ps) = new.tree
rbiom_weighted = rbiom::unifrac(otu_table(water_ps), weighted=TRUE, tree=new.tree)
# Plot
ord = ordinate(water_ps, method = "PCoA", distance = rbiom_weighted)
plot_ordination(water_ps, ord, type='samples',
                          color="LakeOrigin", justDF = F) +
  stat_ellipse() +
  labs(color = "Lake") +
  ggforce::geom_mark_ellipse(aes(color = LakeOrigin))
# Statistical test
water_var_meta <- data.frame(sample_data(water_ps))
adonis2(rbiom_weighted ~ LakeOrigin,
        data = water_var_meta, permutations = 9999)
pairwise.adonis2(rbiom_weighted ~ LakeOrigin, data = water_var_meta, permutations = 9999)
disp = betadisper(rbiom_weighted, water_var_meta$LakeOrigin)
permutest(disp, pairwise = TRUE, permutations = 9999)   

# Differential abundance analysis (same for other pairs of lakes)
BL_LL_unrar_ps = subset_samples(water_ps, (LakeOrigin != "SpeckledTrout"))
BL_LL_unrar_ps = filter_taxa(BL_LL_unrar_ps, function(x) sum(x) > 0, TRUE) 
deseq_env <- phyloseq_to_deseq2(BL_LL_unrar_ps, ~ LakeOrigin)
deseq_env2 = DESeq(deseq_env, test= "Wald", fitType="parametric",sfType="poscounts")
resultsNames(deseq_env2)
res <- results(deseq_env2, name = "LakeOrigin_LostRay_vs_Bat", alpha = 0.01, pAdjustMethod="BH")
summary(res)
res = res[order(res$padj, na.last = NA), ]     
sigtab01 = res[(res$padj < 0.05), ]                 
summary(sigtab01)                                            
sigtab01 = cbind(as(sigtab01, "data.frame"), as(tax_table(BL_LL_unrar_ps)[rownames(sigtab01), ], "matrix"))
sigtab01 = sigtab01[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class","Order", "Family", "Genus", "Species")]
possigtab01 = sigtab01[(sigtab01$log2FoldChange >= 0) & (sigtab01$padj <= 0.05),]
negsigtab01 = sigtab01[(sigtab01$log2FoldChange <= 0) & (sigtab01$padj <= 0.05),]
possigtab01[["change"]] <- (rep(c("Increase"),dim(possigtab01)[1]))
negsigtab01[["change"]] <- (rep(c("Decrease"),dim(negsigtab01)[1]))
contrasts_onto = rbind(possigtab01,negsigtab01)
contrasts_onto
