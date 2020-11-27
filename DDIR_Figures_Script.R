###############################################################################################
##############################    Research Script for    ######################################
#### In-depth clinical and biological exploration of DNA Damage Immune Response (DDIR) as  ####
###############    a biomarker for oxaliplatin use in colorectal cancer    ####################
###########################----------------------------------##################################
###########################      Authour: Sudhir Malla       ##################################
###########################    Email: smalla01@qub.ac.uk     ##################################
##########      Molecular Pathology Research Group (Dunne Lab | dunne-lab.com)       ##########
######################      Lab Email : dunnelabqub@gmail.com       ###########################
###############################################################################################
#                                                                                             #
# Figure 2 : Lines 77 to 387                                                                  #
# Figure 3 : Lines 389 to 706                                                                 #
# Figure 4 : Lines 708 to 1090                                                                #
# Figure 5 : Lines 1092 to 1213                                                               #
#                                                                                             #
#---------------------------------------------------------------------------------------------#
#                                                                                             #
# Supplementary Figure 3 : Lines 1222 to 1489                                                 #
# Supplementary Figure 4 : Lines 1491 to 1758                                                 #
# Supplementary Figure 5 : Lines 1760 to 2164                                                 #
# Supplementary Figure 6 : Lines 2166 to 2656                                                 #
# Supplementary Figure 7 : Lines 2658 to 2753                                                 #
# Supplementary Figure 8 : Lines 2755 to 2849                                                 #
#                                                                                             #
###############################################################################################

## Load libraries
library(ggplot2)
library(beeswarm)
library(CePa) # was used to read .gct file
library(tidyverse) # data manipulation
library(plyr)
library(cowplot)
library(grid)
library(ggpubr)

## Put the Reproducible_script directory here - indicate where you want to save
setwd()

###############################################################################################

## Save default par before manipulating plots
default.par <- par()
options(scipen = 999) ## no scientific decimal system

## Cohorts used ---
cohort.1 <- "FOCUS"
cohort.2 <- "FOxTROT"
cohort.3 <- "TRANSBIG"

## data folder - Create a data folder and move data files here required for this script
data.folder <- file.path(getwd(), "data")
#if(!dir.exists(data.folder)) dir.create(data.folder, showWarnings = F, recursive = T)

focus.file.path <- file.path(data.folder,
                             "FOCUS_DDIR_n361_maxMean_MCP_ESTIMATE_ssGSEA_clino_expression_data.txt")
foxtrot.file.path <- file.path(data.folder,
                               "FOxTROT_DDIR_n97_maxMean_MCP_ESTIMATE_ssGSEA_clino_expression_data.txt")
transbig.file.path <- file.path(data.folder,
                                "TRANSBIG_DDIR_n198_maxMean_MCP_ESTIMATE_ssGSEA_clino_expression_data.txt")


## Import ---
focus.data <- read.delim(focus.file.path, header = T, stringsAsFactors = F)
foxtrot.data <- read.delim(foxtrot.file.path, header = T, stringsAsFactors = F)
transbig.data <- read.delim(transbig.file.path, header = T, stringsAsFactors = F)


## Change the data.1 with the data from three cohort as needed
#data.1 <- focus.data
#data.1 <- foxtrot.data
#data.1 <- transbig.data


###############################################################################################
######################################   Figure 2    ##########################################
###############################################################################################

######################################   Figure 2A   ##########################################

###  Figure 2A) CMS vs DDIR scores - boxplot - FOCUS --->>>

data.1 <- focus.data
data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

## Statistics:

## Tukey's HSD test
## form a model (one-way ANOVA)
model.1 <- lm(DDIR_scores ~ CMS, data = data.1)
summary(model.1)
ANOVA.1 <- aov(model.1, data = data.1)
summary(ANOVA.1)
## posthoc tukey
posthoc.1 <- TukeyHSD(ANOVA.1, 'CMS', conf.level = 0.95)
Tukey.result.1 <- posthoc.1$CMS
## Setting p-value
options(scipen = 999)
pvalue.1 <- ifelse(Tukey.result.1[1,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[1,4], 4)),
                   paste0("P < 0.0001"))
pvalue.2 <- ifelse(Tukey.result.1[2,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[2,4], 4)),
                   paste0("P < 0.0001"))
pvalue.3 <- ifelse(Tukey.result.1[3,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[3,4], 4)),
                   paste0("P < 0.0001"))
pvalue.4 <- ifelse(Tukey.result.1[4,4] > 0.0001,
                   paste0("P = ", round(Tukey.result.1[4,4], 4)),
                   paste0("P < 0.0001"))

kruskal.1 <- kruskal.test(data.1[["DDIR_scores"]] ~ as.factor(data.1[["CMS"]]),
                          data = data.1)
pvalue.5 <- ifelse(kruskal.1$p.value > 0.0001, 
                   paste0("Kruskal-Wallis, p = ", round(kruskal.1$p.value, 4)),
                   paste0("Kruskal-Wallis, p < 0.0001"))

### ---

#tiff("Figure_2A.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(2.1,4.1,4,2.1))

boxplot(data.1[["DDIR_scores"]] ~ data.1[["CMS"]], data = data.1,
        ylim = c(-0.30, 1), cex.axis = 1.2, outline = FALSE,
        ylab = "", xlab = "")
beeswarm(data.1[["DDIR_scores"]] ~ data.1[["CMS"]], data = data.1,
         add = TRUE, pch = 19, cex = 1.6, spacing = 1.2,
         method = "swarm", corral="wrap",
         col = c("#E18C28", "#0B5F9C", "#C46294", "#129063", "#A4A4A4"))

## Add titles
title(main = paste0("DDIR score by CMS subtypes\n",cohort.1," trial cohort"),
      line = 1, cex.main = 1.8)
title(ylab = paste0("DDIR Scores"), line = 2.6, cex.lab = 1.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add stats
segments(1,0.6,2,0.6, lwd = 2)
text(1.5, 0.63, pvalue.1, cex = 1.2)
segments(1,0.7,3,0.7, lwd = 2)
text(2, 0.73, pvalue.2, cex = 1.2)
segments(1,0.8,4,0.8, lwd = 2)
text(2.5, 0.83, pvalue.3, cex = 1.2)
segments(1,0.9,5,0.9, lwd = 2)
text(3, 0.93, pvalue.4, cex = 1.2)
mtext(pvalue.5, side = 3, line = -1.4, adj = 0.98, cex = 1.2)

## dev off
dev.off()


###############################################################################################
######################################   Figure 2B   ##########################################

###  Figure 2B) CMS vs DDIR scores - boxplot - FOxTROT --->>>

data.1 <- foxtrot.data
data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

## Statistics:

## Tukey's HSD test
## form a model (one-way ANOVA)
model.1 <- lm(DDIR_scores ~ CMS, data = data.1)
summary(model.1)
ANOVA.1 <- aov(model.1, data = data.1)
summary(ANOVA.1)
# posthoc tukey
posthoc.1 <- TukeyHSD(ANOVA.1, 'CMS', conf.level = 0.95)
Tukey.result.1 <- posthoc.1$CMS
## Setting p-value
options(scipen = 999)
pvalue.1 <- ifelse(Tukey.result.1[1,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[1,4], 4)),
                   paste0("P < 0.0001"))
pvalue.2 <- ifelse(Tukey.result.1[2,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[2,4], 4)),
                   paste0("P < 0.0001"))
pvalue.3 <- ifelse(Tukey.result.1[3,4] > 0.0001, 
                   paste0("P = ", round(Tukey.result.1[3,4], 4)),
                   paste0("P < 0.0001"))
pvalue.4 <- ifelse(Tukey.result.1[4,4] > 0.0001,
                   paste0("P = ", round(Tukey.result.1[4,4], 4)),
                   paste0("P < 0.0001"))

kruskal.1 <- kruskal.test(data.1[["DDIR_scores"]] ~ as.factor(data.1[["CMS"]]),
                          data = data.1)
pvalue.5 <- ifelse(kruskal.1$p.value > 0.0001, 
                   paste0("Kruskal-Wallis, p = ", round(kruskal.1$p.value, 4)),
                   paste0("Kruskal-Wallis, p < 0.0001"))

### ---

#tiff("Figure_2B.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(2.1,4.1,4,2.1))

boxplot(data.1[["DDIR_scores"]] ~ data.1[["CMS"]], data = data.1,
        ylim = c(-0.30, 1), cex.axis = 1.2, outline = FALSE,
        ylab = "", xlab = "")
beeswarm(data.1[["DDIR_scores"]] ~ data.1[["CMS"]], data = data.1,
         add = TRUE, pch = 19, cex = 1.6, spacing = 1.2,
         method = "swarm", corral="wrap",
         col = c("#E18C28", "#0B5F9C", "#C46294", "#129063", "#A4A4A4"))

## Add titles
title(main = paste0("DDIR score by CMS subtypes\n",cohort.2," trial cohort"),
      line = 1, cex.main = 1.8)
title(ylab = paste0("DDIR Scores"), line = 2.6, cex.lab = 1.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add stats
segments(1,0.66,2,0.66, lwd = 2)
text(1.5, 0.69, pvalue.1, cex = 1.2)
segments(1,0.74,3,0.74, lwd = 2)
text(2, 0.77, pvalue.2, cex = 1.2)
segments(1,0.82,4,0.82, lwd = 2)
text(2.5, 0.85, pvalue.3, cex = 1.2)
segments(1,0.90,5,0.90, lwd = 2)
text(3, 0.93, pvalue.4, cex = 1.2)
mtext(pvalue.5, side = 3, line = -1.4, adj = 0.98, cex = 1.2)

## dev off
dev.off() 


###############################################################################################
##################################   Figure 2C and 2E   #######################################

### Figure 2C)  MSI status vs DDIR subtypes Bargraph - FOCUS --->>>
### Figure 2E)  MSI status vs DDIR subtypes Bargraph - FOxTROT --->>>

## Select ONE: 
data.1 <- focus.data # For Fig 2C
#data.1 <- foxtrot.data # For Fig 2E

## Form table for MSI status in DDIR POS vs DDIR NEG
table.1 <- table(data.1[["MSI"]], data.1[["DDIR"]]) # for fishers test, removing NA
table.2 <- table(data.1[["MSI"]], data.1[["DDIR"]], useNA="ifany") #for proportion and graph
## Proportion table for graph
prop.1 <- prop.table(table.2, 2)
rownames(prop.1)[3] <- "NA" # rename the third row name
## Turn proportion table into long format
prop.long.1 <- as.data.frame(prop.1)


## Statistics:
## Fisher's Exact test (MSI) ## everthing left at default
fisher.1 <- fisher.test(table.1)
pvalue.1 <- ifelse(fisher.1$p.value > 0.0001, 
                   paste0("Fisher's Exact test, p = ",round(fisher.1$p.value, 4)),
                   "Fisher's Exact, p < 0.0001")

### ---

#tiff("Figure_2C.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_2E.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
ggplot(prop.long.1) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1),data = prop.long.1, 
           stat = "identity", width = 0.7,
           position = position_stack(reverse = TRUE))+
  ## Add colour and legend
  scale_fill_manual(values = c("darkorange", "blue4", "darkgrey"),
                    guide=guide_legend(reverse=TRUE, title = paste0('MSI'))) +
  ## Add limit to y-axis
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  ## Add titles
  xlab(" ") +
  ylab("Sample Proportions") +
  labs(title = paste0("MSI status in DDIR\n",
                      ## Select cohort.1 for FOCUS, cohort.2 for FOxTROT
                      cohort.1,
                      #cohort.2,
                      " trial cohort")) +
  
  
  ## Theme for main title
  theme(plot.title = element_text(face = 'bold', size = 22,
                                  margin = margin(15,0,0,0), 
                                  vjust = 4.3, hjust = 0.5, lineheight = 1)) +
  
  ## Add p-value
  labs(caption = pvalue.1)+
  theme(plot.caption = element_text(size = 20, 
                                    margin = margin(15,0,0,0), hjust = 0.5)) +
  
  ## Plot background theme
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, vjust = 3),
        axis.text.y = element_text(size = 18, colour = "black", angle = 90, hjust = 0.5),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white'),
        
        ##legend
        legend.background = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = 'bold'),
        legend.title.align = 0,
        legend.position = c(0.95,0.78),
        legend.key.height = unit(1.0, 'cm'),
        legend.key.width = unit(1.0, 'cm'),
        legend.margin = margin(-25,0,0,35),
        plot.margin = unit(c(5,50,5,5),'mm'))


## dev off
dev.off()


###############################################################################################
##################################   Figure 2D and 2F   #######################################

### Figure 2D)  MSI status vs DDIR score boxplot - FOCUS --->>>
### Figure 2F)  MSI status vs DDIR score boxplot - FOxTROT --->>>

## Select ONE: 
data.1 <- focus.data # For Fig 2D
#data.1 <- foxtrot.data # For Fig 2F

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72

## Statistics:
t.test.1 <- t.test(DDIR_scores ~ MSI, data = data.1)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
wilcox.1 <- wilcox.test(DDIR_scores ~ MSI, data = data.1)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")

### ---

#tiff("Figure_2D.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_2F.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(5.1,4.5,4.1,2.1)+0.5)

boxplot(DDIR_scores ~ MSI, data = data.1, 
        outline = FALSE, ylim = c(y.min.1, y.max.1),
        cex.axis = 1.6,cex.lab = 1.6)
beeswarm(DDIR_scores ~ MSI, data = data.1, 
         method = "swarm", corral = "wrap", pch = 19,
         col = c("darkorange", "blue4"),cex = 1.5, spacing = 1.2,
         add = TRUE)

## Add titles
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("MSI status vs DDIR scores\n",
                    ## Select cohort.1 for FOCUS, cohort.2 for FOxTROT
                    cohort.1,
                    #cohort.2,
                    " trial cohort"), cex.main = 1.8)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add stats
mtext(pvalue.1, side = 1, line = 2.8, cex = 1.6)
mtext(pvalue.2, side = 1, line = 4.2, cex = 1.6)

## dev off
dev.off()


###############################################################################################
######################################   Figure 3    ##########################################
###############################################################################################

## Import GSEA findings ---
## GSEA file paths
focus.gsea.filepath.1 <- file.path(data.folder,
                                   "FOCUS_gsea_report_for_DDRD_POS_1555345550304.xls")
focus.gsea.filepath.2 <- file.path(data.folder,
                                   "FOCUS_gsea_report_for_DDRD_NEG_1555345550304.xls")
foxtrot.gsea.filepath.1 <- file.path(data.folder,
                                     "FOxTROT_gsea_report_for_DDRD_POS_1560787000819.xls")
foxtrot.gsea.filepath.2 <- file.path(data.folder,
                                     "FOxTROT_gsea_report_for_DDRD_NEG_1560787000819.xls")
transbig.gsea.filepath.1 <- file.path(data.folder,
                                      "TRANSBIG_gsea_report_for_DDRD_POS_1579086912188.xls")
transbig.gsea.filepath.2 <- file.path(data.folder,
                                      "TRANSBIG_gsea_report_for_DDRD_NEG_1579086912188.xls")

## import
focus.gsea.pos <- read.delim(focus.gsea.filepath.1, header = T, stringsAsFactors = F)
focus.gsea.neg <- read.delim(focus.gsea.filepath.2, header = T, stringsAsFactors = F)
foxtrot.gsea.pos <- read.delim(foxtrot.gsea.filepath.1, header = T, stringsAsFactors = F)
foxtrot.gsea.neg <- read.delim(foxtrot.gsea.filepath.2, header = T, stringsAsFactors = F)
transbig.gsea.pos <- read.delim(transbig.gsea.filepath.1, header = T, stringsAsFactors = F)
transbig.gsea.neg <- read.delim(transbig.gsea.filepath.2, header = T, stringsAsFactors = F)

###############################################################################################
###############################################################################################

## Plot with combined GSEA results for three cohorts with geneset 
## significant in at least one cohort.

## FOCUS
focus.posneg <- rbind(focus.gsea.pos, focus.gsea.neg)
head(focus.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
focus.posneg$NAME <- focus.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(focus.posneg) # to select column with NAME, NES, FDR, size
FOCUS <- subset(focus.posneg, select = c(1,4,6,8))
head(FOCUS)
FOCUS <- FOCUS %>% mutate(Cohort = "FOCUS") # add column named 'Cohort'
FOCUS.fdrsig <- FOCUS[FOCUS$FDR.q.val < 0.25 ,] # only FDR sig genesets

## FOxTROT
foxtrot.posneg <- rbind(foxtrot.gsea.pos, foxtrot.gsea.neg)
head(foxtrot.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
foxtrot.posneg$NAME <- foxtrot.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(foxtrot.posneg)
FOxTROT <- subset(foxtrot.posneg, select = c(1,4,6,8))
head(FOxTROT)
FOxTROT <- FOxTROT %>% mutate(Cohort = "FOxTROT")
FOxTROT.fdrsig <- FOxTROT[FOxTROT$FDR.q.val < 0.25 ,]

## TRANSBIG
transbig.posneg <- rbind(transbig.gsea.pos, transbig.gsea.neg)
head(transbig.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
transbig.posneg$NAME <- transbig.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(transbig.posneg)
TRANSBIG <- subset(transbig.posneg, select = c(1,4,6,8))
head(TRANSBIG)
TRANSBIG <- TRANSBIG %>% mutate(Cohort = "TRANSBIG")
TRANSBIG.fdrsig <- TRANSBIG[TRANSBIG$FDR.q.val < 0.25 ,]

## Find out common and different pathway between cohorts
## common and different (FOCUS, FOxTROT, and TRANSBIG)
FcFxTr <- intersect(intersect(FOCUS.fdrsig$NAME, FOxTROT.fdrsig$NAME), 
                    TRANSBIG.fdrsig$NAME) # common in all three

FcFx <- intersect(FOCUS.fdrsig$NAME, FOxTROT.fdrsig$NAME) # common in Fc and Fx
## NOTE: FcFx do not have anything in common than FcFxTr

FcTr <- intersect(FOCUS.fdrsig$NAME, TRANSBIG.fdrsig$NAME) # common in Fc and Tr
FcTr <- FcTr[!FcTr %in% FcFxTr] # common in Fc and Tr excluding FcFxTr

FxTr <- intersect(FOxTROT.fdrsig$NAME, TRANSBIG.fdrsig$NAME) # common in Fx and Tr
FxTr <- FxTr[!FxTr %in% FcFxTr] # common in Fx and Tr excluding FcFxTr

## FOCUS only
Fconly <- FOCUS.fdrsig$NAME
Fconly <- Fconly[!Fconly %in% FcFxTr]
Fconly <- Fconly[!Fconly %in% FcTr]

## FOxTROT only
FxOnly <- FOxTROT.fdrsig$NAME
FxOnly <- FxOnly[!FxOnly %in% FcFxTr]
FxOnly <- FxOnly[!FxOnly %in% FxTr]

## TRANSBIG only
TrOnly <- TRANSBIG.fdrsig$NAME
TrOnly <- TrOnly[!TrOnly %in% FcFxTr]
TrOnly <- TrOnly[!TrOnly %in% FcTr]
TrOnly <- TrOnly[!TrOnly %in% FxTr]


## rbind three cohort
all3cohorts <- rbind(FOCUS, FOxTROT, TRANSBIG)
dim(all3cohorts)
head(all3cohorts)


## add new variable column
all3cohort2 <- all3cohorts %>% mutate(Common = as.character(ifelse(all3cohorts$NAME %in% FcFxTr, "FcFxTr",
                                                                   ifelse(all3cohorts$NAME %in% FcTr, "FcTr",
                                                                          ifelse(all3cohorts$NAME %in% FxTr, "FxTr",
                                                                                 ifelse(all3cohorts$NAME %in% Fconly, "Fconly",
                                                                                        ifelse(all3cohorts$NAME %in% FxOnly, "FxOnly",
                                                                                               ifelse(all3cohorts$NAME %in% TrOnly, "TrOnly", NA))))))))


## select fdrsig hallmarks from all cohorts
FDRsig.name <- as.character(all3cohort2[as.numeric(as.character(all3cohort2$FDR.q.val)) < 0.25 ,]$NAME)
FDRsig.name <- FDRsig.name[!duplicated(FDRsig.name)] #remove duplicates

## select these fdr.sig hallmark from all cohorts
all3fdr <- all3cohort2[all3cohort2$NAME %in% FDRsig.name ,]

## remove strings
all3fdr$NAME <- gsub("HALLMARK ", "", gsub("_", " ", all3fdr$NAME))
head(all3fdr)

### ---

### Combined dot plot ---

######################################   Figure 3A   ##########################################

### Figure 3A) Dot plot with common gene sets in all three cohort and only in BC cohort --->>>

par(mar=c(5.1,3.1,4.1,2.1))

## All common
p1a <- all3fdr[all3fdr$Common == "FcFxTr" ,] %>%
  ggplot(aes(x = Cohort,y = reorder(NAME, NES), color = NES, size = 5,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) + # leave blank if ns
  scale_shape_identity()+
  geom_point(stat = "identity")+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 12, angle = 0, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 0, face = "bold",
                                   margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position="right",
        legend.key = element_rect(fill = "white")) +
  theme(legend.position="none",  plot.margin = unit(c(0,0.5,0,6), "lines"))+
  ggtitle("GSEA - Hallmarks") +
  theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, vjust = 0,
                                  margin = ggplot2::margin(0,0,0.5,0,'cm'))) +
  scale_x_discrete(position = "top") 


## adding side line and text
p1b <- p1a +
  annotation_custom(grob = linesGrob(gp = gpar(col='red', lwd = 3)), 
                    ymin = 0.8, ymax = 5.2, xmin = -2, xmax = -2)
txt.g1a <- text_grob(label = "ALL", color = 'red', face = 'bold',
                     rot = 90, size = 20, hjust = 0.5)
p1c <- p1b + 
  annotation_custom(txt.g1a, ymin = 0.8, ymax = 5.2, xmin = -2.2, xmax = -2.2)

## Code to override clipping
gt.p1b <- ggplotGrob(p1c)
gt.p1b$layout[grepl("panel", gt.p1b$layout$name), ]$clip <- "off"

## Draw the plot
grid.newpage()
grid.draw(gt.p1b)

### ---

## TRANSBIG only

p4a <- all3fdr[all3fdr$Common == "TrOnly" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES, size = 5,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity")+
  guides(size = FALSE, color = guide_colourbar(title = 'NES'), shape = FALSE)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 12, angle = 0, face = "bold")) +
  theme(legend.position="none",
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = "white"),
        plot.margin = unit(c(0,0.5,0,6), "lines"))

## adding side line and text
p4b <- p4a +
  annotation_custom(grob = linesGrob(gp = gpar(col='pink2', lwd = 3)),
                    ymin = 0.8, ymax = 8.2, xmin = -2, xmax = -2)
txt.g2a <- text_grob(label = "TRANSBIG", color = 'pink2', face = 'bold',
                     rot = 90, size = 20, hjust = 0.1)
p4c <- p4b + 
  annotation_custom(txt.g2a, ymin = 0.8, ymax = 6.2, xmin = -2.2, xmax = -2.2)

## Code to override clipping
gt.p4b <- ggplotGrob(p4c)
gt.p4b$layout[grepl("panel", gt.p4b$layout$name), ]$clip <- "off"

## Draw the plot
grid.newpage()
grid.draw(gt.p4b)

### ---
## Combine the separate plots into one
## plot separate legend 
legend_p4a <- get_legend(p4a +   theme(legend.position="right",
                                       legend.text = element_text(size = 10),
                                       legend.key = element_rect(fill = "white"),
                                       plot.margin = unit(c(0,0.5,0,6), "lines")))

plot1a <- plot_grid(gt.p1b,gt.p4b, align = "v", nrow = 2, 
                    rel_heights = c(0.4, 0.5))
plot2a <- plot_grid(plot1a, legend_p4a, ncol = 2, 
                    rel_widths = c(1, 0.1))

### ---

#tiff("Figure_3A.tiff", height = 16, width = 23, units = 'cm', res = 800)
plot2a
dev.off()


###############################################################################################
#################################   Figure 3B, 3C and 3D  #####################################


### Figure 3B)  CXCL10 expression vs DDIR score - TRANSBIG --->>>
### Figure 3C)  CXCL10 expression vs DDIR score - FOCUS --->>>
### Figure 3D)  CXCL10 expression vs DDIR score - FOxTROT --->>>

data.1 <- transbig.data
#data.1 <- focus.data
#data.1 <- foxtrot.data

variable.1 <- "CXCL10"

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
## min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72


## Statistics:
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001") 

### ---

#tiff("Figure_3B.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_3C.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_3D.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(5.1,4.1,4.1,2.1))

plot(data.1[["DDIR_scores"]] ~ data.1[[variable.1]], data = data.1,
     ylab = "", xlab = "", 
     ## Select ylim for FOCUS and FOxTROT only
     #ylim = c(y.min.1, y.max.1),
     cex = 1.6, cex.axis = 1.6, pch = 19)

## Add titles
title(xlab = paste0(variable.1, ' expression (maxMean)'), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0(variable.1, " expression vs DDIR scores\n",
                    ## Select as needed for each cohort
                    cohort.3,
                    #cohort.1,
                    #cohort.2,
                    " trial cohort"), line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)), 
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.1), side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (BC = 0.37; CRC = 0.1094)
abline(h=0.37, col = "red", lty = 2, lwd = 3) # TRANSBIG
#abline(h=0.1094, col = "red", lty = 2, lwd = 3) # FOCUS and FOxTROT

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]), lwd = 2, col = "blue")

## dev off
dev.off()


###############################################################################################
######################################   Figure 4    ##########################################
###############################################################################################


###############################################################################################
#################################   Figure 4A, 4B and 4C  #####################################

### Figure 4A)  MCP scores vs DDIR scores - TRANSBIG --->>>
### Figure 4B)  MCP scores vs DDIR scores - FOCUS --->>>
### Figure 4C)  MCP scores vs DDIR scores - FOxTROT --->>>

data.1 <- transbig.data
#data.1 <- focus.data
#data.1 <- foxtrot.data
variable.1 <- "T.cells"
variable.2 <- "B.lineage"
variable.3 <- "Monocytic.lineage"

## Statistics:
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")
pearson.2 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.2]], data = data.1,
                      alternative = "two.sided", method = "pearson")
pearson.3 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.3]], data = data.1,
                      alternative = "two.sided", method = "pearson")


## pvalue
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("p = ",round(pearson.1$p.value, 4)),
                   "p < 0.0001") 
pvalue.2 <- ifelse(pearson.2$p.value > 0.0001,
                   paste0("p = ",round(pearson.2$p.value, 4)),
                   "p < 0.0001") 
pvalue.3 <- ifelse(pearson.3$p.value > 0.0001,
                   paste0("p = ",round(pearson.3$p.value, 4)),
                   "p < 0.0001") 

### ---

## Determine MCP (x-axis) range
empty.list.1 <- list()
for(i in c(variable.1, variable.2, variable.3)){
  
  a <- data.frame(MCP = paste0(i), min = min(data.1[[i]]), max = max(data.1[[i]]))
  empty.list.1[[i]] <- a
  
}
x.min.max.1 <- do.call(rbind, empty.list.1)
x.min.1 <- min(x.min.max.1[["min"]])-0.1
x.max.1 <- max(x.min.max.1[["max"]])+0.1

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72


#tiff("Figure_4A.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_4B.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_4C.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(7.1,4.1,3.5,2.1))

plot(data.1[["DDIR_scores"]] ~ data.1[[variable.1]],
     xlim=c(x.min.1,x.max.1),
     ## Select ylim for FOCUS and FOxTROT only
     #ylim = c(y.min.1, y.max.1),
     col="red", pch = 19,ylab = "", xlab = "", cex.axis = 1.6, cex = 1.6)
points(data.1[["DDIR_scores"]]  ~ data.1[[variable.2]],col="orange", pch = 19, cex = 1.6)
points(data.1[["DDIR_scores"]]  ~ data.1[[variable.3]],col="blue", pch = 19, cex = 1.6)

## Add titles
title(xlab = paste0("MCP Scores"), line = 2.4, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("MCP scores vs DDIR scores\n",
                    ## Select as needed for each cohort
                    cohort.3,
                    #cohort.1,
                    #cohort.2,
                    " trial cohort"),
      line = 0.6, cex.main = 1.8)

## Add legend for four immune cells
legend("bottom", legend = c(paste0('T cells, Pearson r = ',round(pearson.1$estimate, 4),
                                   ' , ', pvalue.1),
                            paste0('B Lineage, Pearson r = ',round(pearson.2$estimate, 4),
                                   ' , ', pvalue.2),
                            paste0('M Lineage, Pearson r = ',round(pearson.3$estimate, 4),
                                   ' , ', pvalue.3)),
       col = c("red", "orange", "blue"), pch = 19, bty = "n", 
       xpd = TRUE, inset = -0.363, ncol=1, cex = 1.2, x.intersp = 1, y.intersp = 1, 
       xjust = 0.5, yjust = 0)

## Add DDIR threshold (BC = 0.37; CRC = 0.1094)
abline(h=0.37, col = "red", lty = 2, lwd = 3) # TRANSBIG
#abline(h=0.1094, col = "red", lty = 2, lwd = 3) # FOCUS and FOxTROT

## as abline length cannot be controlled, clip the sides of the plot to cut lines off
clip(4.4, 11.8, -0.03,0.9) # TRANSBIG
#clip(3, 6, -0.35,0.64) # FOCUS
#clip(2.9, 7.4, -0.35,0.72) # FOxTROT

abline(lm(data.1[["DDIR_scores"]]  ~ data.1[[variable.1]]), col = "red", lwd = 2)
abline(lm(data.1[["DDIR_scores"]]  ~ data.1[[variable.2]]), col = "orange", lwd = 2)
abline(lm(data.1[["DDIR_scores"]]  ~ data.1[[variable.3]]), col = "blue", lwd = 2)

## dev off
dev.off()


###############################################################################################
######################################   Figure 4D   ##########################################

###  Figure 4D) Cytotoxic Lymphocytes MCP score vs DDIR - barplot/scatterplot - TRANSBIG --->>>

data.1 <- transbig.data
variable.1 <-  "Cytotoxic.lymphocytes"

## Statistics:
t.test.1 <-  t.test(data.1[[variable.1]] ~ DDIR, data = data.1)
wilcox.1 <- wilcox.test(data.1[[variable.1]] ~ DDIR, data = data.1)
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")
pvalue.3 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001")       

### ---

#tiff("Figure_4D.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(3.1,3.1,3.1,4.5)+1.5)

## Plot 1 (Scatterplot):
par(fig=c(0,1,0,0.85))
plot(data.1[[variable.1]], data.1[["DDIR_scores"]],
     pch = 19, ylab = "", xlab = "", cex = 1.6, cex.axis = 1.6, 
     ylim = c(min(data.1[["DDIR_scores"]])-0.1, max(data.1[["DDIR_scores"]])+0.1))

## Add titles
title(ylab = paste0("DDIR Scores"), line =2.8, cex.lab = 1.6)
title(xlab = paste0("Cytotoxic Lymphocytes MCP scores"), line = 2.8, cex.lab = 1.6)

## Add stats 1 (pearson)
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.3), side = 3,adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (BC = 0.37)
abline(h=0.37, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]), lwd = 2)

## Plot 2 (boxplot):
par(fig=c(0,1,0.65,1), new=TRUE) # adding another graph to the same plot
boxplot(data.1[[variable.1]] ~ data.1[["DDIR"]],  horizontal=TRUE, axes=FALSE, 
        col = c("blue", "red"))

## Add legend 2 (DDIR)
legend("topright", legend = c("POS", "NEG"), title = list("DDIR"),
       col = c("red", "blue"),pch = 15, title.adj = 0.78, y.intersp = 0.9,
       bty = "n", cex = 1.4, xpd = TRUE,inset = c(-0.22, -1.2), adj = c(0.1,0.5))

## Add titles
title(main =  paste0("Cytotoxic Lymphocytes MCP scores\n",cohort.3," trial cohort"),
      line = 1.4, cex.main = 1.8)

## Add stats 2 (t.test)
mtext(pvalue.1, side = 1, line = 0.5, cex = 1.4, adj = 0.02)
mtext(pvalue.2, side = 1, line = 1.7, cex = 1.4, adj = 0.02)

## dev off
dev.off()


###############################################################################################
###################################   Figure 4E and 4F   ######################################

###  Figure 4E)  Cytotoxic Lymphocytes MCP score vs DDIR - barplot/scatterplot - FOCUS --->>>
###  Figure 4F)  Cytotoxic Lymphocytes MCP score vs DDIR - barplot/scatterplot - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

variable.1 <-  "Cytotoxic.lymphocytes"

data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

## Statistics:
t.test.1 <-  t.test(data.1[[variable.1]] ~ DDIR, data = data.1)
wilcox.1 <- wilcox.test(data.1[[variable.1]] ~ DDIR, data = data.1)
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")
pvalue.3 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001")       


## Set up CMS colour
col.1 <- c("#E18C28","#0B5F9C","#C46294","#129063","#A4A4A4")
names(col.1) <- c("CMS1", "CMS2", "CMS3", "CMS4", "UNK")
col.1 <- col.1[data.1[, 'CMS']]

### ---

#tiff("Figure_4E.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Figure_4F.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(3.1,3.1,3.1,4.5)+1.5)

## Plot 1 (Scatterplot):
par(fig=c(0,1,0,0.85))
plot(data.1[[variable.1]], data.1[["DDIR_scores"]],
     pch = 19, col = col.1, ylab = "", xlab = "", cex = 1.6, cex.axis = 1.45, 
     ylim = c(-0.3, 0.65))

## Add legend 1 (CMS)
legend('topright', legend = c("CMS1", "CMS2", "CMS3", "CMS4", "UNK"), title = list("CMS"),
       col = c("#E18C28","#0B5F9C","#C46294","#129063","#A4A4A4"),pch = 19, title.adj = 0.55,
       bty = "n", cex = 1.4, xpd = TRUE, inset = c(-0.25, -0.03), adj = c(0.1, 0.5), y.intersp = 0.9)

## Add titles
title(ylab = paste0("DDIR Scores"), line =2.8, cex.lab = 1.6)
title(xlab = paste0("Cytotoxic Lymphocytes MCP scores"), line = 2.8, cex.lab = 1.6)

## Add stats 1 (pearson)
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.3), side = 3,adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]), lwd = 2)

## Plot 2 (boxplot):
par(fig=c(0,1,0.65,1), new=TRUE) # adding another graph to the same plot
boxplot(data.1[[variable.1]] ~ data.1[["DDIR"]],  horizontal=TRUE, axes=FALSE, 
        col = c("blue", "red"))

## Add legend 2 (DDIR)
legend("topright", legend = c("POS", "NEG"), title = list("DDIR"),
       col = c("red", "blue"),pch = 15, title.adj = 0.78, y.intersp = 0.9,
       bty = "n", cex = 1.4, xpd = TRUE,inset = c(-0.22, -1.2), adj = c(0.1,0.5))

## Add titles
title(main =  paste0("Cytotoxic Lymphocytes MCP scores\n",
                     ## select cohort.1 for FOCUS, and cohort.2 for FOxTROT
                     cohort.1,
                     #cohort.2,
                     " trial cohort"),
      line = 1.4, cex.main = 1.8)

## Add stats 2 (t.test)
mtext(pvalue.1, side = 1, line = 0.5, cex = 1.4, adj = 0.02)
mtext(pvalue.2, side = 1, line = 1.7, cex = 1.4, adj = 0.02)

## dev off
dev.off()


###############################################################################################
######################################   Figure 4G   ##########################################

###  Figure 4G)  Immunohistochemistry images for CD8 in DDIR subtypes - FOCUS --->>>

## Figure 4G was taken from QuPath (version 0.1.2)

# --- Sample ID --- DDIR --- High/Low  --- TMA   --- Co-ordinates (x,y)
# --- SC00332A  --- NEG  --- High      --- TMA7  --- (3,4)
# --- SC00265A  --- NEG  --- Low       --- TMA5  --- (9,2)
# --- SC00443A  --- POS  --- High      --- TMA10 --- (3,4)
# --- SC00228A  --- POS  --- Low       --- TMA4  --- (8,2)

###############################################################################################
######################################   Figure 4H   ##########################################

###  Figure 4H) IHC CD8 score vs DDIR  - barplot/scatterplot - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy
variable.1 <- "CD8_MEAN"

## Remove missing data
data.2 <- data.2[which(!is.na(data.2[[variable.1]])) ,]
range(data.2[[variable.1]])
dim(data.2)
## For log, onyl select > 0 
data.2 <- data.2[which(data.2[[variable.1]] > 0) ,]
dim(data.2)

## Log2 transformed
## NOTE: no need to add pseudocount (i.e. +1),
## as we excluded all 0 values from the data, only 1 sample
data.2[["CD8_MEAN_LOG"]] <- log2(data.2[[variable.1]])
range(data.2[["CD8_MEAN_LOG"]])

## Statistics:
pearson.1 <- cor.test(data.2[["DDIR_scores"]], data.2[["CD8_MEAN_LOG"]], data = data.2,
                      alternative = "two.sided", method = "pearson")
t.test.1 <-  t.test(data.2[["CD8_MEAN_LOG"]] ~ DDIR, data = data.2)
wilcox.1 <- wilcox.test(data.2[["CD8_MEAN_LOG"]] ~ DDIR, data = data.2)

pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")
pvalue.3 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001")

### ---

#tiff("Figure_4H.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(3.1,3.1,3.1,4.5)+1.5)

## Plot 1 (Scatterplot):
par(fig=c(0,1,0,0.85))
plot(data.2[["DDIR_scores"]] ~ data.2[["CD8_MEAN_LOG"]], data = data.2,
     pch = 19, ylab = "", xlab = "", cex = 1.6, cex.axis = 1.6,
     ylim = c(-0.3, 0.65))

## Add titles 1
title(ylab = paste0("DDIR Scores"), line =2.8, cex.lab = 1.6)
title(xlab = paste0("mean CD8 IHC scores [Log2]"), 
      line = 2.8, cex.lab = 1.6)

## Add stats 1
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(pvalue.3, side = 3,adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.2[["DDIR_scores"]] ~ data.2[["CD8_MEAN_LOG"]]), lwd = 2)

## Plot 2 (boxplot):
par(fig=c(0,1,0.65,1), new=TRUE)
boxplot(data.2[["CD8_MEAN_LOG"]] ~ data.2[["DDIR"]],  horizontal=TRUE, axes=FALSE, 
        col = c("blue", "red"))

## Add legend 2 (DDIR)
legend("topright", legend = c("POS", "NEG"), title = list("DDIR"),
       col = c("red", "blue"),pch = 15, title.adj = 0.78, y.intersp = 0.9,
       bty = "n", cex = 1.4, xpd = TRUE,inset = c(-0.22, -1.2), adj = c(0.1,0.5))

## Add titles 2
title(main =  paste0("CD8 IHC Score vs DDIR scores\n",cohort.1," trial cohort"),
      line = 1.4, cex.main = 1.8)
## Add stats 2
mtext(pvalue.1, side = 1, line = 0.5, cex = 1.4, adj = 0.02)
mtext(pvalue.2, side = 1, line = 1.7, cex = 1.4, adj = 0.02)

## dev off
dev.off()


###############################################################################################
######################################   Figure 5    ##########################################
###############################################################################################


######################################   Figure 5A   ##########################################

### Figure 5A)  Differential Expression Analysis VennDiagram --->>>
### --- Differential expression between DDIR Pos and DDIR Neg was performed within each cohort.
### --- Differentially expressed genes were then overlapped to identify common genes.


###############################################################################################
######################################   Figure 5B   ##########################################

###  Figure 5B)  Activated upstream regulators for 9 genes --->>>

ipa.9gene.filenm <- "FOCUS_DDRD_n361_degs_FDR_0_05_FC_1_5_9genes_upstream_regulators.txt"

ipa.9gene.file <- read.delim(file = file.path(data.folder,ipa.9gene.filenm),
                             header = T, stringsAsFactors = F, skip = 2)

data.1 <- ipa.9gene.file[ which( ipa.9gene.file$Predicted.Activation.State != " ") ,]
head(data.1)

## remove chemical or drug associated molecule type
data.2 <- data.1[ - which( grepl("chemical", data.1$Molecule.Type) |
                             grepl("drug", data.1$Molecule.Type) |
                             data.1$Activation.z.score < 0) ,]
## order based on activation score highest to lowest
data.2 <- data.2[order(data.2$Activation.z.score, decreasing = T) ,]
names(data.2)
data.3 <- data.2[,c(1,4,5,7,8)]


### --------------------------------------------------------------------------------------- ###
### Save the file (Supplementary Figure 8A)
#write.table(data.3, 
#            file = file.path(data.folder,"focus_9gene_activated_upstream_regulator.txt"), 
#            row.names = F, sep = "\t")
### --------------------------------------------------------------------------------------- ###

## Plot:
data.4 <- data.3[,c(1,3)] %>% .[order(.$Activation.z.score, decreasing = T),]  ### order
p1 <- ggplot(data.4, aes(x = reorder(Upstream.Regulator, Activation.z.score),
                         y = Activation.z.score, fill = Activation.z.score))+
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill = NA,colour="gray50"),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right",
        #legend.direction="horizontal",
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 0, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 0),
        axis.title.x = element_text(size = 14, angle = 0)) +
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  scale_x_discrete(limits=rev(data.4$Upstream.Regulator)) + 
  scale_fill_gradient2(name = paste0("Activation \nz-score"), 
                       high = "red", mid = "white", low = "blue")+
  ylab("Activation z-score") +
  ggtitle("Upstream Regulators of 9 genes")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                                  margin = ggplot2::margin(0.5,0.5,0.3,0.5, 'cm')))

#tiff("Figure_5B.tiff", height = 16, width = 16, units = 'cm', res = 800)
p1
dev.off()


###############################################################################################
######################################   Figure 5C   ##########################################

###  Figure 5C)  Almac DDIR vs 9 gene scores --->

data.1 <- focus.data

## Statistics:
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[["sum.score"]], data = data.1,
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001")       

### ---

#tiff("Figure_5C.tiff", height = 16, width = 18, units = 'cm', res = 800)

plot(data.1[["sum.score"]], data.1[["DDIR_scores"]],
     cex = 1.6, cex.axis = 1.6, pch = 19,
     xlab = "", ylab = "",
     ylim = c(-0.3, 0.65))

lines(lowess(data.1[["DDIR_scores"]] ~ data.1[["sum.score"]], f = 1/3), col = "dodgerblue", lwd = 3)

## Add titles
title(xlab = paste0("9 gene scores"), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("DDIR score vs 9 gene score\n", cohort.1, " trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(pvalue.1, side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add Legend
legend('topleft', legend = c("Lowess smoother"), 
       lty = 1, lwd = 3, col = c("dodgerblue"), bty = 'n',
       cex = 1.4, xpd = TRUE,inset = c(-0.017, 0.075))

## dev off
dev.off()


### --------------------------------------------------------------------------------------- ###

###############################################################################################
############################### --- SUPPLEMENTARY FIGURES --- #################################
###############################################################################################

### --------------------------------------------------------------------------------------- ###

###############################################################################################
######################################   Figure S3   ##########################################
###############################################################################################


#####################################   Figure S3A   ##########################################

### Figure s3A)  KEGG Homologous Recombination ssGSEA score - TRANSBIG --->>>

data.1 <- transbig.data
variable.1 <- 'KEGG_HOMOLOGOUS_RECOMBINATION'

## Statistics:
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001") 

### ---

#tiff("Supplementary_Figure_3A.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
plot(data.1[["DDIR_scores"]] ~ data.1[[variable.1]], data = data.1,
     ylab = "", xlab = "", pch = 19,
     cex = 1.6, cex.axis = 1.6)

## Add titles
title(xlab = paste0(variable.1), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("KEGG HR pathway"," ssGSEA scores\n",cohort.3," trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3,adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.1), side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold line [for BC: 0.37, for CRC: 0.109592]
abline(h=0.37, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]),lwd = 2, col = "blue")

## dev off
dev.off()


###############################################################################################
#####################################   Figure S3B   ##########################################


### Figure S3B)  REACTOME Fanconi Anaemia ssGSEA score vs DDIR score - TRANSBIG --->>>

data.1 <- transbig.data
variable.1 <- 'REACTOME_FANCONI_ANEMIA_PATHWAY'

## Statistics:
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001") 

### ---

#tiff("Supplementary_Figure_3B.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
plot(data.1[["DDIR_scores"]] ~ data.1[[variable.1]], data = data.1,
     ylab = "", xlab = "", pch = 19,
     cex = 1.6, cex.axis = 1.6)

## Add titles
title(xlab = paste0(variable.1), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("REACTOME FA pathway"," ssGSEA scores\n",cohort.3," trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3,adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.1), side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold line [for BC: 0.37, for CRC: 0.109592]
abline(h=0.37, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]),lwd = 2, col = "blue")

## dev off
dev.off()


###############################################################################################
##################################   Figure S3C and S3E  ######################################

### Figure S3C)  KEGG Homologous Recombination ssGSEA score - FOCUS --->>>
### Figure S3E)  KEGG Homologous Recombination ssGSEA score - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

variable.1 <- 'KEGG_HOMOLOGOUS_RECOMBINATION'

## We are including MSI for FOCUS/FOxTROT cohort ---
## Remove NA and replace with label 'NA' from MSI column
data.2 <- data.1 # make a copy
sum(is.na(data.2[["MSI"]]))
levels(data.2[["MSI"]]) <- c(levels(data.2[["MSI"]]),"NA")  # Add an extra level to the factor
data.2[["MSI"]][is.na(data.2[["MSI"]])] <- "NA"           # Change NA to "None"
table(data.2[["MSI"]])

## Set Colour for MSI status ---
col.1 <- c("orangered2", "black", "grey")
names(col.1) <- c("MSI", "MSS", "NA")
col.1 <- col.1[data.2[, 'MSI']]

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72

## Statistics:
pearson.1 <- cor.test(data.2[["DDIR_scores"]], data.2[[variable.1]], data = data.2,
                      alternative = "two.sided", method = "pearson")
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001") 


### ---

#tiff("Supplementary_Figure_3C.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Supplementary_Figure_3E.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(4.1,3.1,4.1,4.1)+1)

plot(data.2[["DDIR_scores"]] ~ data.2[[variable.1]], data = data.2,
     ylab = "", xlab = "", 
     ylim = c(y.min.1, y.max.1),
     cex = 1.6, cex.axis = 1.6, pch = 19, col = col.1)

## Add legend
legend("topright", legend = c("MSI", "MSS", "NA"),
       col = c("orangered2", "black", "grey"), pch = 19, bty = 'n',
       title = 'MSI', title.adj = 0.55,
       xpd = TRUE,inset = c(-0.23, -0.03), 
       y.intersp = 0.8, x.intersp = 1,
       adj = 0.12, cex = 1.6)

## Add titles
title(xlab = paste0(variable.1), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("KEGG HR pathway"," ssGSEA scores\n",
                    ## Select cohort.1 for FOCUS, cohort.2 for FOxTROT labels
                    cohort.1,
                    #cohort.2,
                    " trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.1), side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.2[["DDIR_scores"]] ~ data.2[[variable.1]]),
       lwd = 2, col = "blue")

## dev off
dev.off()

###############################################################################################
##################################   Figure S3D and S3F  ######################################

### Figure S3D)  REACTOME Fanconi Anaemia ssGSEA score vs DDIR score - FOCUS --->>>
### Figure S3F)  REACTOME Fanconi Anaemia ssGSEA score vs DDIR score - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

variable.1 <- 'REACTOME_FANCONI_ANEMIA_PATHWAY'

## We are including MSI for FOCUS/FOxTROT cohort ---
# Remove NA and replace with label 'NA' from MSI column
data.2 <- data.1 # make a copy
sum(is.na(data.2[["MSI"]]))
levels(data.2[["MSI"]]) <- c(levels(data.2[["MSI"]]),"NA")  # Add an extra level to the factor
data.2[["MSI"]][is.na(data.2[["MSI"]])] <- "NA"           # Change NA to "None"
table(data.2[["MSI"]])

## Set Colour for MSI status ---
col.1 <- c("orangered2", "black", "grey")
names(col.1) <- c("MSI", "MSS", "NA")
col.1 <- col.1[data.2[, 'MSI']]

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72

## Statistics:
pearson.1 <- cor.test(data.2[["DDIR_scores"]], data.2[[variable.1]], data = data.2,
                      alternative = "two.sided", method = "pearson")
pvalue.1 <- ifelse(pearson.1$p.value > 0.0001,
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001") 

### ---

#tiff("Supplementary_Figure_3D.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Supplementary_Figure_3F.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(4.1,3.1,4.1,4.1)+1)

plot(data.2[["DDIR_scores"]] ~ data.2[[variable.1]], data = data.2,
     ylab = "", xlab = "", 
     ylim = c(y.min.1, y.max.1),
     cex = 1.6, cex.axis = 1.6, pch = 19, col = col.1)

## Add legend
legend("topright", legend = c("MSI", "MSS", "NA"),
       col = c("orangered2", "black", "grey"), pch = 19, bty = 'n',
       title = 'MSI', title.adj = 0.55,
       xpd = TRUE,inset = c(-0.23, -0.03), 
       y.intersp = 0.8, x.intersp = 1,
       adj = 0.12, cex = 1.6)

## Add titles
title(xlab = paste0(variable.1), line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("REACTOME FA pathway"," ssGSEA scores\n",
                    ## Select cohort.1 for FOCUS, cohort.2 for FOxTROT labels
                    cohort.1,
                    #cohort.2,
                    " trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.1), side = 3, adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (CRC = 0.109592)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.2[["DDIR_scores"]] ~ data.2[[variable.1]]),
       lwd = 2, col = "blue")

## dev off
dev.off()


###############################################################################################
######################################   Figure S4   ##########################################
###############################################################################################

###############################################################################################
#####################################   Figure S4A   ##########################################

###  Figure S4A)  Immunohistochemistry images for cGAS in DDIR subtypes - FOCUS --->>>

## Figure S4A was taken from QuPath (version 0.1.2)

# --- Sample ID --- DDIR --- High/Low  --- TMA   --- Co-ordinates (x,y)
# --- SC00150A  --- NEG  --- High      --- TMA1  --- (10,10)
# --- SC00480A  --- NEG  --- Low       --- TMA11 --- (3,7)
# --- SC00301A  --- POS  --- High      --- TMA6  --- (4,5)
# --- SC00202A  --- POS  --- Low       --- TMA3  --- (10,4)


###############################################################################################
#####################################   Figure S4B   ##########################################

###  Figure S4B)  Immunohistochemistry images for STING in DDIR subtypes - FOCUS --->>>

## Figure S4B was taken from QuPath (version 0.1.2)

# --- Sample ID --- DDIR --- High/Low  --- TMA   --- Co-ordinates (x,y)
# --- SC00376A  --- NEG  --- High      --- TMA8  --- (8,4)
# --- SC00478A  --- NEG  --- Low       --- TMA11 --- (13,4)
# --- SC00393A  --- POS  --- High      --- TMA8  --- (8,10)
# --- SC00173A  --- POS  --- Low       --- TMA2  --- (13,7)


###############################################################################################
#####################################   Figure S4C   ##########################################

###  Figure S4C)  IHC cGAS scores in DDIR - boxplot - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy

## Examine number of missing data
sum(is.na(data.2$Tumour_PerCent_cGAS_Positive_Cells))
data.2 <- data.2[which(!is.na(data.2$Tumour_PerCent_cGAS_Positive_Cells)) ,]

## Statistics:
t.test.1 <- t.test(data.2[c(which(data.2$DDIR == "DDIR NEG")) ,]$Tumour_PerCent_cGAS_Positive_Cells,
                   data.2[c(which(data.2$DDIR == "DDIR POS")) ,]$Tumour_PerCent_cGAS_Positive_Cells,
                   data = data.2)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")

wilcox.1 <- wilcox.test(data.2[c(which(data.2$DDIR == "DDIR NEG")) ,]$Tumour_PerCent_cGAS_Positive_Cells,
                        data.2[c(which(data.2$DDIR == "DDIR POS")) ,]$Tumour_PerCent_cGAS_Positive_Cells,
                        data = data.2)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")


## Set Colour for DDIR status ---
col.1 <- c("red", "blue")
names(col.1) <- c("DDIR POS", "DDIR NEG")
col.1 <- col.1[data.2[, 'DDIR']]

### ---

#tiff("Supplementary_Figure_4C.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
boxplot(data.2$Tumour_PerCent_cGAS_Positive_Cells ~ data.2$DDIR, data = data.2,
        ylim = c(0,100), xlab = "", ylab = "",
        cex.axis = 1.6,cex.lab = 1.6)
beeswarm(data.2$Tumour_PerCent_cGAS_Positive_Cells ~ data.2$DDIR, data = data.2,
         add = TRUE, pch = 19, pwcol = col.1, 
         method = "swarm", corral = "wrap",cex = 1.5, spacing = 1.2)

## Add titles
title(ylab = paste0("% of total cGAS-positive cells per tumour sample"), 
      line = 2.8, cex.lab = 1.6)
title(main = paste0("cGAS IHC positive cells in each samples\n", cohort.1," trial cohort"),
      line = 1.2, cex.main = 1.8)

## Add stats
mtext(pvalue.1, side = 1, line = 2.6, cex = 1.6)
mtext(pvalue.2, side = 1, line = 4, cex = 1.6)

## dev off
dev.off()

###############################################################################################
#####################################   Figure S4D   ##########################################

###  Figure S4D)  IHC cGAS scores against DDIR scores - scatterplot - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy

## Examine number of missing data
sum(is.na(data.2$Tumour_PerCent_cGAS_Positive_Cells))
data.2 <- data.2[which(!is.na(data.2$Tumour_PerCent_cGAS_Positive_Cells)) ,]

## Statistics:
pearson.1 <- cor.test(data.2$Tumour_PerCent_cGAS_Positive_Cells, 
                      data.2[["DDIR_scores"]], data = data.2, 
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ", round(pearson.1$p.value,4)),
                   paste0("P (two-tailed) < 0.0001"))

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72

### ---

#tiff("Supplementary_Figure_4D.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
plot(data.2[["DDIR_scores"]] ~ data.2$Tumour_PerCent_cGAS_Positive_Cells, 
     data = data.2,xlab = "", ylab = "", ylim = c(y.min.1, y.max.1),
     pch = 19, cex = 1.6, cex.axis = 1.6)

## Add titles
title(xlab = paste0("% of cGAS IHC positive cells per sample"),line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("% of total cGAS IHC positive cells\n",cohort.1," trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate,4)),
      cex = 1.4, adj = 0.02, line = -1.4)
mtext(pvalue.1, cex = 1.4, adj = 0.02, line = -2.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.2[["DDIR_scores"]] ~ data.2$Tumour_PerCent_cGAS_Positive_Cells),
       col = "blue", lwd=2)

## dev off
dev.off()

###############################################################################################
#####################################   Figure S4E   ##########################################

###  Figure S4E)  IHC STING scores in DDIR - boxplot - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy

## Examine number of missing data
sum(is.na(data.2$Tumour_PerCent_STING_Positive_Cells))
data.2 <- data.2[which(!is.na(data.2$Tumour_PerCent_STING_Positive_Cells)) ,]

## Statistics:
t.test.1 <- t.test(data.2[c(which(data.2$DDIR == "DDIR NEG")) ,]$Tumour_PerCent_STING_Positive_Cells,
                   data.2[c(which(data.2$DDIR == "DDIR POS")) ,]$Tumour_PerCent_STING_Positive_Cells,
                   data = data.2)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")

wilcox.1 <- wilcox.test(data.2[c(which(data.2$DDIR == "DDIR NEG")) ,]$Tumour_PerCent_STING_Positive_Cells,
                        data.2[c(which(data.2$DDIR == "DDIR POS")) ,]$Tumour_PerCent_STING_Positive_Cells,
                        data = data.2)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")


## Set Colour for DDIR status ---
col.1 <- c("red", "blue")
names(col.1) <- c("DDIR POS", "DDIR NEG")
col.1 <- col.1[data.2[, 'DDIR']]

### ---

#tiff("Supplementary_Figure_4E.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
boxplot(data.2$Tumour_PerCent_STING_Positive_Cells ~ data.2$DDIR, data = data.2,
        ylim = c(0,100), xlab = "", ylab = "",
        cex.axis = 1.6,cex.lab = 1.6)
beeswarm(data.2$Tumour_PerCent_STING_Positive_Cells ~ data.2$DDIR, data = data.2,
         add = TRUE, pch = 19, pwcol = col.1, 
         method = "swarm", corral = "wrap",cex = 1.5, spacing = 1.2)

## Add titles
title(ylab = paste0("% of total STING-positive cells per tumour sample"), 
      line = 2.8, cex.lab = 1.6)
title(main = paste0("STING IHC positive cells in each samples\n", cohort.1," trial cohort"),
      line = 1.2, cex.main = 1.8)

## Add stats
mtext(pvalue.1, side = 1, line = 2.6, cex = 1.6)
mtext(pvalue.2, side = 1, line = 4, cex = 1.6)

## dev off
dev.off()

###############################################################################################
#####################################   Figure S4F   ##########################################

###  Figure S4F)  IHC STING scores against DDIR scores - scatterplot - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy

## Examine number of missing data
sum(is.na(data.2$Tumour_PerCent_STING_Positive_Cells))
data.2 <- data.2[which(!is.na(data.2$Tumour_PerCent_STING_Positive_Cells)) ,]

## Statistics:
pearson.1 <- cor.test(data.2$Tumour_PerCent_STING_Positive_Cells, 
                      data.2[["DDIR_scores"]], data = data.2, 
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ", round(pearson.1$p.value,4)),
                   paste0("P (two-tailed) < 0.0001"))

## Determine DDIR (y-axis) range
## Range if y-axis is DDIR score
## NOTE: to keep consistent between FOCUS and FOxTROT, lowest minimun value, 
## and highest maximum value (+/-0.1) of DDIR scores between the two cohort 
## will be used as y-axis.
# min value from FOCUS and max value from FOxTROT DDIR scores (-/+0.1) will be used.
y.min.1 <- -0.36
y.max.1 <- 0.72

### ---

#tiff("Supplementary_Figure_4F.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
plot(data.2[["DDIR_scores"]] ~ data.2$Tumour_PerCent_STING_Positive_Cells, 
     data = data.2, xlab = "", ylab = "", ylim = c(y.min.1, y.max.1),
     pch = 19, cex = 1.6, cex.axis = 1.6)

## Add titles
title(xlab = paste0("% of STING IHC positive cells per sample"),line = 2.8, cex.lab = 1.6)
title(ylab = paste0("DDIR Scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("% of total STING IHC positive cells\n",cohort.1," trial cohort"),
      line = 1, cex.main = 1.8)

## Add stats
mtext(paste0("Pearson r = ", round(pearson.1$estimate,4)),
      cex = 1.4, adj = 0.02, line = -1.4)
mtext(pvalue.1, cex = 1.4, adj = 0.02, line = -2.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.2[["DDIR_scores"]] ~ data.2$Tumour_PerCent_STING_Positive_Cells),
       col = "blue", lwd=2)

## dev off
dev.off()


###############################################################################################
######################################   Figure S5   ##########################################
###############################################################################################

###############################################################################################
###############################    Figure S5A and S5B   #######################################

###  Figure S5A)  CMS vs DDIR - barplot - FOCUS --->>>
###  Figure S5B)  CMS vs DDIR - barplot - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

table.1 <- table(data.1[["CMS"]], data.1[["DDIR"]])
prop.1 <- prop.table(table.1, 2)
prop.long.1 <- as.data.frame(prop.1)
names(prop.long.1) <- c("CMS", "DDIR", "Proportion")

## Statistics:
fisher.1 <- fisher.test(table.1)
pvalue.1 <- ifelse(fisher.1$p.value > 0.0001, 
                   paste0("Fisher's Exact test, p = ",round(fisher.1$p.value, 4)),
                   "Fisher's Exact, p < 0.0001")


### ---

#tiff("Supplementary_Figure_5A.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Supplementary_Figure_5B.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
ggplot(prop.long.1) + 
  geom_bar(aes(x = DDIR, y = Proportion, fill = CMS),data = prop.long.1, 
           stat = "identity", width = 0.7,
           position = position_stack(reverse = TRUE))+
  ## Add colour and legend
  scale_fill_manual(values = c("#E18C28", "#0B5F9C", "#C46294", "#129063", "#A4A4A4"),
                    guide=guide_legend(reverse=TRUE, title = paste0('CMS'))) +
  ## Add limit to y-axis
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  ## Add titles
  xlab(" ") +
  ylab("Sample Proportions") +
  labs(title = paste0("CMS subtypes in DDIR\n",
                      ## Select cohort.1 for FOCUS, and cohort.2 for FOxTROT
                      cohort.1,
                      #cohort.2,
                      " trial cohort")) +
  
  
  ## Theme for main title
  theme(plot.title = element_text(face = 'bold', size = 22,
                                  margin = margin(15,0,0,0), 
                                  vjust = 4.3, hjust = 0.5, lineheight = 1)) +
  
  ## Add p-value
  labs(caption = pvalue.1)+
  theme(plot.caption = element_text(size = 20, 
                                    margin = margin(15,0,0,0), hjust = 0.5)) +
  
  ## Plot background theme
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, vjust = 3),
        axis.text.y = element_text(size = 18, colour = "black", angle = 90, hjust = 0.5),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white'),
        
        ## legend
        legend.background = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = 'bold'),
        legend.title.align = 0,
        legend.position = c(0.95,0.765),
        legend.key.height = unit(0.7, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.margin = margin(-25,0,0,35),
        plot.margin = unit(c(5,50,5,5),'mm'))


## dev off
dev.off()


###############################################################################################
###############################    Figure S5C and S5D   #######################################

###  Figure S5C)  CRIS vs DDIR - barplot - FOCUS --->>>
###  Figure S5D)  CRIS vs DDIR - barplot - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

data.1[["CRIS"]] <- plyr::revalue(data.1[["CRIS"]], c("Unclassified"= "UNK"))

table.1 <- table(data.1[["CRIS"]], data.1[["DDIR"]])
prop.1 <- prop.table(table.1, 2)
prop.long.1 <- as.data.frame(prop.1)
names(prop.long.1) <- c("CRIS", "DDIR", "Proportion")

## Statistics:
fisher.1 <- fisher.test(table.1)
pvalue.1 <- ifelse(fisher.1$p.value > 0.0001, 
                   paste0("Fisher's Exact test, p = ",round(fisher.1$p.value, 4)),
                   "Fisher's Exact, p < 0.0001")


### ---

#tiff("Supplementary_Figure_5C.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Supplementary_Figure_5D.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
ggplot(prop.long.1) + 
  geom_bar(aes(x = DDIR, y = Proportion, fill = CRIS),data = prop.long.1, 
           stat = "identity", width = 0.7,
           position = position_stack(reverse = TRUE))+
  ## Add colour and legend
  scale_fill_manual(values = c("#FA230A", "#A00013", "#111C54", "#0F6F22", "#1AB897", "#A4A4A4"),
                    guide=guide_legend(reverse=TRUE, title = paste0('CRIS'))) +
  ## Add limit to y-axis
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  ## Add titles
  xlab(" ") +
  ylab("Sample Proportions") +
  labs(title = paste0("CRIS subtypes in DDIR\n",
                      ## Select cohort.1 for FOCUS, and cohort.2 for FOxTROT
                      cohort.1,
                      #cohort.2,
                      " trial cohort")) +
  
  
  ## Theme for main title
  theme(plot.title = element_text(face = 'bold', size = 22,
                                  margin = margin(15,0,0,0), 
                                  vjust = 4.3, hjust = 0.5, lineheight = 1)) +
  
  ## Add p-value
  labs(caption = pvalue.1)+
  theme(plot.caption = element_text(size = 20, 
                                    margin = margin(15,0,0,0), hjust = 0.5)) +
  
  ## Plot background theme
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, vjust = 3),
        axis.text.y = element_text(size = 18, colour = "black", angle = 90, hjust = 0.5),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'white'),
        
        ## legend
        legend.background = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = 'bold'),
        legend.title.align = 0,
        legend.position = c(0.95,0.73),
        legend.key.height = unit(0.7, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.margin = margin(-25,0,0,35),
        plot.margin = unit(c(5,50,5,5),'mm'))


## dev off
dev.off()


###############################################################################################
#####################################   Figure S5E   ##########################################

###  Figure S5E)  CRIS vs DDIR scores - boxplot - FOCUS --->>>

data.1 <- focus.data
data.1[["CRIS"]] <- plyr::revalue(data.1[["CRIS"]], c("Unclassified"= "UNK"))

## Statistics:
kruskal.1 <- kruskal.test(data.1[["DDIR_scores"]] ~ as.factor(data.1[["CRIS"]]), data = data.1)
pvalue.1 <- ifelse(kruskal.1$p.value > 0.0001, 
                   paste0("Kruskal-Wallis, p = ", round(kruskal.1$p.value, 4)),
                   paste0("Kruskal-Wallis, p < 0.0001"))

## Tukey's HSD test
## form a model (one-way ANOVA)
model.1 <- lm(DDIR_scores ~ CRIS, data = data.1)
summary(model.1)
ANOVA.1 <- aov(model.1, data = data.1)
summary(ANOVA.1)
## posthoc tukey
posthoc.1 <- TukeyHSD(ANOVA.1, 'CRIS', conf.level = 0.95)
Tukey.result.1 <- posthoc.1$CRIS
## NOTE: none of the stats are significant, so will not be included in the figure


### ---

#tiff("Supplementary_Figure_5E.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(2.1,4.1,4,2.1))

boxplot(data.1[["DDIR_scores"]] ~ data.1[["CRIS"]], data = data.1,
        ylim = c(-0.30, 1), cex.axis = 1.2, outline = FALSE, xlab = "")
beeswarm(data.1[["DDIR_scores"]] ~ data.1[["CRIS"]], data = data.1,
         add = TRUE, pch = 19, cex = 1.6, spacing = 1.2,
         method = "swarm",corral="wrap",
         col = c("#FA230A", "#A00013", "#111C54", "#0F6F22", "#1AB897", "#A4A4A4"))

## Add titles
title(main = paste0("DDIR score by CRIS subtypes\n",cohort.1," trial cohort"),
      line = 1, cex.main = 1.8)
title(ylab = paste0("DDIR Scores"), line = 2.6, cex.lab = 1.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add stats
mtext(pvalue.1, side = 3, line = -1.4, adj = 0.98, cex = 1.2)

## dev off
dev.off()

###############################################################################################
#####################################   Figure S5F   ##########################################

###  Figure S5F)  CRIS vs DDIR scores - boxplot - FOxTROT --->>>

data.1 <- foxtrot.data
data.1[["CRIS"]] <- plyr::revalue(data.1[["CRIS"]], c("Unclassified"= "UNK"))

## Statistics:
kruskal.1 <- kruskal.test(data.1[["DDIR_scores"]] ~ as.factor(data.1[["CRIS"]]), data = data.1)
pvalue.1 <- ifelse(kruskal.1$p.value > 0.0001, 
                   paste0("Kruskal-Wallis, p = ", round(kruskal.1$p.value, 4)),
                   paste0("Kruskal-Wallis, p < 0.0001"))

## Tukey's HSD test
## form a model (one-way ANOVA)
model.1 <- lm(DDIR_scores ~ CRIS, data = data.1)
summary(model.1)
ANOVA.1 <- aov(model.1, data = data.1)
summary(ANOVA.1)
## posthoc tukey
posthoc.1 <- TukeyHSD(ANOVA.1, 'CRIS', conf.level = 0.95)
Tukey.result.1 <- posthoc.1$CRIS
## NOTE: none of the stats are significant, so will not be included in the figure


### ---

#tiff("Supplementary_Figure_5F.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(2.1,4.1,4,2.1))

boxplot(data.1[["DDIR_scores"]] ~ data.1[["CRIS"]], data = data.1,
        ylim = c(-0.30, 1), cex.axis = 1.2, outline = FALSE, xlab = "")
beeswarm(data.1[["DDIR_scores"]] ~ data.1[["CRIS"]], data = data.1,
         add = TRUE, pch = 19, cex = 1.6, spacing = 1.2,
         method = "swarm",corral="wrap",
         col = c("#FA230A", "#A00013", "#111C54", "#0F6F22", "#1AB897", "#A4A4A4"))

## Add titles
title(main = paste0("DDIR score by CRIS subtypes\n",cohort.2," trial cohort"),
      line = 1, cex.main = 1.8)
title(ylab = paste0("DDIR Scores"), line = 2.6, cex.lab = 1.6)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add stats
mtext(pvalue.1, side = 3, line = -1.4, adj = 0.98, cex = 1.2)

## dev off
dev.off()


###############################################################################################
####################################    Figure S5H    #########################################

### Figure S5H)  Mutational Burden vs DDIR - FOCUS --->>>

data.1 <- focus.data
data.2 <- data.1 # make a copy

data.2 <- data.2[which(!is.na(data.2$Num_Muts)) ,]

## Statistics:
t.test.1 <- t.test(Num_Muts ~ DDIR, data = data.2)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
wilcox.1 <- wilcox.test(Num_Muts ~ DDIR, data = data.2)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")

## determine the y-axis range
min.1 <- min(data.2$Num_Muts) -0.1
max.1 <- max(data.2$Num_Muts) +0.1

## Set Colour for MSI status ---
col.1 <- c("darkorange", "blue4", "grey")
names(col.1) <- c("MSI", "MSS", "NA")
col.1 <- col.1[data.2[, 'MSI']]

### ---

#tiff("Supplementary_Figure_5G.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(5.1,4.5,4.1,2.1)+0.5)

boxplot(Num_Muts ~ DDIR, data = data.2, 
        outline = FALSE, ylim = c(min.1, max.1),
        cex.axis = 1.6,cex.lab = 1.6)
beeswarm(Num_Muts ~ DDIR, data = data.2, 
         method = "swarm", corral = "wrap", pch = 19,
         cex = 1.5, spacing = 1.2,add = TRUE, pwcol = col.1)

## Add Legend
legend("topright", legend = c("MSI", "MSS"),col = c("darkorange", "blue4"),
       pch = 19, bty = "n", cex = 1.6,
       title = 'MSI', title.adj = 0.67,
       y.intersp = 0.8, x.intersp = 1)

## Add titles
title(main = paste0("Mutation burden in DDIR\n",cohort.1," trial cohort"), cex.main = 1.8)
title(ylab = paste0("Number of mutations"), line = 2.8, cex.lab = 1.6)

## Add stats
mtext(pvalue.1, side = 1, line = 2.8, cex = 1.6)
mtext(pvalue.2, side = 1, line = 4.2, cex = 1.6)

## dev off
dev.off()


###############################################################################################
####################################    Figure S5H    #########################################

### Figure S5H)  Mutational Burden vs DDIR - FOxTROT --->>>

data.1 <- foxtrot.data
data.2 <- data.1 # make a copy

data.2 <- data.2[- c( which(data.2$Total_coding_muts == "NULL")) ,]
data.2$Total_coding_muts <- as.numeric(data.2$Total_coding_muts)

## Statistics:
t.test.1 <- t.test(Total_coding_muts ~ DDIR, data = data.2)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
wilcox.1 <- wilcox.test(Total_coding_muts ~ DDIR, data = data.2)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")

# determine the y-axis range
min.1 <- min(data.2$Total_coding_muts) -0.1
max.1 <- max(data.2$Total_coding_muts) +0.1

## Set Colour for MSI status ---
col.1 <- c("darkorange", "blue4", "grey")
names(col.1) <- c("MSI", "MSS", "NA")
col.1 <- col.1[data.2[, 'MSI']]

### ---

#tiff("Supplementary_Figure_5H.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(5.1,4.5,4.1,2.1)+0.5)

boxplot(Total_coding_muts ~ DDIR, data = data.2, 
        outline = FALSE, ylim = c(min.1, max.1),
        cex.axis = 1.6,cex.lab = 1.6)
beeswarm(Total_coding_muts ~ DDIR, data = data.2, 
         method = "swarm", corral = "wrap", pch = 19,
         cex = 1.5, spacing = 1.2,add = TRUE, pwcol = col.1)

## Add Legend
legend("topright", legend = c("MSI", "MSS"),col = c("darkorange", "blue4"),
       pch = 19, bty = "n", cex = 1.6,
       title = 'MSI', title.adj = 0.67,
       y.intersp = 0.8, x.intersp = 1)

## Add titles
title(main = paste0("Mutation burden in DDIR\n",cohort.2," trial cohort"), cex.main = 1.8)
title(ylab = paste0("Number of mutations"), line = 2.8, cex.lab = 1.6)

## Add stats
mtext(pvalue.1, side = 1, line = 2.8, cex = 1.6)
mtext(pvalue.2, side = 1, line = 4.2, cex = 1.6)

## dev off
dev.off()


###############################################################################################
######################################   Figure S6    #########################################
###############################################################################################

###############################################################################################
####################################    Figure S6A    #########################################

###  Figure S6A) GSEA Hallmarks DDIR POS vs DDIR NEG - dotplot - FOCUS --->>>

pos.1 <- focus.gsea.pos
neg.1 <- focus.gsea.neg

## Chose the desired columns for pos and neg gsea df
pos.1 <- data.frame(pos.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
pos.1 <- pos.1[order(pos.1[['FDR.q.val']], decreasing = FALSE),]# rearrange based on "FDR q-val"
head(pos.1)
neg.1 <- data.frame(neg.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
neg.1 <- neg.1[order(neg.1[['FDR.q.val']], decreasing = FALSE),]
head(neg.1)

## Remove the string characters (i.e. "HALLMARK_", and "_")
pos.1[,1] <- gsub("HALLMARK_|_", " ", pos.1[,1])
neg.1[,1] <- gsub("HALLMARK_|_", " ", neg.1[,1])

## Order them based on column 'Name' and merge
pos.1 <- pos.1[order(pos.1$NAME, decreasing = FALSE),]
neg.1 <- neg.1[order(neg.1$NAME, decreasing = FALSE),]
gsea.pos.neg <- rbind(pos.1, neg.1)
gsea.pos.neg.nes.ordered <- gsea.pos.neg[order(gsea.pos.neg[['NES']], decreasing = TRUE) ,]
## Select FDR significant rows (<0.25)
fdr.sig.1 <- subset(gsea.pos.neg.nes.ordered, gsea.pos.neg.nes.ordered[['FDR.q.val']] < 0.25)

### ---

## Plot:
p1 <- ggplot(fdr.sig.1, aes(NES, NAME))+ 
  geom_point(aes(colour=FDR.q.val, size=SIZE)) +
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  expand_limits(x=c(-2.25,2.25))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2)) +
  scale_y_discrete(limits=rev(fdr.sig.1[["NAME"]])) + 
  labs(list(colour = "FDR q-value", size = "Gene set size")) +
  xlab("Normalized enrichment score") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  ggtitle(paste0("Gene Set Enrichment Analysis (FDR significant)\n",cohort.1," trial cohort")) # change with cohort

p1a <- p1 + theme(panel.background=element_blank(),
                  panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
                  panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
                  panel.border=element_rect(fill = NA, linetype = 'solid', colour="gray20"),
                  axis.title.y=element_blank(),
                  axis.text.x = element_text(size = 10, angle = 0, face = "bold"),
                  axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
                  axis.line =  element_blank(),
                  axis.ticks = element_line(colour = "gray20"),
                  legend.position="right",
                  legend.direction="vertical",
                  legend.key = element_rect(fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5, margin = ggplot2::margin(0,0.5,0.8,0.5, 'cm')))



## Add DDRD labels on top of each panel
t1 <- text_grob(label = "DDIR NEG",color = 'blue', face = 'bold')
t2 <- text_grob(label = "DDIR POS",color = 'red', face = 'bold')

## Add text on graph (may need position adjustment)
## FOCUS: p2a -> ymax = 46, xmin = -2.4, xmax = -1.25, p2b -> xmin = 1.7, xmax = 2.1
## FOxTROT: p2a -> ymax = 18.2, xmin = -1.7, xmax = -1, p2b -> xmin = 1.25, xmax = 1.6
## TRANSBG: p2a -> ymax = 58.5, xmin = -2.3, xmax = -1.2, p2b -> xmin = 1.25, xmax = 2.3
p1b <- p1a +
  annotation_custom(t1, xmin = -2.30, xmax = -1.20, ymax = 46)+
  coord_cartesian(clip = 'off')
p1c <- p1b +
  annotation_custom(t2, xmin = 1.77, xmax = 1.80, ymax = 46)+ 
  coord_cartesian(clip = 'off')

#tiff("Supplementary_Figure_6A.tiff", height = 16, width = 18, units = 'cm', res = 800)
p1c
dev.off()

###############################################################################################
####################################    Figure S5B    #########################################

###  Figure S5B) GSEA Hallmarks DDIR POS vs DDIR NEG - dotplot - FOxTROT --->>>

pos.1 <- foxtrot.gsea.pos
neg.1 <- foxtrot.gsea.neg

## Chose the desired columns for pos and neg gsea df
pos.1 <- data.frame(pos.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
pos.1 <- pos.1[order(pos.1[['FDR.q.val']], decreasing = FALSE),]# rearrange based on "FDR q-val"
head(pos.1)
neg.1 <- data.frame(neg.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
neg.1 <- neg.1[order(neg.1[['FDR.q.val']], decreasing = FALSE),]
head(neg.1)

## Remove the string characters (i.e. "HALLMARK_", and "_")
pos.1[,1] <- gsub("HALLMARK_|_", " ", pos.1[,1])
neg.1[,1] <- gsub("HALLMARK_|_", " ", neg.1[,1])

## Order them based on column 'Name' and merge
pos.1 <- pos.1[order(pos.1$NAME, decreasing = FALSE),]
neg.1 <- neg.1[order(neg.1$NAME, decreasing = FALSE),]
gsea.pos.neg <- rbind(pos.1, neg.1)
gsea.pos.neg.nes.ordered <- gsea.pos.neg[order(gsea.pos.neg[['NES']], decreasing = TRUE) ,]
## Select FDR significant rows (<0.25)
fdr.sig.1 <- subset(gsea.pos.neg.nes.ordered, gsea.pos.neg.nes.ordered[['FDR.q.val']] < 0.25)

### ---

## Plot:
p1 <- ggplot(fdr.sig.1, aes(NES, NAME))+ 
  geom_point(aes(colour=FDR.q.val, size=SIZE)) +
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  expand_limits(x=c(-2.25,2.25))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2)) +
  scale_y_discrete(limits=rev(fdr.sig.1[["NAME"]])) + 
  labs(list(colour = "FDR q-value", size = "Gene set size")) +
  xlab("Normalized enrichment score") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  ggtitle(paste0("Gene Set Enrichment Analysis (FDR significant)\n",cohort.2," trial cohort")) # change with cohort

p1a <- p1 + theme(panel.background=element_blank(),
                  panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
                  panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
                  panel.border=element_rect(fill = NA,linetype = 'solid', colour="gray20"),
                  axis.title.y=element_blank(),
                  axis.text.x = element_text(size = 10, angle = 0, face = "bold"),
                  axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
                  axis.line =  element_blank(),
                  axis.ticks = element_line(colour = "gray20"),
                  legend.position="right",
                  legend.direction="vertical",
                  legend.key = element_rect(fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5, margin = ggplot2::margin(0,0.5,0.8,0.5, 'cm')))


## Add DDIR labels on top of each panel
t1 <- text_grob(label = "DDIR NEG",color = 'blue', face = 'bold')
t2 <- text_grob(label = "DDIR POS",color = 'red', face = 'bold')

## Add text on graph (may need position adjustment)
## FOCUS: p2a -> ymax = 46, xmin = -2.4, xmax = -1.25, p2b -> xmin = 1.7, xmax = 2.1
## FOxTROT: p2a -> ymax = 18.2, xmin = -1.7, xmax = -1, p2b -> xmin = 1.25, xmax = 1.6
## TRANSBG: p2a -> ymax = 58.5, xmin = -2.3, xmax = -1.2, p2b -> xmin = 1.25, xmax = 2.3
p1b <- p1a +
  annotation_custom(t1, xmin = -2.05, xmax = -1.7, ymax = 17.3)+
  coord_cartesian(clip = 'off')
p1c <- p1b +
  annotation_custom(t2, xmin = 1.7, xmax = 2.05, ymax = 17.3)+ 
  coord_cartesian(clip = 'off')

#tiff("Supplementary_Figure_3B.tiff", height = 12, width = 18, units = 'cm', res = 800)
p1c
dev.off()

###############################################################################################
####################################    Figure S6C    #########################################

###  Figure S6C) GSEA Hallmarks DDIR POS vs DDIR NEG - dotplot - TRANSBIG --->>>

pos.1 <- transbig.gsea.pos
neg.1 <- transbig.gsea.neg

## Chose the desired columns for pos and neg gsea df
pos.1 <- data.frame(pos.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
pos.1 <- pos.1[order(pos.1[['FDR.q.val']], decreasing = FALSE),]# rearrange based on "FDR q-val"
head(pos.1)
neg.1 <- data.frame(neg.1[,c("NAME", "SIZE", "NES", "FDR.q.val")])
neg.1 <- neg.1[order(neg.1[['FDR.q.val']], decreasing = FALSE),]
head(neg.1)

## Remove the string characters (i.e. "HALLMARK_", and "_")
pos.1[,1] <- gsub("HALLMARK_|_", " ", pos.1[,1])
neg.1[,1] <- gsub("HALLMARK_|_", " ", neg.1[,1])

## Order them based on column 'Name' and merge
pos.1 <- pos.1[order(pos.1$NAME, decreasing = FALSE),]
neg.1 <- neg.1[order(neg.1$NAME, decreasing = FALSE),]
gsea.pos.neg <- rbind(pos.1, neg.1)
gsea.pos.neg.nes.ordered <- gsea.pos.neg[order(gsea.pos.neg[['NES']], decreasing = TRUE) ,]
## Select FDR significant rows (<0.25)
fdr.sig.1 <- subset(gsea.pos.neg.nes.ordered, gsea.pos.neg.nes.ordered[['FDR.q.val']] < 0.25)

### ---

## Plot:
p1 <- ggplot(fdr.sig.1, aes(NES, NAME))+ 
  geom_point(aes(colour=FDR.q.val, size=SIZE)) +
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  expand_limits(x=c(-2.25,2.25))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2)) +
  scale_y_discrete(limits=rev(fdr.sig.1[["NAME"]])) + 
  labs(list(colour = "FDR q-value", size = "Gene set size")) +
  xlab("Normalized enrichment score") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  ggtitle(paste0("Gene Set Enrichment Analysis (FDR significant)\n",cohort.3," trial cohort")) # change with cohort

p1a <- p1 + theme(panel.background=element_blank(),
                  panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
                  panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
                  panel.border=element_rect(fill = NA,colour="gray20", linetype='solid'),
                  axis.title.y=element_blank(),
                  axis.text.x = element_text(size = 10, angle = 0, face = "bold"),
                  axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
                  axis.line =  element_blank(),
                  axis.ticks = element_line(colour = "gray20"),
                  legend.position="right",
                  legend.direction="vertical",
                  legend.key = element_rect(fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5, margin = ggplot2::margin(0,0.5,0.8,0.5, 'cm')))


## Add DDIR labels on top of each panel
t1 <- text_grob(label = "DDIR NEG",color = 'blue', face = 'bold')
t2 <- text_grob(label = "DDIR POS",color = 'red', face = 'bold')

## Add text on graph (may need position adjustment)
## FOCUS: p2a -> ymax = 46, xmin = -2.4, xmax = -1.25, p2b -> xmin = 1.7, xmax = 2.1
## FOxTROT: p2a -> ymax = 18.2, xmin = -1.7, xmax = -1, p2b -> xmin = 1.25, xmax = 1.6
## TRANSBG: p2a -> ymax = 58.5, xmin = -2.3, xmax = -1.2, p2b -> xmin = 1.25, xmax = 2.3
p1b <- p1a +
  annotation_custom(t1, xmin = -2.16, xmax = -1.32, ymax = 58.5)+
  coord_cartesian(clip = 'off')
p1c <- p1b +
  annotation_custom(t2, xmin = 1.25, xmax = 2.3, ymax = 58.5)+ 
  coord_cartesian(clip = 'off')

#tiff("Supplementary_Figure_6C.tiff", height = 16, width = 18, units = 'cm', res = 800)
p1c
dev.off()


###############################################################################################
####################################    Figure S6D    #########################################

###  Figure S3D) GSEA FDR sig - dotplot - FOCUS, FOxTROT and TRANSBIG --->>>

## Plot with combined GSEA results for three cohorts with geneset significant in at least one cohort.

## FOCUS
focus.posneg <- rbind(focus.gsea.pos, focus.gsea.neg)
head(focus.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
focus.posneg$NAME <- focus.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(focus.posneg) # to select column with NAME, NES, FDR, size
FOCUS <- subset(focus.posneg, select = c(1,4,6,8))
head(FOCUS)
FOCUS <- FOCUS %>% mutate(Cohort = "FOCUS") # add column named 'Cohort'
FOCUS.fdrsig <- FOCUS[FOCUS$FDR.q.val < 0.25 ,] # only FDR sig genesets

## FOxTROT
foxtrot.posneg <- rbind(foxtrot.gsea.pos, foxtrot.gsea.neg)
head(foxtrot.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
foxtrot.posneg$NAME <- foxtrot.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(foxtrot.posneg)
FOxTROT <- subset(foxtrot.posneg, select = c(1,4,6,8))
head(FOxTROT)
FOxTROT <- FOxTROT %>% mutate(Cohort = "FOxTROT")
FOxTROT.fdrsig <- FOxTROT[FOxTROT$FDR.q.val < 0.25 ,]

## TRANSBIG
transbig.posneg <- rbind(transbig.gsea.pos, transbig.gsea.neg)
head(transbig.posneg)
## NOTE: "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" has been misspelled in the GSEA result.
## Amend it so it does not impact result downstream
transbig.posneg$NAME <- transbig.posneg$NAME %>% 
  revalue(c("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY" = 
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"))
names(transbig.posneg)
TRANSBIG <- subset(transbig.posneg, select = c(1,4,6,8))
head(TRANSBIG)
TRANSBIG <- TRANSBIG %>% mutate(Cohort = "TRANSBIG")
TRANSBIG.fdrsig <- TRANSBIG[TRANSBIG$FDR.q.val < 0.25 ,]

## Find out common and different pathway between cohorts
## common and different (FOCUS, FOxTROT, and TRANSBIG)
FcFxTr <- intersect(intersect(FOCUS.fdrsig$NAME, FOxTROT.fdrsig$NAME), 
                    TRANSBIG.fdrsig$NAME) # common in all three

FcFx <- intersect(FOCUS.fdrsig$NAME, FOxTROT.fdrsig$NAME) # common in Fc and Fx
## NOTE: FcFx do not have anything in common than FcFxTr

FcTr <- intersect(FOCUS.fdrsig$NAME, TRANSBIG.fdrsig$NAME) # common in Fc and Tr
FcTr <- FcTr[!FcTr %in% FcFxTr] # common in Fc and Tr excluding FcFxTr

FxTr <- intersect(FOxTROT.fdrsig$NAME, TRANSBIG.fdrsig$NAME) # common in Fx and Tr
FxTr <- FxTr[!FxTr %in% FcFxTr] # common in Fx and Tr excluding FcFxTr

## FOCUS only
Fconly <- FOCUS.fdrsig$NAME
Fconly <- Fconly[!Fconly %in% FcFxTr]
Fconly <- Fconly[!Fconly %in% FcTr]

## FOxTROT only
FxOnly <- FOxTROT.fdrsig$NAME
FxOnly <- FxOnly[!FxOnly %in% FcFxTr]
FxOnly <- FxOnly[!FxOnly %in% FxTr]

## TRANSBIG only
TrOnly <- TRANSBIG.fdrsig$NAME
TrOnly <- TrOnly[!TrOnly %in% FcFxTr]
TrOnly <- TrOnly[!TrOnly %in% FcTr]
TrOnly <- TrOnly[!TrOnly %in% FxTr]


## rbind three cohort
all3cohorts <- rbind(FOCUS, FOxTROT, TRANSBIG)
dim(all3cohorts)
head(all3cohorts)


## add new variable column
all3cohort2 <- all3cohorts %>% mutate(Common = as.character(ifelse(all3cohorts$NAME %in% FcFxTr, "FcFxTr",
                                                                   ifelse(all3cohorts$NAME %in% FcTr, "FcTr",
                                                                          ifelse(all3cohorts$NAME %in% FxTr, "FxTr",
                                                                                 ifelse(all3cohorts$NAME %in% Fconly, "Fconly",
                                                                                        ifelse(all3cohorts$NAME %in% FxOnly, "FxOnly",
                                                                                               ifelse(all3cohorts$NAME %in% TrOnly, "TrOnly", NA))))))))


## select fdrsig hallmarks from all cohorts
FDRsig.name <- as.character(all3cohort2[as.numeric(as.character(all3cohort2$FDR.q.val)) < 0.25 ,]$NAME)
FDRsig.name <- FDRsig.name[!duplicated(FDRsig.name)] #remove duplicates

## select these fdr.sig hallmark from all cohorts
all3fdr <- all3cohort2[all3cohort2$NAME %in% FDRsig.name ,]

## remove strings
all3fdr$NAME <- gsub("HALLMARK ", "", gsub("_", " ", all3fdr$NAME))
head(all3fdr)

### ---

### Combined dot plot ---


p1 <- all3fdr[all3fdr$Common == "FcFxTr" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) + # leave blank if ns
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')),# add ns to blanks
            colour = "black") + 
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "FcFxTr" ,]$SIZE))+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 0, face = "bold"),
        legend.position="right",
        legend.key = element_rect(fill = "white")) +  
  ggtitle("Gene Set Enrichement Analysis - Hallmarks") +
  theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, 
                                  margin = ggplot2::margin(0,0,0.5,0,'cm'))) +
  scale_x_discrete(position = "top")

p2 <- all3fdr[all3fdr$Common == "FcTr" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "FcTr" ,]$SIZE))+
  guides(size = guide_legend(title = 'Size\n(FDR Sig)', override.aes = list(linetype = "blank")), 
         color = guide_colourbar(title = 'NES'), shape = guide_legend())+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +
  labs(shape = "FDR.q.val") +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"))


p3 <- all3fdr[all3fdr$Common == "FxTr" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "FxTr" ,]$SIZE))+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        legend.position="right",
        legend.key = element_rect(fill = "white"))

p4 <- all3fdr[all3fdr$Common == "TrOnly" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "TrOnly" ,]$SIZE))+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        legend.position="right",
        legend.key = element_rect(fill = "white"))

p5 <- all3fdr[all3fdr$Common == "Fconly" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "Fconly" ,]$SIZE))+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        legend.position="right",
        legend.key = element_rect(fill = "white"))

p6 <- all3fdr[all3fdr$Common == "FxOnly" ,] %>%
  ggplot(aes(x = Cohort, y = reorder(NAME, NES), color = NES,
             shape = ifelse(FDR.q.val < 0.25, 19, 144))) +
  scale_shape_identity()+
  geom_text(aes(label = ifelse(FDR.q.val < 0.25, ' ', 'ns')), colour = "black")+
  geom_point(stat = "identity", aes(size = all3fdr[all3fdr$Common == "FxOnly" ,]$SIZE))+
  guides(size = FALSE, color = FALSE, shape = FALSE)+
  scale_size(limits = c(25, 200)) +
  xlab("")+
  ylab("")+
  theme_minimal()+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")+
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        limits = c(-2.5,2.5)) +  labs(shape = "FDR.q.val")+
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        legend.position="right",
        legend.key = element_rect(fill = "white"))


# plot separate legend
legend_p2 <- get_legend(p2 + theme(legend.position="right"))

plot1 <- plot_grid(p1,p2+theme(legend.position="none"),p3,p4,p5,p6, 
                   align = "v", nrow = 6, rel_heights = c(0.7, 1.2,0.3, 0.7,0.4, 0.2))
plot2 <- plot_grid(plot1, legend_p2, ncol = 2, rel_widths = c(1, 0.1))

### ---

#tiff("Supplementary_Figure_6D.tiff", height = 25, width = 20, units = 'cm', res = 800)
plot2
dev.off()


###############################################################################################
######################################   Figure S7    #########################################
###############################################################################################

###############################################################################################
###############################    Figure S7A and S7B   #######################################

###  Figure S7A) Fibroblasts MCP score vs DDIR - barplot/scatterplot - FOCUS --->>>
###  Figure S7B) Fibroblasts MCP score vs DDIR - barplot/scatterplot - FOxTROT --->>>

data.1 <- focus.data
#data.1 <- foxtrot.data

variable.1 <-  "Fibroblasts"

data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

## Statistics:
t.test.1 <-  t.test(data.1[[variable.1]] ~ DDIR, data = data.1)
wilcox.1 <- wilcox.test(data.1[[variable.1]] ~ DDIR, data = data.1)
pearson.1 <- cor.test(data.1[["DDIR_scores"]], data.1[[variable.1]], data = data.1,
                      alternative = "two.sided", method = "pearson")

pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")
pvalue.3 <- ifelse(pearson.1$p.value > 0.0001, 
                   paste0("P (two-tailed) = ",round(pearson.1$p.value, 4)),
                   "P (two-tailed) < 0.0001")       


## Set up CMS colour
col.1 <- c("#E18C28","#0B5F9C","#C46294","#129063","#A4A4A4")
names(col.1) <- c("CMS1", "CMS2", "CMS3", "CMS4", "UNK")
col.1 <- col.1[data.1[, 'CMS']]

### ---

#tiff("Supplementary_Figure_7A.tiff", height = 16, width = 18, units = 'cm', res = 800)
#tiff("Supplementary_Figure_7B.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(3.1,3.1,3.1,4.5)+1.5)

## Plot 1 (Scatterplot):
par(fig=c(0,1,0,0.85))
plot(data.1[[variable.1]], data.1[["DDIR_scores"]],
     pch = 19, col = col.1, ylab = "", xlab = "", cex = 1.6, cex.axis = 1.6, 
     ylim = c(-0.3, 0.65))

## Add legend 1 (CMS)
legend('topright', legend = c("CMS1", "CMS2", "CMS3", "CMS4", "UNK"), title = list("CMS"),
       col = c("#E18C28","#0B5F9C","#C46294","#129063","#A4A4A4"),pch = 19, title.adj = 0.55,
       bty = "n", cex = 1.4, xpd = TRUE, inset = c(-0.25, -0.03), adj = c(0.1, 0.5), y.intersp = 0.9)

## Add titles
title(ylab = paste0("DDIR Scores"), line =2.8, cex.lab = 1.6)
title(xlab = paste0("Fibroblasts MCP scores"), line = 2.8, cex.lab = 1.6)

## Add stats 1 (pearson)
mtext(paste0("Pearson r = ", round(pearson.1$estimate, digits = 4)),
      side = 3, adj = 0.02, line = -1.4, cex = 1.4)
mtext(paste0(pvalue.3), side = 3,adj = 0.02, line = -2.6, cex = 1.4)

## Add DDIR threshold (CRC = 0.1094)
abline(h=0.1094, col = "red", lty = 2, lwd = 3)

## Add line of best fit
abline(lm(data.1[["DDIR_scores"]] ~ data.1[[variable.1]]), lwd = 2)

## Plot 2 (boxplot):
par(fig=c(0,1,0.65,1), new=TRUE) # adding another graph to the same plot
boxplot(data.1[[variable.1]] ~ data.1[["DDIR"]],  horizontal=TRUE, axes=FALSE, 
        col = c("blue", "red"))

## Add legend 2 (DDIR)
legend("topright", legend = c("POS", "NEG"), title = list("DDIR"),
       col = c("red", "blue"),pch = 15, title.adj = 0.78, y.intersp = 0.9,
       bty = "n", cex = 1.4, xpd = TRUE,inset = c(-0.22, -1.2), adj = c(0.1,0.5))

## Add titles
title(main =  paste0("Fibroblasts MCP scores\n",
                     ## select cohorts for label
                     cohort.1,
                     #cohort.2,
                     " trial cohort"),
      line = 1.4, cex.main = 1.8)

## Add stats 2 (t.test)
mtext(pvalue.1, side = 1, line = 0.5, cex = 1.4, adj = 0.02)
mtext(pvalue.2, side = 1, line = 1.7, cex = 1.4, adj = 0.02)

## dev off
dev.off()


###############################################################################################
######################################   Figure S8    #########################################
###############################################################################################

###############################################################################################
######################################   Figure S8A   #########################################

###  Figure S8A) IPA table that indicates activated upstream regulators of the 9 genes

###  Go to Supplementary Figure 4C for the table


###############################################################################################
######################################   Figure S8B   #########################################

###  Figure S8B) 9gene score against CMS Subtypes --->>>

data.1 <- focus.data
data.1[["CMS"]] <- plyr::revalue(data.1[["CMS"]], c("Unclassified"= "UNK"))

## Statistics:
kruskal.1 <- kruskal.test(data.1[["sum.score"]] ~ as.factor(data.1[["CMS"]]), data = data.1)
pvalue.1 <- ifelse(kruskal.1$p.value > 0.0001, 
                   paste0("Kruskal-Wallis, p = ", round(kruskal.1$p.value, 4)),
                   paste0("Kruskal-Wallis, p < 0.0001"))

### ---

#tiff("Supplementary_Figure_8B.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar = c(2.1,4.1,4,2.1))

boxplot(data.1[["sum.score"]] ~ data.1[["CMS"]], data = data.1,
        #ylim = c(0, 4500),
        cex.axis = 1.2, outline = FALSE,
        ylab = "", xlab = "")
beeswarm(data.1[["sum.score"]] ~ data.1[["CMS"]], data = data.1,
         add = TRUE, pch = 19, cex = 1.6, spacing = 1.2,
         method = "swarm", corral="wrap",
         col = c("#E18C28", "#0B5F9C", "#C46294", "#129063", "#A4A4A4"))

## Add titles
title(main = paste0("9 gene scores by CMS subtypes\n",cohort.1," trial cohort"), line = 1, cex.main = 1.8)
title(ylab = paste0("9 gene scores"), line = 2.6, cex.lab = 1.6)

## dev off
dev.off()


###############################################################################################
######################################   Figure S8C   #########################################

###  Figure S8C) 9gene score against MSI --->>>

data.1 <- focus.data

## Statistics:
t.test.1 <- t.test(sum.score ~ MSI, data = data.1)
pvalue.1 <- ifelse(t.test.1$p.value > 0.0001, 
                   paste0("t-test, p = ",round(t.test.1$p.value, 4)),
                   "t-test, p < 0.0001")
wilcox.1 <- wilcox.test(sum.score ~ MSI, data = data.1)
pvalue.2 <- ifelse(wilcox.1$p.value > 0.0001, 
                   paste0("Wilcoxon test, p = ",round(wilcox.1$p.value, 4)),
                   "Wilcoxon test, p < 0.0001")

### ---

#tiff("Supplementary_Figure_8C.tiff", height = 16, width = 18, units = 'cm', res = 800)

## Plot:
par(mar=c(4.1,4.1,4.1,2.1)+0.5)

boxplot(sum.score ~ MSI, data = data.1, 
        outline = FALSE, 
        yaxt = 'n',
        cex.axis = 1.6,cex.lab = 1.6)
axis(side = 2, labels = T, cex.axis = 1.3)
beeswarm(sum.score ~ MSI, data = data.1, 
         method = "swarm", corral = "wrap", pch = 19,
         col = c("darkorange", "blue4"),cex = 1.5, spacing = 1.2,
         add = TRUE)

## Add titles
title(ylab = paste0("9 gene scores"), line = 2.8, cex.lab = 1.6)
title(main = paste0("MSI status vs 9 gene scores\n",cohort.1," trial cohort"), cex.main = 1.6)

## Add stats
mtext(pvalue.1, side = 1, line = 2.4, cex = 1.4)
mtext(pvalue.2, side = 1, line = 3.6, cex = 1.4)

## dev off
dev.off()


###############################################################################################
####################################   Session Info   #########################################
###############################################################################################

sessionInfo()

#---------------------------------------------------------------------------------------------#

R version 3.4.0 (2017-04-21)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
  [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] ggpubr_0.2.3    magrittr_1.5    cowplot_0.9.4   plyr_1.8.4      forcats_0.4.0   stringr_1.4.0   dplyr_0.8.3    
[8] purrr_0.3.2     readr_1.3.1     tidyr_0.8.3     tibble_2.1.3    tidyverse_1.2.1 CePa_0.6        beeswarm_0.2.3 
[15] ggplot2_3.2.1  

loaded via a namespace (and not attached):
  [1] tidyselect_0.2.5    haven_2.1.1         lattice_0.20-38     colorspace_1.4-1    generics_0.0.2      vctrs_0.2.0        
[7] stats4_3.4.0        yaml_2.2.0          rlang_0.4.0         pillar_1.4.2        glue_1.3.1          withr_2.1.2        
[13] Rgraphviz_2.22.0    BiocGenerics_0.24.0 modelr_0.1.5        readxl_1.3.1        ggsignif_0.6.0      munsell_0.5.0      
[19] gtable_0.3.0        cellranger_1.1.0    rvest_0.3.4         labeling_0.3        parallel_3.4.0      broom_0.5.2        
[25] Rcpp_1.0.2          backports_1.1.4     scales_1.0.0        graph_1.56.0        jsonlite_1.6        digest_0.6.21      
[31] hms_0.5.1           stringi_1.4.3       cli_1.1.0           tools_3.4.0         lazyeval_0.2.2      crayon_1.3.4       
[37] pkgconfig_2.0.3     zeallot_0.1.0       xml2_1.2.2          lubridate_1.7.4     assertthat_0.2.1    httr_1.4.1         
[43] rstudioapi_0.10     R6_2.4.0            igraph_1.2.4.1      nlme_3.1-131        compiler_3.4.0  

###############################################################################################
###############################################################################################


