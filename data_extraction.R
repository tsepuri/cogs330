# Code implemented based on https://raw.githubusercontent.com/dcdace/Domain-general/master/Analysis_Notebook_v14.ipynb
# Script used to do some data mangling for diagrams and to create CSV files in a format usable by the Python Bayesian modeling file bayes.ipnyb
library(plyr)
library(dplyr)
# data contains the change in activations of rHC and lM1 for think/no-think and go/stop tasks
data <-
  read.csv('https://raw.githubusercontent.com/dcdace/Domain-general/master/data/PSC_targetROIs.csv')
# dataTnt contains the learned word information relating to how many words a participant is able to remember
dataTNT <- read.csv('https://raw.githubusercontent.com/dcdace/Domain-general/master/data/behavioural_TNT.csv')
# ccdata contains the slopes for decrease in inhibitory control/activation detected in rlPFC, rHC and lM1
ccdata <- read.csv('https://raw.githubusercontent.com/dcdace/Domain-general/master/data/decoding_perRun.csv')
# df contains the classification accuracies of the multivoxel pattern analyzer on the rlPFC, rHC and lM1
df <- read.csv('https://raw.githubusercontent.com/dcdace/Domain-general/master/data/decoding.csv', fileEncoding ="UTF-8-BOM")
data[, 4:9] = data[, 4:9]-data[, 4]
# Have a look at what it contains
head(data,1)
cat('rois: ')
cat(unique(data$roi))
cat('\nconditions: ')
cat(unique(data$condition))
dataInh <- subset(data, (condition == 'nt' | condition == 's'))
dataExp <- subset(data, (condition == 't' | condition == 'g'))
dataModality$PSC <- rowMeans(dataModality[, 5:7])
dataDiffModality <- data
# create a distribution of inhibition and expression
# condition as percentage of 'g' and 't'
# Subtract Express from Inhibit
dataModality <- cbind(dataInh[, 1:3], dataInh[, 4:9] - dataExp[, 4:9])
dataModality$condition <- as.factor(dataModality$condition)
# Change conditions to the corresponding Modality
dataModality$condition <- revalue(dataModality$condition, c("nt" = "nt_t", "s" = "s_g"))
names(dataModality)[names(dataModality) == 'condition'] <- 'modality'
# Average PSC within a timespan of 2-6s and substracting the onset value to account for pretrial variability.
# Similar as in Levy&Anderson, 2012 (they did 4-8s)
dataDiffModality$PSC <- rowMeans(dataDiffModality[, 5:7])
dataModality$PSC <- rowMeans(dataModality[, 5:7]) 
# Change the Type levels to better reflect the data content
dataModality$modality <- factor(dataModality$modality, levels = c("nt", "s"), labels = c("No-Think", "Stop"))
# See what dataModality contains now
head(dataModality,1)
levels(dataModality$modality)
# Write change in activation x modality data
write.csv(dataModality, "modality.csv", row.names=FALSE)
# Using rm_2by2_anova.R function from https://raw.githubusercontent.com/dcdace/R_functions/
# Specify columns of interest
columns <- list(
  sID = "sid",     # subject ID
  DV  = "PSC", # dependent variable
  Fc1 = "modality",      # within-subj factor 1
  Fc2 = "roi"   # within-subj factor 2
)
# Define plot label names, tilte, colors and error bar type
param <- list(
  y.label    = "% signal change difference",
  Fc1.label  = "Modality",
  Fc2.label  = "Target region",
  cat.color  = c('red', '#ffc000', 'green', 'blue'), 
  errorbar   = 'se'   
)
# the plot title
param$title <- sprintf('%s x %s interaction', param$Fc2.label,
                       param$Fc1.label)

# Run the ANOVA function but don't display the printouts
anova.results <- (rm_2by2_anova(dataDiffModality, columns, param))
p.Fc1.L1 <-
  ifelse(anova.results$pwc1$p[1] < 0.05, stars.pval(anova.results$pwc1$p[1]), "n.s.")
p.Fc1.L2 <-
  ifelse(anova.results$pwc1$p[2] < 0.05, stars.pval(anova.results$pwc1$p[2]), "n.s.")

# add horizontal line behind the plot
anova.results$plot.anova$layers <- c(geom_hline(aes(yintercept = 0), size = 0.1), # the new layer
                                     anova.results$plot.anova$layers) # original layers
plot.targetrois <- anova.results$plot.anova + 
  # Add post-hoc pairwise comparisons
  # might need to adjust the x and y positions
  annotate(
    "text",
    x = 1.38,
    y = -0.06,
    label = p.Fc1.L1,
    color = param$cat.color[1],
    size = 8
  ) +
  annotate(
    "text",
    x = 1.6,
    y = -0.15,
    label = p.Fc1.L2,
    color = param$cat.color[2],
    size = 8
  ) +
  # Change ROI captions
  scale_x_discrete(limits = c("rHC", "lM1"), label = c('Hippocampus', 'M1')) 
options(repr.plot.width = 5, repr.plot.height = 4, repr.plot.res = 200) # change plot size
plot.targetrois

df_filtered <- subset(df, type == "ntnt" | type == "ss")
df_filtered_new <- df_filtered[df_filtered$roi %in% c("rdlpfc", "rvlpfc"),]
# merge the rows where roi is rdlpfc and rvlpfc
df_merged <- df_filtered_new %>%
  group_by(sid, type) %>%
  summarize(acc = mean(acc),
            roi = ifelse(all(roi %in% c("rdlpfc", "rvlpfc")), "rlpfc", unique(acc))) %>%
  bind_rows(df_filtered[!df_filtered$roi %in% c("rvlpfc", "rdlpfc"),])
write.csv(df_merged, "accs.csv", row.names=FALSE)
df_filtered_new <- ccdata[ccdata$ROI %in% c("rDLPFC", "rVLPFC"),]
dataTNT[, 10:ncol(dataTNT)] <- dataTNT[, 10:ncol(dataTNT)] * 100
options(width = 114)
print(cbind(dataTNT[, 4:10], round(dataTNT[, 11:ncol(dataTNT)], 0)))
# merge the rows where roi is rdlpfc and rvlpfc
df_merged <- df_filtered_new %>%
  group_by(sNR) %>%
  summarize(SIF = mean(SIF),
            SSRT = mean(SSRT),
            propslope = mean(propslope),
            ROI = ifelse(all(ROI %in% c("rDLPFC", "rVLPFC")), "rlpfc", unique(SIF, SSRT))) %>%
  bind_rows(ccdata[!ccdata$ROI %in% c("rVLPFC", "rDLPFC"),])
# Using sid to mark participant number
df_merged$sid = df_merged$sNR
df_filtered_merged <- distinct(df_merged, sid, .keep_all = TRUE)
write.csv(df_merged, "slopes.csv", row.names=FALSE)
