rm(list = ls())
library(plyr)
library(survival)
library(gdata)
library(survival)
library(reshape2)
library(ggplot2); library(reshape2)
library(limma)
library(ggsignif)

## working directory assumed to be inside a GitHub folder on Mac Desktop
setwd("~/Desktop/TN-05_GitHub/data_files/")

## set defaults
theme_set(theme_bw(26) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("point", aes(size = 3))

###### set variables

cut = c(50)
visit.time = c(26)
auc.time = c(26)
NV = paste("NV", auc.time, sep = "")
top = c(0.75)
bot = c(0.75)

## load counts
counts <- read.delim("tmm_norm_Ritux_counts_outliers_removed_collapsed.txt")

#### load annotation
anno <- read.delim("TN Unblinded Specimen List_Linsley 20140924_trimmed.txt")
anno = na.omit(anno)
colnames(anno)[10] = c("id")
anno = anno[c(10, 1,2,12, 13,15,16)]

anno.sub = subset(anno, anno$Visit == visit.time & !anno$Treatment == "Placebo") # for active only
#anno.sub = subset(anno, anno$Visit == visit.time) # for all subjects
anno.sub = anno.sub[,-(5:7)]

## load demographcs
demo <- read.delim("TN05 Data_Linsley 20141126_trimmed.txt")

## merge anno and demographcs
anno = merge(anno, demo, by = "id")

## load and reformat C-peptide data
mmtt <- read.delim("TN05 MMTT Data_Linsley 20141020_full_AUC_redo.txt")
mmtt = mmtt[,-(16:ncol(mmtt))]
colnames(mmtt)[1] = c("id")
m = mmtt[,4:ncol(mmtt)]
m[is.na(m)] = 0
mmtt.2 = data.frame(mmtt[,1:3], m)
mmtt.2$Treatment = as.character(mmtt.2$Treatment)
mmtt.2$Responder. = as.character(mmtt.2$Responder.)
melt.mmtt = melt(mmtt.2, id=c("id","Treatment", "Responder."))
melt.mmtt = na.omit(melt.mmtt)
colnames(melt.mmtt)[5] = c("AUC")
melt.mmtt = subset(melt.mmtt, melt.mmtt$variable %in% NV)

anno.sub = subset(anno, libID %in% colnames(counts))
anno.auc = merge(anno.sub, melt.mmtt, by = "id")
anno.auc = subset(anno.auc, !Visit == "M30")

anno.auc = anno.auc[(c("id", "Treatment.x", "Responder..x", "age", "Gender", "Race"))]
colnames(anno.auc) = c("id",  "treatment", "responder", "age", "gender", "race")
anno.auc = anno.auc[!(duplicated(anno.auc$id)),]

## load flow data
flow <- read.delim("TN05 Flow Data_Linsley_20150112_trimmed.txt")
flow = flow[,-(14:23)]
colnames(flow)[4] = c("id")
flow.use = subset(flow, id %in% anno.auc$id)

flow.anno = merge(flow.use, anno.auc, by = "id")
flow.anno = subset(flow.anno, Result_type == "RPTD")

## count samples at each visit

temp = subset(flow.anno, Result_name %in% c("CD19_PCT", "CD3_PCT", "CD4_PCT", "CD8_PCT"))
temp$Result_name = as.character(temp$Result_name)
table(temp$VisitNo, temp$treatment, temp$Result_name)

## load curvefit survival data

surv.cf <- read.csv("TN05.study_days_to_threshold_percent_AUC_from_models.csv", stringsAsFactors=FALSE)
surv.cf[,1]= NULL
colnames(surv.cf) = c("id", "separate", "combined")
surv.cf$id = gsub("TN05_", "", surv.cf$id)
surv.cf = data.frame(id = surv.cf$id, survival = surv.cf$combined) # use combined model
surv.cf$survival = gsub("Inf", 7567, surv.cf$survival)
surv.cf[is.na(surv.cf)] <- 7567
surv.cf$survival = as.numeric(surv.cf$survival)

surv.cf$survival = ifelse(surv.cf$survival>365*2, 365*2, surv.cf$survival)
surv.cf$status = ifelse(surv.cf$survival== (365*2), 0, 1)
surv.cf$status = ifelse(surv.cf$survival >= (365*2), 0, 
                        ifelse(surv.cf$survival < (365*2), 1, 1e6))
#
## calculate p-values for active versus placebo by test by visit
xnames = as.character(unique(flow.anno$VisitNo))
ncx = length(xnames)

ynames = unique(flow.use$Result_name)
ncy = length(ynames)

DF <- data.frame(visit = as.numeric(), flow.name = as.character(), p.val = as.numeric(), stringsAsFactors = F)

for (i in 1:ncx) { 
  q = xnames[i]
  qsub = subset(flow.anno, VisitNo == q)
  
  for (j in 1:ncy) { 
    
    r = as.character(ynames[j])
    rsub = subset(qsub, Result_name == r)
    rsub = subset(rsub, !Result == c("."))
    rsub$Result = as.character(rsub$Result)
    rsub$Result = as.numeric(rsub$Result)
    
    qnt.top <- subset(flow.anno, treatment == "Active")
    qnt.bot <- subset(flow.anno, treatment == "Placebo")
    
    top.set = subset(rsub, id %in% qnt.top$id)
    bot.set = subset(rsub, id %in% qnt.bot$id)
    
    test = wilcox.test(top.set$Result, bot.set$Result)
    pval = test$p.value
    
    result = as.character(c(q, r, pval))
    result = t(data.frame(result, stringsAsFactors = F))
    colnames(result) = c("visit", "flow.name", "p.val")			
    DF = rbind(result, DF, stringsAsFactors = F)
    colnames(DF) = c("visit", "flow.name", "p.val")		
  }
}
DF$p.val = as.numeric(DF$p.val)
DF.s = DF[order(DF$p.val),]
p.bh = as.data.frame(p.adjust(DF.s$p.val, method = c("BH")))
DF.s.corr = data.frame(DF.s, p.bh)
colnames(DF.s.corr) = c("visit", "flow.name", "p.value", "adj.p.value")
DF.s.corr$visit = as.numeric(DF.s.corr$visit)

## sanity check
temp = subset(flow.anno, VisitNo == 26 & Result_name == "CD3_PCT") # gave 53 subjects

## normalize flow values
xnames = unique(flow.use$Result_name)
ynames = as.character(unique(flow.anno$VisitNo))
ncx = length(xnames)
ncy = length(ynames)

DF <- data.frame(id = as.numeric(), Result_name = as.character(), VisitNo = as.numeric(), p.val = as.numeric(), stringsAsFactors = F)

for (i in 1:ncx) { 
  q = as.character(xnames[i])
  qsub = subset(flow.anno, Result_name == q)
  qsub = subset(qsub, !Result == c("."))
  qsub$Result = as.character(qsub$Result)
  qsub$Result = as.numeric(qsub$Result)
  qsub$zResult = (qsub$Result-mean(qsub$Result))/sd(qsub$Result)
  
  for (j in 1:ncy) { 
    
    r = as.character(ynames[j])
    rsub = subset(qsub, VisitNo == r)
    
    qnt.top <- subset(flow.anno, treatment == "Active")
    qnt.bot <- subset(flow.anno, treatment == "Placebo")
    
    top.set = subset(rsub, id %in% qnt.top$id)
    bot.set = subset(rsub, id %in% qnt.bot$id)
    
    test = wilcox.test(top.set$zResult, bot.set$zResult)
    pval = test$p.value
    
    result = data.frame(id = rsub$id, Result_name = q, VisitNo = r, p.value = pval, treatment = rsub$treatment, zResult =  rsub$zResult)
    
    DF = rbind(result, DF, stringsAsFactors = F)
    colnames(DF) = c("id", "Result_name", "VisitNo", "p.value", "treatment", "zResult")		
    
  }
}
DF$p.val = as.numeric(DF$p.val)
DF.s = DF[order(DF$p.val),]
DF.s$p.bh = p.adjust(DF.s$p.val, method = c("BH"))
colnames(DF.s)[ncol(DF.s)] = c("adj.p.value")
DF.s$VisitNo = as.character(DF.s$VisitNo)
DF.s$VisitNo = as.numeric(DF.s$VisitNo)

## Fig S4
temp = subset(flow.anno, VisitNo == 26 & Result_name == "CD3_PCT") # gave 53 subjects

## plot z normalized results

to.plot = subset(DF, Result_name %in% c("CD19_PCT", "CD3_PCT", "CD4_PCT", "CD8_PCT"))
to.plot$Result = as.character(to.plot$Result)
to.plot$Result = as.numeric(to.plot$Result)
to.plot$treatment = gsub("Active", "R", to.plot$treatment)
to.plot$treatment = gsub("Placebo", "P", to.plot$treatment)
to.plot$VisitNo = as.character(to.plot$VisitNo)
to.plot$VisitNo = as.numeric(to.plot$VisitNo)

dev.off()
quartz(height =10, width =14, dpi =72);
ggplot(to.plot, aes(x= treatment, y = zResult)) + geom_boxplot(fill = "#8da0cb") + facet_grid(Result_name~VisitNo)
last_plot() + geom_dotplot(data = to.plot, aes(x= treatment, y = zResult), position = "identity", method = "histodot", binaxis = "y", stackdir = "center", dotsize = 0.5) 
last_plot() + scale_y_continuous(limits = c(-4, 8))
last_plot() + geom_signif(comparisons = list(c("R", "P")), y_position = 5, textsize = 6,  map_signif_level = T)
last_plot() + labs(y = "% marker positive cells, normalized", x = "Visit, weeks")

## Fig 2B

####
tgc <- ddply(DF.s, c("VisitNo", "Result_name", "treatment"), summarise,
             N    = length(zResult),
             mean = mean(zResult),
             sd   = sd(zResult),
             se   = sd / sqrt(N)
)
tgc = na.omit(tgc)

to.plot = subset(tgc, Result_name %in% c("CD19_PCT", "CD3_PCT", "CD4_PCT", "CD8_PCT"))
#to.plot$treatment = gsub("Active", "A", to.plot$treatment)
#to.plot$treatment = gsub("Placebo", "P", to.plot$treatment)
to.plot$VisitNo = as.character(to.plot$VisitNo)
to.plot$VisitNo = as.numeric(to.plot$VisitNo)

dev.off()
quartz(height =10, width =14, dpi =72);
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("errorbar", aes(size = 2))

theme_set(theme_bw(30) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
cbPalette = c("#ca0020", "#0571b0") # 

ggplot(to.plot, aes(x= VisitNo, y = mean, group = treatment)) + geom_line(aes(colour = treatment)) + facet_wrap(~Result_name)
last_plot() + scale_colour_manual(values=cbPalette)
last_plot() + geom_errorbar(data = to.plot, aes(ymin=mean-se, ymax=mean+se, width=2, colour = treatment))
last_plot() + labs(y = "Mean % marker positive cells, normalized", x = "Visit, weeks")





