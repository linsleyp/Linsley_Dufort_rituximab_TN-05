### scripts to analyze C-peptide, RNAseq, and other data for the TN-05 rituximab trial

##### set up environment #####

rm(list=ls())

## load general packages
library(xlsx)
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1)))
library(ggthemes)

## load analysis-specific packages
library(limma)
library(edgeR)
library(lmerTest)
library(annotables)

# load useful functions (packages installable from github.com/mjdufort)
library(miscHelpers)
library(RNAseQC)
library(countSubsetNorm)
library(limmaTools)
library(geneSetTools)

## working directory assumed to be inside a GitHub folder on Mac Desktop
setwd("~/Desktop/TN-05_GitHub/data_files/")

##### set up color palette #####

color.placebo <- "#0571b0"
color.ritux <- "#ca0020"


##### import detailed study schedule, with visit numbers, names, etc. #####

study_schedule <-
  read.csv(file="TN05_trial_schedule.csv", header=TRUE) %>%
  standardize_dimnames()
for (i in
     colnames(study_schedule)[
       !(colnames(study_schedule) %in% c("study", "notes")) &
       !str_detect(colnames(study_schedule), "name|scheduled")]) {
  if (is.character(study_schedule[[i]]))
    study_schedule[[i]] <- study_schedule[[i]] %>%
      str_replace_all("<0", "NA") %>%
      as.numeric()
}
glimpse(study_schedule)


##### import CBC data, and transform values #####

## import cleaned CBC data
cbc_data.untransformed <-
  read.csv(file="TN05_CBC_data_cleaned.csv", header=TRUE)

## create version to transform
cbc_data <- cbc_data.untransformed

## transform percentages and raw numbers from CBCs
raw.tmp <-
  c(paste0(
    c("basophils", "eosinophils", "lymphocytes", "monocytes", "neutrophils"), "_abs"),
    "rbc", "wbc")
percent.tmp <-
  c("basophils", "eosinophils", "lymphocytes", "monocytes", "neutrophils")

# log-transform raw values
for (i in raw.tmp) {
  if (any(na.omit(cbc_data.untransformed[,i])==0)) {
    cbc_data[,i] <- log(cbc_data.untransformed[,i]+0.01)
  } else {
    cbc_data[,i] <- log(cbc_data.untransformed[,i])
  }
}

# arcsin-transform percent values
for (i in percent.tmp) {
  cbc_data[,i] <- asin(sqrt(cbc_data.untransformed[,i]/100))
}

## sort the data by patient and visit number
cbc_data <-
  cbc_data %>%
  arrange(participant_id, cbc_visit_number)

rm_tmp(ask=FALSE)


##### import patient data (age, sex, visit dates, etc.) #####

patient_data <-
  read.csv(file="TN05_patient_data.csv", header=TRUE)


##### import C-peptide data #####

## import cleaned C-peptide data
cpeptide_data <-
  read.csv(file="TN05_Cpeptide_data_cleaned.csv", header=TRUE)

### determine visits with c-peptide previously at detection threshold, to exclude them
## this is necessary to exclude the "bottomed-out" points from the curve fitting
# otherwise they would incorrectly bias the lowest slopes upward

# lower limit of detection in nmol/L/min
ng_ml_to_nmol_l_min <- (1000 / 3020.29) * (1 / 120)
low_lim_detect <- 3 * ng_ml_to_nmol_l_min

cpeptide_data$prev_at_cpep_detect_thresh <- as.logical(NA)
cpeptide_data <-
  cpeptide_data %>%
  arrange(participant_id, cpeptide_study_day)

for (i in unique(cpeptide_data$participant_id)) { # iterate over patient
  rows.tmp <-
    which(
      with(cpeptide_data,
           participant_id==i & !is.na(auc2hr))) # determine all rows for current patient
  cpeptide_data$prev_at_cpep_detect_thresh[rows.tmp[1]] <- FALSE # set first visit to FALSE
  if (length(rows.tmp) > 1)
    for (j in (2:length(rows.tmp))) {
      cpeptide_data$prev_at_cpep_detect_thresh[rows.tmp[j]] <-
        (cpeptide_data$auc2hr[rows.tmp[j-1]] < low_lim_detect * 1.01) |  # previous visit is at detection limit
        cpeptide_data$prev_at_cpep_detect_thresh[rows.tmp[j-1]]  # previous visit had prev_at_cpep_detect_thresh
    }
}

rm_tmp(ask=FALSE)


##### fit curves to c-peptide AUC by cpeptide_study_day #####

## scale cpeptide_study_day so that it works better for model fitting
# this is more important for higher-order models, which are not used here
days.sd <- sd(cpeptide_data$cpeptide_study_day, na.rm=TRUE)
cpeptide_data$cpeptide_study_day_scaled <-
  scale(cpeptide_data$cpeptide_study_day, center=FALSE,
        scale=days.sd) %>%
  as.vector()


## fit mixed-effects models to log_auc2hr and log_prop_auc2hr
# with terms to allow different fixed-effect slopes by treatment

# log_auc2hr with random intercept by patient
lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random <-
  lmer(
    log_auc2hr ~
      (1|participant_id) +
      cpeptide_study_day_scaled +
      cpeptide_study_day_scaled:treatment +
      (cpeptide_study_day_scaled-1|participant_id),
    data=cpeptide_data[!cpeptide_data$prev_at_cpep_detect_thresh,])

# log_prop_auc2hr with random intercept by patient
lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random <-
  lmer(
    log_prop_auc2hr ~
      1 + (1|participant_id) +
      cpeptide_study_day_scaled +
      cpeptide_study_day_scaled:treatment +
      (cpeptide_study_day_scaled-1|participant_id),
    data=cpeptide_data[!cpeptide_data$prev_at_cpep_detect_thresh,])

# log_prop_auc2hr with random intercept by patient (used for plotting Fig. S5D)
lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random <-
  lmer(
    log_prop_auc2hr ~
      0 +
      # cpeptide_study_day_scaled +
      cpeptide_study_day_scaled:treatment +
      (cpeptide_study_day_scaled-1|participant_id),
    data=cpeptide_data[!cpeptide_data$prev_at_cpep_detect_thresh,])


##### plot curves fit to log_auc2hr by cpeptide_study_day #####

x_min.tmp <- -60  # set lower boundary for cpeptide_study_day
x_max.tmp <- 1025 # set upper boundary for cpeptide_study_day

## random effect linear with intercept, with fixed effect for separate slopes for placebo and rituximab
pdf("TN05.log_auc2hr_vs_cpeptide_study_day.fitted_curve.linear_random_with_intercept.color_by_treatment.pdf",
    w=8, h=6)
par(mar=c(5,6.5,4,2)+0.1)
plot(x=0, type="n", xlim=c(x_min.tmp, x_max.tmp)/days.sd, ylim=c(-4.8,0.9),
     xlab="Years in study", ylab="C-peptide 2hr AUC\n(nmol / L / min)", xaxt="n", yaxt="n",
     cex.lab=1.5)
axis(1, at=seq(0,3,by=0.5)*365.25/days.sd,
     labels=seq(0,3,by=0.5),
     cex.axis=1.5)
# axis(1, at = seq(0,1000, by=250)/days.sd,
#      labels = seq(0,1000, by=250),
#      cex.axis=1.5)
axis(2, at = log(10^(-2:0)), labels=10^(-2:0),
     cex.axis=1.5, las=1)
for (i in unique(cpeptide_data$participant_id)) {
  if (sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),]))==3) {
    b0.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'(Intercept)'[
        match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))]
    b1.tmp <-
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] +
      coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
        match(i, rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] *
      (cpeptide_data$treatment[match(i, cpeptide_data$participant_id)]=="rituximab")
    plot(
      function(x) {b0.tmp + (x*b1.tmp)},
      from=x_min.tmp/days.sd,
      to=max(cpeptide_data$cpeptide_study_day[cpeptide_data$participant_id==i], na.rm=TRUE)/days.sd,
      col=
        setNames(c(color.placebo,color.ritux), c("placebo", "rituximab"))[
          unique(cpeptide_data$treatment[cpeptide_data$participant_id ==i])],
      lwd=1.5,
      add=TRUE)
  }
}
dev.off()

rm_tmp(ask=FALSE)


##### extract slopes from curves fit to log_prop_auc2hr and log_auc2hr by cpeptide_study_day #####

## extract slope for each patient
cpeptide_data$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept <-
  as.numeric(NA)
cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept <-
  as.numeric(NA)
cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random <-
  as.numeric(NA)

for (i in 1:nrow(cpeptide_data)) {
  cat("Starting row", i, "of", nrow(cpeptide_data), "\n")
  
  ## extract slopes for log_prop_auc2hr_cpeptide_study_day
  # with separate fixed effects for placebo and rituximab
  ## units are log units per year
  
  # with random intercept by patient
  if (sum(!is.na(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(
      cpeptide_data$participant_id[i],
      rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),])) == 3)
    cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[i] <-
      (coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(cpeptide_data$participant_id[i],
              rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] +
         coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
           match(cpeptide_data$participant_id[i],
                 rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] *
         (cpeptide_data$treatment[i]=="rituximab")) /
      days.sd * 365.25

  # without fitting intercept (used for plotting only)
  if (sum(!is.na(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id[
    match(
      cpeptide_data$participant_id[i],
      rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id)),])) == 3) {
    
    treatment.tmp <- 
      cpeptide_data$treatment[i]=='rituximab'
    cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random[i] <-
      (coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(
          cpeptide_data$participant_id[i],
          rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))] +
      treatment.tmp *
      coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
        match(
          cpeptide_data$participant_id[i],
          rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))] +
      (!treatment.tmp) *
      coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentplacebo'[
        match(
          cpeptide_data$participant_id[i],
          rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))]) /
      days.sd * 365.25
  }
  
    
  ## extract slope for log_auc2hr_cpeptide_study_day
  # with separate fixed effects for placebo and rituximab
  if (sum(!is.na(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id[
    match(
      cpeptide_data$participant_id[i],
      rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id)),])) == 3)
    cpeptide_data$log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept[i] <-
      (coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(cpeptide_data$participant_id[i],
              rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] +
         coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
           match(cpeptide_data$participant_id[i],
                 rownames(coef(lmer.log_auc2hr_linear_cpeptide_study_day_with_intercept.patient_random)$participant_id))] *
         (cpeptide_data$treatment[i]=="rituximab")) /
      days.sd * 365.25
  
}

rm_tmp(ask=FALSE)


##### calculate time to threshold AUC from models #####

## if the value of the predictor is < 0, that means the AUC is trending up, and there is no timepoint
# set those to Inf

## units are days

cpeptide_data[
  , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random_with_intercept"] <-
  log(0.5) /
  cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept *
  365.25

# without a fitted intercept (use only for plotting)
cpeptide_data[
  , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random"] <-
  log(0.5) /
  cpeptide_data$log_prop_auc2hr_slope_cpeptide_study_day_linear_random *
  365.25


# set values <= 0 to infinity (because, by the models, these patients will never hit the threshold)
cpeptide_data[
  !is.na(cpeptide_data[
    , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random_with_intercept"]) &
    cpeptide_data[
      , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random_with_intercept"] <= 0,
  "study_days_to_0.5prop_baseline_auc2hr.model_linear_random_with_intercept"] <-
  Inf

cpeptide_data[
  !is.na(cpeptide_data[
    , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random"]) &
    cpeptide_data[
      , "study_days_to_0.5prop_baseline_auc2hr.model_linear_random"] <= 0,
  "study_days_to_0.5prop_baseline_auc2hr.model_linear_random"] <-
  Inf

## export file with study_days_to_threshold, for use in downstream analyses
write.csv(
  unique(
    cpeptide_data[
      ,c("participant_id",
         "study_days_to_0.5prop_baseline_auc2hr.model_linear_random_with_intercept",
         "study_days_to_0.5prop_baseline_auc2hr.model_linear_random")]),
  file="TN05.study_days_to_threshold_prop_AUC_from_models.csv")

rm_tmp(ask=FALSE)


##### plot time to threshold C-peptide, colored by treatment #####

x_min.tmp <- 0  # set lower boundary for cpeptide_study_day
x_max.tmp <- 1025 # set upper boundary for cpeptide_study_day

# log_prop_auc2hr, random effect linear with intercept
pdf("TN05.log_prop_auc2hr_vs_cpeptide_study_day.fitted_curve.linear_random.colored_by_treatment.pdf",
    w=8.3, h=6, useDingbats=FALSE)
par(mar=c(5,7.5,4,2)+0.1)
plot(x=0, type="n",
     xlim=c(x_min.tmp, x_max.tmp)/days.sd, ylim=c(-1.4, 0.7),
     xlab="Years in study", ylab="",
     xaxt="n", yaxt="n", cex.lab=1.5)
axis(1, at=seq(0,3,by=0.5)*365.25/days.sd,
     labels=seq(0,3,by=0.5),
     cex.axis=1.5)
axis(2, at = log(c(0.25, 0.5, 1, 2)),
     labels = paste0(c(25, 50, 100, 200), "%"),
     # labels = seq(0, 150, by=50),
     las=1, cex.axis=1.5)
mtext("C-peptide 2hr AUC\n(percent of baseline)", side=2, line=4.5, cex=1.5) # add y-axis
abline(h = log(0.5), lty=2, lwd=1.5)
# for (i in "TN05_757136") { # test one individual
for (i in unique(cpeptide_data$participant_id)) {
  if (sum(!is.na(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id[
    match(i, rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id)),]))==3) {
    treatment.tmp <- filter(cpeptide_data, participant_id==i)[1,"treatment"]=='rituximab'
    col.tmp <- if (is.na(treatment.tmp)) "grey50" else if (treatment.tmp) color.ritux else color.placebo
    b0.tmp <- 0
    b1.tmp <-
      coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$cpeptide_study_day_scaled[
        match(i, rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))] +
      treatment.tmp *
      coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentrituximab'[
        match(i, rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))] +
      (!treatment.tmp) *
      coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id$'cpeptide_study_day_scaled:treatmentplacebo'[
        match(i, rownames(coef(lmer.log_prop_auc2hr_linear_cpeptide_study_day.patient_random)$participant_id))]
    # print(paste(i, b0.tmp, b1.tmp))
    plot(function(x) {b0.tmp + (x*b1.tmp)},
         from=0/days.sd,
         to=(x_max.tmp*1.25/days.sd),
         col=col.tmp, add=TRUE)
    points( # add points for time to 50% of baseline?
      x=unique(
        cpeptide_data$study_days_to_0.5prop_baseline_auc2hr.model_linear_random[
          cpeptide_data$participant_id == i]) / days.sd,
      y=log(0.5),
      col=col.tmp, pch=16, cex=1, add=TRUE)
  }
}
segments(
  x0=min(cpeptide_data$study_days_to_0.5prop_baseline_auc2hr.model_linear_random, na.rm=T)/days.sd,
  y0=log(0.5), y1=log(0.1),
  col=color.placebo, lty=3, lwd=3)
mtext(
  paste0(
    round(
      min(cpeptide_data$study_days_to_0.5prop_baseline_auc2hr.model_linear_random, na.rm=T),0),
    "\ndays"),
  at=min(cpeptide_data$study_days_to_0.5prop_baseline_auc2hr.model_linear_random, na.rm=T)/days.sd,
  side=1, line=1.1, cex=1, col=color.placebo) # add days to threshold for one patient
dev.off()

rm_tmp(ask=FALSE)


##### plot basic C-peptide variables #####

## plot auc2hr by individual vs cpeptide_study_day
pdf("TN05.auc2hr_vs_cpeptide_study_day.by_treatment.pdf", w=9, h=6)
ggplot(
  data=cpeptide_data[!is.na(cpeptide_data$auc2hr),],
  mapping=aes(
    x=cpeptide_study_day, y=auc2hr, group=participant_id,
    color=treatment)) +
  geom_line(size=0.3) +
  scale_color_manual(
    "Treatment\ngroup",
    values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
  labs(x="Days in study", y="C-peptide 2hr AUC\n(nmol/L/min)") +
  geom_hline(yintercept=low_lim_detect, linetype="dashed")
dev.off()

## plot log AUC by individual vs cpeptide_study_day
pdf(
  "TN05.log_auc2hr_vs_cpeptide_study_day.diagnosis_window.colored_by_treatment.no_legend.pdf",
  w=8, h=5.3)
ggplot(
  data=cpeptide_data[!is.na(cpeptide_data$log_auc2hr),],
  mapping=aes(
    x=cpeptide_study_day, y=log_auc2hr, group=participant_id,
    color=treatment)) +
  geom_polygon(
    inherit.aes=FALSE,
    data=data.frame(
      cpeptide_study_day=
        rep(c(-90,-10), each=2),
      log_auc2hr=
        c(range(cpeptide_data$log_auc2hr, na.rm=TRUE) + c(-1,1),
          rev(range(cpeptide_data$log_auc2hr, na.rm=TRUE) + c(-1,1)))),
    mapping=aes(x=cpeptide_study_day, y=log_auc2hr),
    fill="gray80") +
  geom_line(size=0.5) +
  scale_color_manual(
    # "Treatment\ngroup",
    values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
  guides(color=FALSE) +
  labs(x="Years in study", y="C-peptide 2hr AUC\n(nmol / L / min)") +
  scale_x_continuous(
    breaks=seq(0,3,by=0.5)*365.25,
    labels=seq(0,3,by=0.5)) +
  scale_y_continuous(breaks=log(10^(-2:0)), labels=10^(-2:0)) +
  coord_cartesian(
    xlim=as.numeric(range(cpeptide_data$cpeptide_study_day, na.rm=TRUE) + c(30,0)),
    ylim=range(cpeptide_data$log_auc2hr, na.rm=TRUE)) +
  geom_hline(yintercept=log(low_lim_detect), linetype="dashed")
dev.off()


##### plot C-peptide levels over time, with individuals highlighted to illustrate variation #####

## plot AUC by individual vs cpeptide_study_day, highlighting individuals, with detection limit and diagnosis window
# identify individuals to highlight differences
cpeptide_data[
  with(cpeptide_data,
       age_years >= 12 &
         age_years <= 14 &
         cpeptide_study_day %between% c(-30, 30)),
  c("participant_id",
    "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept",
    "auc2hr", "log_auc2hr")] %>%
  plyr::rename(
    replace=c(
      "log_auc2hr_slope_cpeptide_study_day_linear_random_with_intercept"=
        "cpeptide_change")) %>%
  arrange(cpeptide_change)

## pair to use:
# TN05_757136, TN05_924617

# plot log AUC, with points and fit lines for those two patients only
pdf("TN05.log_auc2hr_vs_cpeptide_study_day.colored_by_treatment.det_limit.diagnosis_window.highlight_individuals.points_and_smooths.no_legend.pdf",
    w=8, h=5.3,
    useDingbats=FALSE)
ggplot(
  data=cpeptide_data,
  mapping=aes(x=cpeptide_study_day, y=log_auc2hr, group=participant_id, color=treatment)) +
  geom_polygon(
    inherit.aes=FALSE,
    data=data.frame(
      cpeptide_study_day=
        rep(c(-90,-10), each=2),
      log_auc2hr=
        c(range(cpeptide_data$log_auc2hr, na.rm=TRUE) + c(-1,1),
          rev(range(cpeptide_data$log_auc2hr, na.rm=TRUE) + c(-1,1)))),
    mapping=aes(x=cpeptide_study_day, y=log_auc2hr),
    fill="gray80") +
  scale_color_manual("Treatment\ngroup", values=c("placebo"=color.placebo, "rituximab"=color.ritux)) +
  guides(color=FALSE) +
  geom_line(size=0.3, alpha=0.4) +
  geom_point(
    data=cpeptide_data[cpeptide_data$participant_id=="TN05_757136",],
    color=color.placebo, size=4) + # plot one fast progressor
  geom_smooth(
    data=cpeptide_data[
      with(cpeptide_data, participant_id=="TN05_757136" & !prev_at_cpep_detect_thresh),],
    se=FALSE, method="lm", color=color.placebo) +
  geom_point(
    data=cpeptide_data[cpeptide_data$participant_id=="TN05_924617",],
    color=color.ritux, size=4) + # plot one slow progressor
  geom_smooth(
    data=cpeptide_data[
      with(cpeptide_data, participant_id=="TN05_924617" & !prev_at_cpep_detect_thresh),],
    se=FALSE, method="lm", color=color.ritux) +
  scale_x_continuous(
    breaks=seq(0,3,by=0.5)*365.25,
    labels=seq(0,3,by=0.5)) +
  scale_y_continuous(breaks=log(10^(-2:0)), labels=10^(-2:0)) +
  labs(y="C-peptide 2hr AUC\n(nmol / L / min)", x="Years in study") +
  geom_hline(yintercept=log(low_lim_detect), linetype="dashed") +
  coord_cartesian(
    xlim=as.numeric(range(cpeptide_data$cpeptide_study_day, na.rm=TRUE) + c(30,0)),
    ylim=range(cpeptide_data$log_auc2hr, na.rm=TRUE))
dev.off()


##### import RNAseq data, and merge with other data #####

# read in counts (already QC'd and TMM-normalized)
counts <- # read in counts
  read.table(
    "TN05_RNAseq_counts_normalized_QCd.txt",
    sep="\t", header=T)

## read in sample annotation
rnaseq_annotation <- # read in annotation
  read.csv("TN05_RNAseq_annotation.csv", header=T)
rnaseq_annotation$treatment <-
  factor(
    rnaseq_annotation$treatment,
    levels=c("placebo", "rituximab"))

# merge patient_data into rnaseq_annotation data
rnaseq_annotation <-
  merge(rnaseq_annotation,
        patient_data[
          c("participant_id",
            setdiff(colnames(patient_data), colnames(rnaseq_annotation)))],
        by="participant_id", all.x=T)

# generate RNAseq and C-peptide visit IDs for matching across data types
rnaseq_annotation$rnaseq_visit_id <-
  with(rnaseq_annotation,
       paste(participant_id, rnaseq_visit_number, sep="_"))
rnaseq_annotation$cpeptide_visit_id <-
  with(rnaseq_annotation,
       paste(participant_id, cpeptide_visit_number, sep="_"))

# merge cpeptide_data with rnaseq_annotation
rnaseq_annotation <-
  merge(rnaseq_annotation,
        cpeptide_data[
          c("cpeptide_visit_id",
            setdiff(colnames(cpeptide_data), colnames(rnaseq_annotation)))],
        by="cpeptide_visit_id", all.x=T)


##### filter RNAseq annotation and counts to include only shared libraries, and put them in the same order #####

setdiff(rnaseq_annotation$libid, colnames(counts))
setdiff(colnames(counts), rnaseq_annotation$libid)
# 8 libraries are present in annotation, but were cut from counts for QC reasons

rnaseq_annotation.merged <-
  rnaseq_annotation[rnaseq_annotation$libid %in% colnames(counts),] %>%
  arrange(participant_id, rnaseq_visit_number)
counts.merged <- counts[,match(rnaseq_annotation.merged$libid, colnames(counts))]
glimpse(rnaseq_annotation.merged)
glimpse(counts.merged)


##### plot saturation curve of RNAseq data #####

saturation.merged <-
  estimate_saturation(counts.merged, verbose=TRUE, method="division", ndepths=20)
pdf("saturation_curve_all_libs_by_division.pdf", w=9, h=6)
plot_saturation_curve(saturation.merged, plot_points=FALSE, plot_lines=TRUE)
dev.off()
# looks good; problematic libraries already removed


##### generate simple combined objects with all libraries passing QC cuts #####

## these objects have had low-quality libraries removed
## the data were normalized previously
counts.final <- counts.merged
master.final <- rnaseq_annotation.merged


##### merge CBC data into master.final #####

master.final <-
  base::merge(
    master.final,
    cbc_data[,setdiff(colnames(cbc_data), colnames(master.final))],
    by.x="rnaseq_visit_id",
    by.y="cbc_visit_id",
    all.x=TRUE) %>%
  arrange(rnaseq_visit_id)


##### log-transform the normalized counts, and calculate median gene set expression for some gene modules #####

## generate vwts.all object
vwts.all <-
  calc_norm_counts(
    counts=counts.final, design=master.final, libID_col="libid",
    min_cpm=1, min_libs_perc=0.15,
    normalize=FALSE,
    return_DGEcounts=TRUE) %>%
  voom(plot=TRUE)

## read in Speake et al. gene sets
gene_sets.Speake_Linsley <-
  read_delim(
    file="Speake_Linsley_2014_gene_modules.gmx",
    delim="\t") %>%
  as.list() %>% # convert to a list of character vectors of genes
  lapply(FUN=na.omit)

## calculate median expression in log2(cpm+1) of CD19.mod
master.final$median_CD19.mod <-
  gene_set_median_count(
    gene_sets.Speake_Linsley[["CD19.mod"]], vwts.all,
    remove_low_count_genes=TRUE)[
      match(master.final$libid, colnames(vwts.all$E))]


##### fit limma models for treatment effect (rituximab- vs. placebo-treated) at each time point #####

# generate lists to store results by visit
for (i in c("master.treatment.by_visit", "DesignMat.treatment.by_visit", "vwts.treatment.by_visit",
            "vfit.treatment.by_visit", "topGenes.treatment.by_visit", "vwts_sig.treatment.by_visit"))
  assign(i, list())

for (i in sort(unique(master.final$rnaseq_visit_weeks))) {
  j <- as.character(i)
  condition.tmp <-
    with(master.final,
         !is.na(treatment) &
           (rnaseq_visit_weeks %in% i))
  master.treatment.by_visit[[j]] <-
    master.final %>%
    filter(condition.tmp) %>%
    fix_factors()
  DGECounts.treatment.by_visit.tmp <-
    calc_norm_counts(
      counts=counts.final,
      design=master.treatment.by_visit[[j]],
      libID_col="libid",
      min_cpm=1, min_libs_perc=0.15,
      normalize=FALSE,
      return_DGEcounts=TRUE,
      group=master.treatment.by_visit[[j]]$participant_id)
  master.treatment.by_visit[[j]] <-
    master.treatment.by_visit[[j]][
      match(colnames(DGECounts.treatment.by_visit.tmp),
            master.treatment.by_visit[[j]][,"libid"]),]
  
  DesignMat.treatment.by_visit[[j]] <-
    model.matrix(
      ~ treatment + sex,
      data=master.treatment.by_visit[[j]])
  vwts.treatment.by_visit[[j]] <-
    voomWithQualityWeights(
      DGECounts.treatment.by_visit.tmp,
      design=DesignMat.treatment.by_visit[[j]], plot=TRUE)
  vfit.treatment.by_visit[[j]] <-
    lmFit(vwts.treatment.by_visit[[j]]) %>%
    eBayes()
  topGenes.treatment.by_visit[[j]] <-
    topTable(
      vfit.treatment.by_visit[[j]],
      coef = 2,
      number=Inf, sort.by="P")
}

rm_tmp(ask=FALSE)


##### export limma model results #####

## combine all results into a data frame
topGenes.treatment.by_visit.combined <-
  data.frame(
    gene=sort(unique(unlist(sapply(topGenes.treatment.by_visit, rownames)))))
topGenes.treatment.by_visit.combined <-
  topGenes.treatment.by_visit.combined[
    !str_detect(topGenes.treatment.by_visit.combined$gene, "^[0-9]"),,drop=FALSE]
for (i in names(topGenes.treatment.by_visit)) {
  topGenes.treatment.by_visit.combined[
    , c(paste("week", i, "logFC", sep="_"), paste("week", i, "adj.P.Val", sep="_"))] <-
    topGenes.treatment.by_visit[[i]][
      match(topGenes.treatment.by_visit.combined$gene, rownames(topGenes.treatment.by_visit[[i]])),
      c("logFC", "adj.P.Val")]
  }

topGenes.treatment.by_visit.combined <-
  topGenes.treatment.by_visit.combined %>%
  dplyr::arrange(
    week_26_adj.P.Val, week_52_adj.P.Val, week_78_adj.P.Val, week_104_adj.P.Val, week_130_adj.P.Val, week_0_adj.P.Val)

write.table(
  topGenes.treatment.by_visit.combined,
  file="ritux_vs_placebo_DE_genes_by_visit.txt",
  row.names=FALSE, quote=FALSE, sep="\t")


##### fit limma models for treatment effect (rituximab- vs. placebo-treated) at each time point, adjusted for age #####

for (i in c("master.treatment.age_years.by_visit", "DesignMat.treatment.age_years.by_visit",
            "vwts.treatment.age_years.by_visit", "vfit.treatment.age_years.by_visit",
            "topGenes.treatment.age_years.by_visit", "vwts_sig.treatment.age_years.by_visit"))
  assign(i, list())

for (i in sort(unique(master.final$rnaseq_visit_weeks))) {
  j <- as.character(i)
  condition.tmp <-
    with(master.final,
         !is.na(treatment) &
           !is.na(age_years) &
           (rnaseq_visit_weeks %in% i))
  master.treatment.age_years.by_visit[[j]] <-
    master.final %>%
    filter(condition.tmp) %>%
    fix_factors()
  DGECounts.treatment.age_years.by_visit.tmp <-
    calc_norm_counts(
      counts=counts.final,
      design=master.treatment.age_years.by_visit[[j]],
      libID_col="libid",
      min_cpm=1, min_libs_perc=0.15,
      normalize=FALSE,
      return_DGEcounts=TRUE,
      group=master.treatment.age_years.by_visit[[j]]$participant_id)
  master.treatment.age_years.by_visit[[j]] <-
    master.treatment.age_years.by_visit[[j]][
      match(colnames(DGECounts.treatment.age_years.by_visit.tmp),
            master.treatment.age_years.by_visit[[j]][,"libid"]),]
  
  DesignMat.treatment.age_years.by_visit[[j]] <-
    model.matrix(
      ~ treatment + age_years + sex,
      data=master.treatment.age_years.by_visit[[j]])
  vwts.treatment.age_years.by_visit[[j]] <-
    voomWithQualityWeights(
      DGECounts.treatment.age_years.by_visit.tmp,
      design=DesignMat.treatment.age_years.by_visit[[j]], plot=TRUE)
  vfit.treatment.age_years.by_visit[[j]] <-
    lmFit(vwts.treatment.age_years.by_visit[[j]]) %>%
    eBayes()
  topGenes.treatment.age_years.by_visit[[j]] <-
    topTable(
      vfit.treatment.age_years.by_visit[[j]],
      coef = 2,
      number=Inf, sort.by="P")
}
# this excludes the 130 week visit, because there aren't enough samples


##### fit limma models for treatment effect (rituximab- vs. placebo-treated) at each time point, adjusted for CBC numbers #####

for (i in c("master.treatment.cbc_percents.by_visit", "DesignMat.treatment.cbc_percents.by_visit",
            "vwts.treatment.cbc_percents.by_visit", "vfit.treatment.cbc_percents.by_visit",
            "topGenes.treatment.cbc_percents.by_visit", "vwts_sig.treatment.cbc_percents.by_visit"))
  assign(i, list())

for (i in sort(unique(master.final$rnaseq_visit_weeks))) {
  j <- as.character(i)
  condition.tmp <-
    with(master.final,
         !is.na(treatment) &
           !is.na(basophils) &
           !is.na(eosinophils) &
           !is.na(lymphocytes) &
           !is.na(monocytes) &
           !is.na(neutrophils) &
           (rnaseq_visit_weeks %in% i))
  master.treatment.cbc_percents.by_visit[[j]] <-
    master.final %>%
    filter(condition.tmp) %>%
    fix_factors()
  DGECounts.treatment.cbc_percents.by_visit.tmp <-
    calc_norm_counts(
      counts=counts.final,
      design=master.treatment.cbc_percents.by_visit[[j]],
      libID_col="libid",
      min_cpm=1, min_libs_perc=0.15,
      normalize=FALSE,
      return_DGEcounts=TRUE,
      group=master.treatment.cbc_percents.by_visit[[j]]$participant_id)
  master.treatment.cbc_percents.by_visit[[j]] <-
    master.treatment.cbc_percents.by_visit[[j]][
      match(colnames(DGECounts.treatment.cbc_percents.by_visit.tmp),
            master.treatment.cbc_percents.by_visit[[j]][,"libid"]),]
  
  DesignMat.treatment.cbc_percents.by_visit[[j]] <-
    model.matrix(
      ~ treatment +
        lymphocytes + neutrophils + monocytes + eosinophils + basophils +
        sex,
      data=master.treatment.cbc_percents.by_visit[[j]])
  vwts.treatment.cbc_percents.by_visit[[j]] <-
    voomWithQualityWeights(
      DGECounts.treatment.cbc_percents.by_visit.tmp,
      design=DesignMat.treatment.cbc_percents.by_visit[[j]], plot=TRUE)
  vfit.treatment.cbc_percents.by_visit[[j]] <-
    lmFit(vwts.treatment.cbc_percents.by_visit[[j]]) %>%
    eBayes()
  topGenes.treatment.cbc_percents.by_visit[[j]] <-
    topTable(
      vfit.treatment.cbc_percents.by_visit[[j]],
      coef = 2,
      number=Inf, sort.by="P")
}
# this one excludes the 130 week visit, because there aren't enough samples

rm_tmp(ask=FALSE)


##### calculate median expression of "persistently down-regulated" genes #####

# this set uses the genes that are significantly DE at 78 weeks (FDR < 0.01)
gene_set.TN05_persistent_78weeks_standard <-
  rownames(topGenes.treatment.by_visit$`78`)[
    with(topGenes.treatment.by_visit$`78`,
         adj.P.Val < 0.01 &
           logFC < -log2(1.5))]
master.final$median_gene_set.TN05_persistent_78weeks_standard <-
  gene_set_median_count(
    gene_set.TN05_persistent_78weeks_standard,
    vwts.all[,match(master.final$libid, colnames(vwts.all$E))],
    min_median_gene_expression=0.25)


##### plot my persistently downregulated gene set vs CD19.mod over time, averaged across patients within treatment groups #####


## from my gene_set.TN05_persistent_78weeks_standard gene set, using only genes significantly DE at 78 weeks in the vanilla model
master.tmp <-
  master.final %>%
  tidyr::gather(
    key=gene_set,
    value=median_count,
    median_CD19.mod, median_gene_set.TN05_persistent_78weeks_standard) %>%
  filter(gene_set %in% c("median_CD19.mod", "median_gene_set.TN05_persistent_78weeks_standard") &
           rnaseq_visit_weeks < 130) %>%
  group_by(treatment, rnaseq_visit_weeks, gene_set) %>%
  summarise(mean_median_count=mean(median_count), se_median_count=se(median_count))
master.tmp$gene_set <-
  factor(master.tmp$gene_set,
         levels=c("median_gene_set.TN05_persistent_78weeks_standard",
                  "median_CD19.mod"))

pdf("gene_set_median_CD19.mod_vs_TN05_persistent_78weeks_standard.vs_weeks.by_treatment.pdf",
    w=8,h=8)
ggplot(
  master.tmp,
  mapping=aes(
    x=rnaseq_visit_weeks,
    y=mean_median_count,
    group=treatment,
    color=treatment)) +
  geom_line(size=1.5) +
  scale_color_manual("Treatment", values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
  scale_x_continuous(breaks=seq(0,104,length.out=5)) +
  geom_errorbar(
    mapping=aes(
      ymin = mean_median_count - se_median_count,
      ymax = mean_median_count + se_median_count),
    color="black", size=2, width=0) +
  facet_wrap(~gene_set, nrow=2, scales="free_y") +
  labs(x="Visit Week", y="Median gene set expression\nlog (cpm + 0.5), mean +/- SE")
dev.off()

rm_tmp(ask=FALSE)


##### plot means and SE of individual genes from these gene sets #####

genes.tmp <- c("NETO1", "MS4A1")

master.tmp <- master.final
for (gene.tmp in genes.tmp)
  master.tmp[,gene.tmp] <-
  vwts.all$E[
    gene.tmp, match(master.tmp$libid, colnames(vwts.all$E))]

master.tmp <-
  master.tmp %>%
  tidyr::gather_(
    key_col="gene",
    value_col="counts",
    gather_cols=genes.tmp) %>%
  filter(rnaseq_visit_weeks < 130) %>%
  group_by(treatment, rnaseq_visit_weeks, gene) %>%
  summarise(mean_counts=mean(counts), se_counts=se(counts))
master.tmp$gene <-
  factor(master.tmp$gene,
         levels=genes.tmp)

pdf(
  paste0(
    genes.tmp[1], "_vs_", genes.tmp[2], ".vs_weeks.by_treatment.pdf"),
  w=8,h=8)
ggplot(
  master.tmp,
  mapping=aes(
    x=rnaseq_visit_weeks,
    y=mean_counts,
    group=treatment,
    color=treatment)) +
  geom_line(size=1.5) +
  scale_color_manual("Treatment", values=c("rituximab"=color.ritux, "placebo"=color.placebo)) +
  scale_x_continuous(breaks=seq(0,104,length.out=5)) +
  geom_errorbar(
    mapping=aes(
      ymin = mean_counts - se_counts,
      ymax = mean_counts + se_counts),
    color="black", size=2, width=0) +
  facet_wrap(~gene, nrow=2, scales="free_y") +
  labs(x="Visit Week", y="Gene expression\nlog (cpm + 0.5), mean +/- SE")
dev.off()

rm_tmp(ask=FALSE)


##### make volcano plot with genes labelled, at 78 week visit #####

plot_volcano_2var(
  topGenes.treatment.by_visit[["78"]],
  plotdims=c(7,9),
  file_prefix="volcano.cbc_percents.by_visit.week_78.with_gene_labels",
  my_cols=c("gray40", color.placebo),
  fc_cut=log2(1.5), p_cut=0.01,
  gene_labs="threshold", y_cut=5, x_cut= -log2(1.5),
  x_cut_direction="lower", gene_labs_repel=FALSE, gene_lab_size=5)


