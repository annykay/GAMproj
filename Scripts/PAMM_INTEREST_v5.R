library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(zoo)
library(survival)
library(mgcv)
library(purrr)
library(pammtools)
Set1 <- RColorBrewer::brewer.pal(9, "Set1")



#######  to transform NSCLC dataset

data <- read.csv("SourceData/megred_table_JM.csv")
head(data)

data1 <- data[data$STUDY == "INTEREST",]
head(data1)

data1_tte <- data1[!duplicated(data1$USUBJID),]
head(data1_tte)
fit1 <- survfit(Surv(time=EVENT_TIME, event=EVENT) ~ 1, data=data1_tte)
print(fit1)
plot(x=c(0,850), y=c(0,1), type = 'n', xlab = "Time,days", ylab = "PFS" )
lines(fit1, col="red", lwd = 2)
######################

##### Make dataset with TDCs

data1s <- data1[, c("USUBJID", "TIME","STUDY","LINE")]
data1hr <- data1s[!duplicated(data1s),]


for (i in 1:7){
  data1si <- data1[data1$YTYPE == i, c("USUBJID", "TIME","DV","YTYPE_DESC")]
  names(data1si)[3] <- data1si$YTYPE_DESC[1]
  data1hr <- merge(data1hr, data1si[,1:3], by = c("USUBJID", "TIME") , all = T)
}
# head(data1hr)
subjects <- unique(data1hr$USUBJID)
data1hr <- data1hr[
  with(data1hr, order(USUBJID, TIME)),
]
for (i in c(1:(nrow(data1hr) - 1))) {
  if (is.na(data1hr$SLD[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$SLD[i+1] <- data1hr$SLD[i]
  }
  if (is.na(data1hr$AST[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$AST[i+1] <- data1hr$AST[i]
  }
  if (is.na(data1hr$ALP[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$ALP[i+1] <- data1hr$ALP[i]
  }
  if (is.na(data1hr$ALT[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$ALT[i+1] <- data1hr$ALT[i]
  }
  if (is.na(data1hr$CREAT[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$CREAT[i+1] <- data1hr$CREAT[i]
  }
  # if (is.na(data1hr$NEUT[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
  #   data1hr$NEUT[i+1] <- data1hr$NEUT[i]
  # }
  if (is.na(data1hr$WBC[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$WBC[i+1] <- data1hr$WBC[i]
  }
}
BM_by_id <-  data1hr %>% group_by(USUBJID) %>%
  summarise( SLDav = mean(SLD, na.rm = T), ALPav = mean(ALP, na.rm = T), 
             CREATav = mean(CREAT, na.rm = T), ALPmax =max(ALP, na.rm = T), 
             CREATmax = max(CREAT, na.rm = T), ASTav = mean(AST, na.rm = T),
             ALTav = mean(ALT, na.rm = T))

#### combine TTE and TDC data
data1_tte <-  data1[!duplicated(data1$USUBJID),c(1:5,9:11,18)]
data1_tte$SLD <- 0
data1_tte$WBC <- 0
#data1_tte$NEUT <- 0
data1_tte$CREAT <- 0
data1_tte$ALP <- 0
data1_tte$ALT <- 0
data1_tte$AST <- 0
head(data1_tte)

ped_data1 <- as_ped(
  data    = list(data1_tte, data1hr),
  cut = seq(0,1000, by =7),
  formula = Surv(time=EVENT_TIME, event=EVENT) ~ . + concurrent(SLD, CREAT,ALP, ALT, AST, WBC, tz_var = "TIME"), 
  id      = "USUBJID")

head(ped_data1)

summary(ped_data1$ALP)
summary(ped_data1$ALT)
summary(ped_data1$AST)
summary(ped_data1$CREAT)
summary(ped_data1$SLD)
summary(ped_data1$WBC)
##############


####### "best" model - biomarker effects  on individual survival functions prediction distributions

ped_data1upd <- merge(ped_data1,BM_by_id, by= "USUBJID", all = T)

stat <- ped_data1upd %>% group_by(USUBJID) %>% summarize(SLDna = sum(is.na(SLD)), ASTna = sum(is.na(AST)), 
                                                         ALPna = sum(is.na(ALP)), ALTna = sum(is.na(ALT)), 
                                                         CREATna = sum(is.na(CREAT)), count = n())
completely_missing <- unique(stat$USUBJID[stat$CREATna == stat$count])
completely_missing <- c(completely_missing, unique(stat$USUBJID[stat$ALPna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$ALTna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$ASTna == stat$count]))
completely_missing <- unique(completely_missing)
ped_data1upd <- ped_data1upd[!(ped_data1upd$USUBJID) %in% completely_missing,  ]

for (i in c(1:nrow(ped_data1upd))) {
  if (is.na(ped_data1upd$AST[i])) {
    ped_data1upd$AST[i] <- ped_data1upd$ASTav[i]
  }
  if (is.na(ped_data1upd$ALP[i])) {
    ped_data1upd$ALP[i] <- ped_data1upd$ALPav[i]
  }
  if (is.na(ped_data1upd$ALT[i])) {
    ped_data1upd$ALT[i] <- ped_data1upd$ALTav[i]
  }
  if (is.na(ped_data1upd$CREAT[i])) {
    ped_data1upd$CREAT[i] <- ped_data1upd$CREATav[i]
  }
}
pam_strata1 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + ALP + SLD + ALT + s(CREAT),
                   data = ped_data1upd,
                   family = poisson(), offset = offset, 
                   method = 'REML')

summary(pam_strata1)
AIC(pam_strata1)

tidy_smooth(pam_strata1) %>%
  ggplot(aes(x = x, y = fit)) +
  geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3) +
  geom_step() + geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~ylab) +
  xlab(expression(t)) + ylab(expression(f[p](t) %.% x[p]))


SurvFunc_strat_plot <- function(ped_data,  pam_model, strata_name,
                                 strata_type = "cat", timecat = 500, prob_range= c(0.05, 0.95),
                                legend_title, legend_names){
  ped_data_IDs <- unique(ped_data$USUBJID)
  print(length(ped_data_IDs))
  if (strata_type == "cont"){
    BMmedian <- median(ped_data[,strata_name])
    strata_name_cat <- paste0(strata_name,"_cat")
    ped_data[,strata_name_cat] <- 0
    ped_data[,strata_name_cat] <- ifelse(ped_data[,strata_name] < BMmedian, 1, 2)
  }
  else{
    strata_name_cat <- strata_name
  }
  
  surf_sum_df <- data.frame()
  for (i in c(1:length(ped_data_IDs))){
     ped_data_i <- ped_data[ped_data$USUBJID==ped_data_IDs[i],]
     if (nrow(ped_data_i) > 0) {
         ped_data_i$intlen <- ped_data_i$tend - ped_data_i$tstart
         ped_df <- ped_data_i  %>%  add_surv_prob(pam_model) 
         surf_sum_df <- rbind(surf_sum_df, ped_df[,c("USUBJID","tend","surv_prob", strata_name_cat )] )
     }
  }
  surf_sum_df <- surf_sum_df[surf_sum_df$tend <= timecat ,]
  names(surf_sum_df)[4] <- "STRATA"
  
  surf_PI <- surf_sum_df %>% group_by(tend, STRATA) %>%
      summarise(surv_median = median(surv_prob, na.rm = T), 
                surv_lower = quantile(surv_prob, probs=prob_range[1], na.rm = T),
                surv_upper = quantile(surv_prob, probs=prob_range[2], na.rm = T))
  
  p1 <- ggplot(data=surf_PI, aes(x = tend, y = surv_median, ymax = surv_lower, ymin = surv_upper, group = as.factor(STRATA))) + 
    geom_line(aes(col=as.factor(STRATA)),size=2) + scale_color_discrete(name = legend_title, labels = legend_names) + 
    geom_ribbon(aes(fill=as.factor(STRATA)),alpha = 0.2)  + scale_fill_discrete(guide="none") +
    ylab(expression(hat(S)(t))) + xlab(paste0(expression(t),", days")) 
    
  return(p1)
}


survF_plot <- SurvFunc_strat_plot(ped_data=ped_data1upd, pam_model=pam_strata1,strata_name="SMKSTAT",
                                  prob_range= c(0.10, 0.90), legend_title = 'Smoking Status',
                                  legend_names= c('Never Smoker', 'Former Smoker', 'Habital Smoker'))

ggsave(survF_plot, file =  "Results/survf_predDistr_SMKSTAT.png", width = 5, height = 4)

survF_plot <- SurvFunc_strat_plot(ped_data=ped_data1upd, pam_model=pam_strata1,strata_name="WHOSTATN",
                                  timecat =500, prob_range= c(0.10, 0.90), 
                                  legend_title = 'Performance Status',legend_names= c(0,1,2))
ggsave(survF_plot, file =  "Results/survf_predDistr_WHOSTATN.png", width = 5, height = 4)

survF_plot <- SurvFunc_strat_plot(ped_data=ped_data1upd, pam_model=pam_strata1,strata_name="SLDav",strata_type = "cont",
                                   timecat =500, prob_range= c(0.10, 0.90), 
                                  legend_title = 'Tumor Size',legend_names= c('< population average', '> population average' ))
ggsave(survF_plot, file =  "Results/survf_predDistr_SLDav.png", width = 5, height = 4)


survF_plot <- SurvFunc_strat_plot(ped_data=ped_data1upd, pam_model=pam_strata1, strata_name="ALPav",strata_type = "cont",
                                   timecat =500, prob_range= c(0.10, 0.90), legend_title = 'Average ALP level',legend_names= c('< population average', '> population average' ))
ggsave(survF_plot, file =  "Results/survf_predDistr_ALPav.png", width = 5, height = 4)

survF_plot <- SurvFunc_strat_plot(ped_data=ped_data1upd, pam_model=pam_strata1, strata_name="CREATmax",strata_type = "cont",
                                   timecat =500, prob_range= c(0.10, 0.90), legend_title = 'Max CREAT level',legend_names= c('< population average', '> population average' ))
ggsave(survF_plot, file =  "Results/survf_predDistr_CREATmax.png", width =5, height = 4)
##############################
n_train   <- 300
train_idx <- sample(data1_tte$USUBJID, n_train)
# train_ped <- ped_data1upd[ped_data1upd$USUBJID %in% train_idx, ]
# test_ped <- ped_data1upd[!(ped_data1upd$USUBJID %in% train_idx), ]
train_ped <- as_ped(
  data    = list(data1_tte[data1_tte$USUBJID %in% train_idx, ], data1hr[data1hr$USUBJID %in% train_idx, ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID",
  cut = seq(0, 1000, 200))
test_ped <- as_ped(
  data    = list(data1_tte[!(data1_tte$USUBJID %in% train_idx), ], data1hr[!(data1hr$USUBJID %in% train_idx), ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID",
  cut = seq(0, 1000, 200))
# pam1 <- gam(
#   formula = ped_status~s(tend) + WHOSTATN + SMKSTAT + ALP + ALT + SLD + s(CREAT),
#   data = train_ped)
# pam2 <- pamm(
#   formula = ped_status ~ s(tend) + SMKSTAT + WHOSTATN,
#   data = train_ped)
pam1 <- gam(ped_status ~ s(tend) + SMKSTAT + WHOSTATN,
                   data = train_ped,
                   family = poisson(), offset = offset, 
                   method = 'REML')
pam2 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + ALP + ALT + SLD + s(CREAT),
                   data = train_ped,
                   family = poisson(), offset = offset, 
                   method = 'REML')
test_ped$intlen <- test_ped$tend - test_ped$tstart
# time_points <- test_ped$tend[test_ped$ped_status == 1]
time_point <- 200
test_ped <- test_ped %>% add_hazard(pam1, overwrite = T)
test_ped_timepoint <- test_ped[test_ped$tend < time_point,]
for_auc <-  test_ped_timepoint
for_auc$predicts <- for_auc$hazard
for_auc$status <- for_auc$ped_status

# for_auc <- test_ped_timepoint %>% 
#   group_by(USUBJID) %>% 
#   summarize(status = sum(ped_status), 
#             predicts = max(hazard))
roc_curve <- roc(for_auc$status, for_auc$predicts, auc.polygon=TRUE)


test_ped <- test_ped %>% add_hazard(pam2, overwrite = T)
test_ped_timepoint <-  test_ped[test_ped$tend < time_point,]
for_auc <-  test_ped_timepoint
for_auc$predicts <- for_auc$hazard
for_auc$status <- for_auc$ped_status

# for_auc <- test_ped_timepoint %>%
#   group_by(USUBJID) %>%
#   summarize(status = sum(ped_status),
#             predicts = max(hazard, na.rm = T))
roc_curve1 <- roc(for_auc$status, for_auc$predicts, auc.polygon=TRUE)
plot(roc_curve)
plot(roc_curve1)
for_plotting <- data.frame(Model = 'Base Model',
                           Sensitivity = roc_curve[["sensitivities"]], 
                           Specificity = roc_curve[["specificities"]])
for_plotting1 <- data.frame(Model = 'Best Model',
                           Sensitivity = roc_curve1[["sensitivities"]], 
                           Specificity = roc_curve1[["specificities"]])
for_plotting <- rbind(for_plotting, for_plotting1)
for_plotting1 <- data.frame(Model = 'Random Classifier',
                            Sensitivity = c(1, 0), 
                            Specificity = c(0, 1))
for_plotting <- rbind(for_plotting, for_plotting1)
ggplot(aes(x=Specificity, y = Sensitivity, color = Model), data = for_plotting) + geom_line()
ggsave('Results/ROC_all.png')
merged <- merge(data1_tte[, c(1:9)], data1hr, by = c('USUBJID', 'LINE', 'STUDY'), all =  F)
merged <- merge(merged, BM_by_id, by = 'USUBJID')

merged <- merged[!(merged$USUBJID) %in% completely_missing,  ]
for (i in c(1:nrow(merged))) {
  if (is.na(merged$AST[i])) {
    merged$AST[i] <- merged$ASTav[i]
  }
  if (is.na(merged$ALP[i])) {
    merged$ALP[i] <- merged$ALPav[i]
  }
  if (is.na(merged$ALT[i])) {
    merged$ALT[i] <- merged$ALTav[i]
  }
  if (is.na(merged$CREAT[i])) {
    merged$CREAT[i] <- merged$CREATav[i]
  }
}
colnames(merged)[1] <- 'id'
pec <- pec(
  list(pam1 = pam1, pam2 = pam2),
  Surv(EVENT_TIME, EVENT) ~ 1, # formula for IPCW
  data = merged[(merged$id %in% train_idx), ],
  times = seq(.01, 400, by = 10),
  start = .01,
  exact = FALSE,
  reference = T,
  cens.model = 'aalen'
)
plot(pec)
ggsave('Results/PEC.png')


unique_measurments <- ped_data1upd %>% group_by(USUBJID) %>% summarize(SLDuni = length(unique(SLD)), 
                                              ALTuni = length(unique(ALT)),
                                              ASTuni = length(unique(AST)),
                                              ALPuni = length(unique(ALP)),
                                              CREATuni = length(unique(CREAT)),
                                              tenduni = length(unique(tend)), 
                                              SMKSTAT = unique(SMKSTAT), 
                                              WHOSTATN = unique(WHOSTATN))
patients_sld <- unique_measurments$USUBJID[unique_measurments$SLDuni > 6]
patients_sld <- c("S7U52DXE")
ped_data1upd$WHOSTATN <- as.numeric(ped_data1upd$WHOSTATN)
ped_data1upd$SMKSTAT <- as.numeric(ped_data1upd$SMKSTAT)
ped_data1upd <- ped_data1upd %>% add_hazard(pam2, overwrite = T)
for (SLD_SUBj in patients_sld) {
plotting_ped <- ped_data1upd[ped_data1upd == SLD_SUBj, ]
melted <- melt(plotting_ped[, c('USUBJID', 'tend', 'SLD', 'hazard')], id = c('USUBJID', 'tend'))
coeff <- max(plotting_ped$hazard, na.rm = T) / max(plotting_ped$SLD, na.rm = T)


Hazardcolor <- "#F8766D"
SLDcolor <- '#00BFC4'

ggplot(plotting_ped, aes(x=tend)) +
  
  geom_line( aes(y=hazard / coeff, color=SLDcolor)) + 
  geom_line( aes(y=SLD, color=Hazardcolor), size = 1.5) + # Divide by 10 to get the same range than the temperature
  geom_ribbon(aes(ymin=ci_lower/ coeff, ymax=ci_upper/ coeff, fill= SLDcolor), alpha = 0.15) + 
  scale_y_continuous(
    
    # Features of the first axis
    name = "SLD, sm",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Hazard")) + 
      
    theme_minimal() +
    xlab('Time') +
    theme(
        axis.title.y = element_text(color = SLDcolor, size=16),
        axis.title.y.right = element_text(color = Hazardcolor, size=16),
        title = element_text(size = 18)
    ) +
    theme(legend.position = "none")
filename <- paste0('Results/SLD_Hazard_pictures/SLD_', SLD_SUBj, '.png')
ggsave(filename)
}


patients_alt <- unique_measurments$USUBJID[unique_measurments$ALTuni > 10]
for (ALT_SUBj in patients_alt) {
  plotting_ped <- ped_data1upd[ped_data1upd == ALT_SUBj, ]
  melted <- melt(plotting_ped[, c('USUBJID', 'tend', 'ALT', 'hazard')], id = c('USUBJID', 'tend'))
  coeff <- max(plotting_ped$hazard, na.rm = T) / max(plotting_ped$ALT, na.rm = T)
  
  
  Hazardcolor <- "#F8766D"
  SLDcolor <- '#00BFC4'
  
  ggplot(plotting_ped, aes(x=tend)) +
    
    geom_line( aes(y=hazard / coeff, color=SLDcolor )) + 
    geom_line( aes(y=ALT, color=Hazardcolor), size = 1.5) + # Divide by 10 to get the same range than the temperature
    geom_ribbon(aes(ymin=ci_lower/ coeff, ymax=ci_upper/ coeff, fill= SLDcolor), alpha = 0.15) + 
    scale_y_continuous(
      
      # Features of the first axis
      name = "ALT, U/ml",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Hazard")) + 
    
    theme_minimal() +
    xlab('Time') +
    theme(
      axis.title.y = element_text(color = SLDcolor, size=16),
      axis.title.y.right = element_text(color = Hazardcolor , size=16),
      title = element_text(size = 18)
    ) +
    theme(legend.position = "none")
  filename <- paste0('Results/ALT_Hazard_pictures/ALT_', ALT_SUBj, '.png')
  ggsave(filename)
}

patients_alp <- unique_measurments$USUBJID[unique_measurments$ALPuni > 10]
for (ALP_SUBj in patients_alp) {
  plotting_ped <- ped_data1upd[ped_data1upd == ALP_SUBj, ]
  melted <- melt(plotting_ped[, c('USUBJID', 'tend', 'ALP', 'hazard')], id = c('USUBJID', 'tend'))
  coeff <- max(plotting_ped$hazard, na.rm = T) / max(plotting_ped$ALP, na.rm = T)
  
  
  Hazardcolor <- "#F8766D"
  SLDcolor <- '#00BFC4'
  
  ggplot(plotting_ped, aes(x=tend)) +
    
    geom_line( aes(y=hazard / coeff, color=SLDcolor)) + 
    geom_line( aes(y=ALP, color=Hazardcolor), size = 1.5) + # Divide by 10 to get the same range than the temperature
    geom_ribbon(aes(ymin=ci_lower/ coeff, ymax=ci_upper/ coeff, fill= SLDcolor), alpha = 0.15) + 
    scale_y_continuous(
      
      # Features of the first axis
      name = "ALP, U/ml",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Hazard")) + 
    
    theme_minimal() +
    xlab('Time') +
    theme(
      axis.title.y = element_text(color = SLDcolor, size=16),
      axis.title.y.right = element_text(color = Hazardcolor, size=16),
      title = element_text(size = 18)
    ) +
    theme(legend.position = "none")
  filename <- paste0('Results/ALP_Hazard_pictures/ALP_', ALP_SUBj, '.png')
  ggsave(filename)
}
patients_creat <- unique_measurments$USUBJID[unique_measurments$CREATuni > 10]
for (CREAT_SUBj in patients_creat) {
  plotting_ped <- ped_data1upd[ped_data1upd == CREAT_SUBj, ]
  melted <- melt(plotting_ped[, c('USUBJID', 'tend', 'CREAT', 'hazard')], id = c('USUBJID', 'tend'))
  coeff <- max(plotting_ped$hazard, na.rm = T) / max(plotting_ped$CREAT, na.rm = T)
  
  
  Hazardcolor <- "#F8766D"
  SLDcolor <- '#00BFC4'
  
  ggplot(plotting_ped, aes(x=tend)) +
    
    geom_line( aes(y=hazard / coeff, color=SLDcolor)) + 
    geom_line( aes(y=CREAT, color=Hazardcolor ), size = 1.5) + # Divide by 10 to get the same range than the temperature
    geom_ribbon(aes(ymin=ci_lower/ coeff, ymax=ci_upper/ coeff, fill= SLDcolor), alpha = 0.15) + 
    scale_y_continuous(
      
      # Features of the first axis
      name = "CREAT, U/ml",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Hazard")) + 
    
    theme_minimal() +
    xlab('Time') +
    theme(
      axis.title.y = element_text(color = SLDcolor, size=16),
      axis.title.y.right = element_text(color = Hazardcolor , size=16),
      title = element_text(size = 18)
    ) +
    theme(legend.position = "none")
  filename <- paste0('Results/CREAT_Hazard_pictures/CREAT_', CREAT_SUBj, '.png')
  ggsave(filename)
}

patints_smk_0 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$SMKSTAT == 0]
patints_smk_1 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$SMKSTAT == 1]
patints_smk_2 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$SMKSTAT == 2]
ped_data1upd$SNKSTAT <- as.factor(ped_data1upd$SMKSTAT)
for (SMK0_SUBj in patints_smk_0){
  for (SMK1_SUBj in patints_smk_1) {
    for (SMK2_SUBj in patints_smk_2) {
      plotting_ped <- ped_data1upd[ped_data1upd$USUBJID %in% c(SMK0_SUBj, SMK1_SUBj, SMK2_SUBj), ]
      plotting_ped$SMKSTAT <- as.factor(plotting_ped$SMKSTAT)
      ggplot(plotting_ped, aes(x=tend)) +
        geom_line( aes(y=hazard, color=SMKSTAT)) +
        geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper, fill= SMKSTAT), alpha = 0.15) +
        theme_minimal() +
        xlab('Time') +
        ylab('Hazard') + 
        theme(
          axis.title.y = element_text(size=16),
          axis.title.y.right = element_text(size=16),
          title = element_text(size = 18)
        ) + scale_color_discrete(name = 'Smoking Status', labels =  c('Never Smoker', 'Former Smoker', 'Habital Smoker'))+
        scale_fill_discrete(guide="none")
      filename <- paste0('Results/SMKSTAT_Hazard_pictures/SMKSTAT_', SMK0_SUBj, SMK1_SUBj, SMK2_SUBj, '.png')
      ggsave(filename)
    }
  }
}

patints_who_0 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$WHOSTATN == 0]
patints_who_1 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$WHOSTATN == 1]
patints_who_2 <- unique_measurments$USUBJID[unique_measurments$tenduni > 450 & unique_measurments$WHOSTATN == 2]
ped_data1upd$WHOSTATN <- as.factor(ped_data1upd$WHOSTATN)
for (WHO0_SUBj in patints_who_0){
  for (WHO1_SUBj in patints_who_1) {
    for (WHO2_SUBj in patints_who_2) {
      plotting_ped <- ped_data1upd[ped_data1upd$USUBJID %in% c(WHO0_SUBj, WHO1_SUBj, WHO2_SUBj), ]
      plotting_ped$WHOSTATN <- as.factor(plotting_ped$WHOSTATN) 
      ggplot(plotting_ped, aes(x=tend)) +
        geom_line(aes(y=hazard, color=WHOSTATN)) +
        geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper, fill= WHOSTATN), alpha = 0.15) +
        theme_minimal() +
        xlab('Time') +
        ylab('Hazard') + 
        theme(
          axis.title.y = element_text(size=16),
          axis.title.y.right = element_text(size=16),
          title = element_text(size = 18)
        ) + scale_color_discrete(name = 'Performance Status') +
        scale_fill_discrete(guide="none")
      filename <- paste0('Results/WHOSTAT_Hazard_pictures/WHOSTAT_', WHO0_SUBj, WHO1_SUBj, WHO2_SUBj, '.png')
      ggsave(filename)
    }
  }
}
