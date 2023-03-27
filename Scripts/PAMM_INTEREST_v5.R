setwd('D:/study/MSD/GAMproj/')
source('Scripts/dependences.r')
theme_set(theme_bw())
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
  if (is.na(data1hr$NEUT[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
     data1hr$NEUT[i+1] <- data1hr$NEUT[i]
  }
  if (is.na(data1hr$WBC[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$WBC[i+1] <- data1hr$WBC[i]
  }
}
BM_by_id <-  data1hr %>% group_by(USUBJID) %>%
  summarise( SLDav = mean(SLD, na.rm = T), ALPav = mean(ALP, na.rm = T), 
             CREATav = mean(CREAT, na.rm = T), ALPmax =max(ALP, na.rm = T), 
             CREATmax = max(CREAT, na.rm = T), ASTav = mean(AST, na.rm = T),
             ALTav = mean(ALT, na.rm = T), NEUTav = mean(NEUT, na.rm = T))

#### combine TTE and TDC data
data1_tte <-  data1[!duplicated(data1$USUBJID),c(1:5,9:11,18)]

ped_data1 <- as_ped(
  data    = list(data1_tte, data1hr),
  cut = seq(0,1000, by =7),
  formula = Surv(time=EVENT_TIME, event=EVENT) ~ . + concurrent(SLD, CREAT,ALP, ALT, AST, NEUT, WBC, tz_var = "TIME"), 
  id      = "USUBJID")


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
  if (is.na(ped_data1upd$NEUT[i])) {
    ped_data1upd$NEUT[i] <- ped_data1upd$NEUTav[i]
  }
}
pam_strata1 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + s(SLD) + WBC + s(LDH) + AST + s(CREAT) + s(ALT),
                   data = ped_data1upd,
                   family = poisson(), offset = offset, 
                   method = 'REML')
summary(pam_strata1)
AIC(pam_strata1)


##############################
set.seed(1)
n_train   <- 300
train_idx <- sample(data1_tte$USUBJID, n_train)
train_ped <- ped_data1upd[ped_data1upd$USUBJID %in% train_idx, ]

d1 <- ped_data1upd[ped_data1upd$tend == max(ped_data1upd$tend),]
ped_long <-  ped_data1upd[ped_data1upd$USUBJID == d1[1,"USUBJID"],]
ped_long[,7:19] <- NA
ped_data1_ids <- ped_data1upd[!duplicated(ped_data1upd$USUBJID),]
ped_data1_long <- data.frame()
for (i in 1:nrow(ped_data1_ids)){
  idi <-  ped_data1_ids[i,"USUBJID"]
  pdati <-  ped_data1upd[ped_data1upd$USUBJID == idi,]
  if (max(pdati$tend) < max(ped_long$tend)){
    pdati_upd <- ped_long[ped_long$tend > max(pdati$tend),]
    pdati_upd$USUBJID <- idi
    pdati_upd$ped_status <- pdati$ped_status[pdati$tend == max(pdati$tend)]
    
    pdati_long <- rbind(pdati, pdati_upd)
    pdatilong <- na.locf(na.locf(pdati_long),fromLast=TRUE)
    
  }
  else{
    pdatilong <- pdati
  }
  ped_data1_long <- rbind(ped_data1_long, pdatilong)
}

test_ped <- ped_data1_long[!(ped_data1_long$USUBJID %in% train_idx),]
pam1 <- gam(ped_status ~ log(tend) + SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam2 <- gam(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam3 <- gam(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP + SLD,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam4 <- gam(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP + ALT + SLD,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam5 <- gam(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP + ALT + SLD + s(CREAT),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
test_ped$intlen <- test_ped$tend - test_ped$tstart
test_ped_new1 <- predict_surv(test_ped, pam1)
test_ped_new2 <- predict_surv(test_ped, pam2)
test_ped_new3 <- predict_surv(test_ped, pam3)
test_ped_new4 <- predict_surv(test_ped, pam4)
test_ped_new5 <- predict_surv(test_ped, pam5)
AUCs <- c()
times <- c(100, 200, 300, 400, 500, 600, 700, 800)
models <- c('hazard~s(tend) + SMKSTAT + WHOSTATN\n',
            'hazard~s(tend) + SMKSTAT + WHOSTATN \n+ ALT',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD + ALP',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD + ALP + s(CREAT)'
)
colors <- hue_pal()(5)
i_models <- c(1:5)
for (timepoint in times){
  mainLab <- sprintf("Time %d", timepoint)
  filename <- paste0('Results/ROC_curves/INTEREST/INTEREST_cv_ZODIAC', timepoint, length(i_models), '.png')
  png(file=filename)
  plot(x = c(0, 1), y = c(0, 1), type = 'n', xlab = '1 - Sensetivity', ylab = 'Specificity', 
       main=mainLab)
  j <- 0
  for (i in i_models) {
    ds_name <- paste0('test_ped_new', i)
    AUCs <- c(AUCs, plot_roc(timepoint, get(ds_name), x = 0.6, y = 0.45 - (j - 1) * 0.1, color = colors[i], 
                             model = models[i]))
    j <- j + 1
  }
  dev.off()
}

ROC_AUCs <- data.frame(AUC = AUCs, Time = rep(times, each = length(i_models)), 
                       models = rep(models[i_models], length(times)))
ROC_AUCs$models <- as.factor(ROC_AUCs$models)
ggplot(data =ROC_AUCs) + geom_line(aes(x = Time, y = AUC,color = models), size = 1) +theme_bw()
ggsave('Results/ROC_curves/INTEREST/all_in_time_INTEREST_cv_ZODIAC.png')


BS_t_sum <- c()
for (timepoint in times){  
  for (i in i_models) {
    ds_name <- paste0('test_ped_new', i)
    BS_t_sum <- c(BS_t_sum, calc_BS(timepoint, get(ds_name)))
  }
}
BS_ds <- data.frame(BS = BS_t_sum, Time = rep(times, each = length(i_models)), 
                    models = rep(models[i_models], length(times)))
p_bs <- ggplot(BS_ds, aes(Time, BS, col=models)) +  geom_line()+
  ggtitle("BS dynamics")

p_bs
ggsave(p_bs, file =  "Results/BS/INTEREST_cv_ZODIAC_BS_time.png", width = 6, height = 4)


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
  if (is.na(merged$NEUT[i])) {
    merged$NEUT[i] <- merged$NEUTav[i]
  }
}
colnames(merged)[1] <- 'id'
pamm1 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN,
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam2 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP,
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam3 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP + ALT,
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam4 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP + ALT + SLD,
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam5 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP + ALT + SLD + s(CREAT),
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')


pec <- pec(
  list(pam1 = pam1),
  Surv(EVENT_TIME, EVENT) ~ 1, # formula for IPCW
  data = merged[!(merged$id %in% train_idx), ],
  times = seq(.1, 20, by = 10),
  start = .01,
  exact = FALSE,
  reference = T,
  cens.model = 'aalen'
)
plot(pec)
ggsave('Results/PEC.png')

ped_data_long_INTEREST <- ped_data1_long
ped_data1_INTEREST <- ped_data1upd
# write.csv(merged_INTEREST, 'DerivedData/merged_INTEREST.csv')
write.csv(ped_data_long_INTEREST, 'DerivedData/ped_long_INTEREST.csv')
write.csv(ped_data1_INTEREST, 'DerivedData/ped_INTEREST.csv')

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
