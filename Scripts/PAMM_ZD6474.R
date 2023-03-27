setwd('D:/study/MSD/GAMproj/')
source('Scripts/dependences.r')
theme_set(theme_bw())
Set1 <- RColorBrewer::brewer.pal(9, "Set1")



#######  Read Data

data <- read.csv("SourceData/megred_table_JM.csv")
head(data)

data1 <- data[data$STUDY == "ZD6474",]
head(data1)
data1 <- get_zero_time(data1)

data1_tte <- data1[!duplicated(data1$USUBJID),]
head(data1_tte)

#######  Plot Survival
fit1 <- survfit(Surv(time=EVENT_TIME, event=EVENT) ~ 1, data=data1_tte)
print(fit1)
plot(x=c(0,850), y=c(0,1), type = 'n', xlab = "Time,days", ylab = "PFS" )
lines(fit1, col="red", lwd = 2)
######################

##### Make dataset with TDCs

data1s <- data1[, c("USUBJID", "TIME","STUDY","LINE")]
data1hr <- data1s[!duplicated(data1s),]


for (i in 1:9){
  data1si <- data1[data1$YTYPE == i, c("USUBJID", "TIME","DV","YTYPE_DESC")]
  names(data1si)[3] <- data1si$YTYPE_DESC[1]
  data1hr <- merge(data1hr, data1si[,1:3], by = c("USUBJID", "TIME") , all = T)
}
# head(data1hr)
subjects <- unique(data1hr$USUBJID)
data1hr <- data1hr[
  with(data1hr, order(USUBJID, TIME)),
]
###### Fill NAs wit previous mesurment of the same subject
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
  if (is.na(data1hr$LDH[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$LDH[i+1] <- data1hr$LDH[i]
  }
  if (is.na(data1hr$NLR[i+1]) & (data1hr$USUBJID[i] == data1hr$USUBJID[i+1])) {
    data1hr$NLR[i+1] <- data1hr$NLR[i]
  }
}


#### combine TTE and TDC data
data1_tte <-  data1[!duplicated(data1$USUBJID),c(1:5,9:11,18)]


ped_data1 <- as_ped(
  data    = list(data1_tte, data1hr),
  
  formula = Surv(time=EVENT_TIME, event=EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, AST, LDH, NEUT, NLR, WBC, tz_var = "TIME"), 
  id      = "USUBJID")

head(ped_data1)

summary(ped_data1)
##############

####### Exclude subjects with copletely missing biomarkers

stat <- ped_data1 %>% group_by(USUBJID) %>% summarize(SLDna = sum(is.na(SLD)), ASTna = sum(is.na(AST)), 
                                                      ALPna = sum(is.na(ALP)), ALTna = sum(is.na(ALT)), 
                                                      CREATna = sum(is.na(CREAT)), WBCna = sum(is.na(WBC)), 
                                                      NEUTna = sum(is.na(NEUT)), LDHna = sum(is.na(LDH)),
                                                      NLRna = sum(is.na(NLR)), count = n())
completely_missing <- unique(stat$USUBJID[stat$CREATna == stat$count])
completely_missing <- c(completely_missing, unique(stat$USUBJID[stat$ALPna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$ALTna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$ASTna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$SLDna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$WBCna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$NEUTna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$NLRna == stat$count]))
completely_missing <- c(completely_missing,  unique(stat$USUBJID[stat$LDHna == stat$count]))
completely_missing <- unique(completely_missing)
ped_data1upd <- ped_data1
ped_data1upd <- ped_data1upd[!(ped_data1upd$USUBJID) %in% completely_missing,  ]
ped_data1upd$intlen <- ped_data1upd$tend - ped_data1upd$tstart
BM_by_id <-  ped_data1upd %>% group_by(USUBJID) %>%
  summarise( SLDav = weighted.mean(SLD,intlen,  na.rm = T), ALPav = weighted.mean(ALP, intlen, na.rm = T), 
             CREATav = weighted.mean(CREAT,intlen, na.rm = T), ALPmax =max(ALP, na.rm = T), 
             CREATmax = max(CREAT, na.rm = T), ASTav = weighted.mean(AST,intlen, na.rm = T),
             ALTav = weighted.mean(ALT, intlen, na.rm = T), WBCav = weighted.mean(WBC, intlen, na.rm = T),
             NEUTav = weighted.mean(NEUT, intlen, na.rm = T), NLRav = weighted.mean(NLR, intlen, na.rm = T),
             LDHav = weighted.mean(LDH, intlen, na.rm = T))
ped_data1upd <- merge(ped_data1upd,BM_by_id, by= "USUBJID", all = T)
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
  if (is.na(ped_data1upd$WBC[i])) {
    ped_data1upd$WBC[i] <- ped_data1upd$WBCav[i]
  }
  if (is.na(ped_data1upd$WBC[i])) {
    ped_data1upd$NLR[i] <- ped_data1upd$NLRav[i]
  }
  if (is.na(ped_data1upd$WBC[i])) {
    ped_data1upd$LDH[i] <- ped_data1upd$LDHav[i]
  }
}

##############################
###### ROC calculation for best model 
set.seed(1)
n_train   <- 300
data1_tte <- data1_tte[data1_tte$USUBJID %in% unique(ped_data1upd$USUBJID), ]
train_idx <- sample(data1_tte$USUBJID, n_train)

train_ped <- as_ped(
  data    = list(data1_tte[data1_tte$USUBJID %in% train_idx, ], data1hr[data1hr$USUBJID %in% train_idx, ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, NLR, LDH, NEUT, tz_var = "TIME"),
  id      = "USUBJID")



###### Alongation of the dataset for test
d1 <- ped_data1upd[ped_data1upd$tend == max(ped_data1upd$tend),]
ped_long <-  ped_data1upd[ped_data1upd$USUBJID == d1[1,"USUBJID"],]
ped_long[,7:21] <- NA
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
    
  } else{
    pdatilong <- pdati
  }
  ped_data1_long <- rbind(ped_data1_long, pdatilong)
}

test_ped <- ped_data1_long[!(ped_data1_long$USUBJID %in% train_idx),]
# test_ped <- ped_data1_long[!(ped_data1_long$USUBJID %in% completely_missing), ]
test_ped$intlen <- test_ped$tend - test_ped$tstart

pam1 <- gam(ped_status ~ s(tend) + SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam2 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + s(NLR),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam3 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + s(NLR) + s(SLD) ,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam4 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + s(NLR) + s(SLD) + s(LDH),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam5 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + s(NLR) + s(SLD) + s(LDH) + s(NEUT),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
test_ped_new1 <- predict_surv(test_ped, pam1)
test_ped_new2 <- predict_surv(test_ped, pam2)
test_ped_new3 <- predict_surv(test_ped, pam3)
test_ped_new4 <- predict_surv(test_ped, pam4)
test_ped_new5 <- predict_surv(test_ped, pam5)

AUCs <- c()
times <- c(100, 200, 300, 400, 500, 600, 700, 800)
models <- c('hazard~s(tend) + SMKSTAT + WHOSTATN\n',
            'hazard~s(tend) + SMKSTAT + WHOSTATN \n+ s(NLR)',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n s(NLR) + s(SLD)',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n s(NLR) + s(SLD) + s(LDH)',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n s(NLR) + s(SLD) + s(LDH) + s(NEUT)'
)
colors <- hue_pal()(5)
i_models <- c(1:5)
for (timepoint in times){
  mainLab <- sprintf("Time %d", timepoint)
  filename <- paste0('Results/ROC_curves/ZD6474/ZD6474', timepoint, length(i_models), '.png')
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
                       models = rep(models[i_models], length(i_models)))
ROC_AUCs$models <- as.factor(ROC_AUCs$models)
ggplot(data =ROC_AUCs) + geom_line(aes(x = Time, y = AUC,color = models), size = 1) +theme_bw()
ggsave('Results/ROC_curves/ZD6474/all5_in_time_ZD6474.png')

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
ggsave(p_bs, file =  "Results/BS/ZD6474_BS_time.png", width = 6, height = 4)

###### PEC calculation 


merged <- merge(data1_tte[, c(1:9)], data1hr, by = c('USUBJID', 'LINE', 'STUDY'), all =  F)
merged <- merge(merged, BM_by_id, by = 'USUBJID')

merged <- merged[!(merged$USUBJID) %in% completely_missing,  ]
for (i in c(1:nrow(merged))) {
  if (is.na(merged$NEUT[i])) {
    merged$NEUT[i] <- merged$NEUTav[i]
  }
  if (is.na(merged$SLD[i])) {
    merged$SLD[i] <- merged$SLDav[i]
  }
  if (is.na(merged$NLR[i])) {
    merged$NLR[i] <- merged$NLRav[i]
  }
  if (is.na(merged$LDH[i])) {
    merged$LDH[i] <- merged$LDHav[i]
  }
}
summary(merged)
pam6 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN,
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam7 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + s(NLR),
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam3 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + s(NLR) + s(SLD),
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam4 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + s(NLR) + s(SLD) + s(LDH),
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')
pam5 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + s(NLR) + s(SLD) + s(LDH) + s(NEUT),
             data = train_ped,
             family = poisson(), offset = offset, 
             method = 'REML')

colnames(merged)[1] <- 'id'

pec <- pec(
  list(pam1 = pam6, pam2 = pam7),
  Surv(EVENT_TIME, EVENT) ~ 1, # formula for IPCW
  data = merged,
  times = seq(.01, 400, by = 10),
  start = .01,
  exact = FALSE,
  reference = T,
  cens.model = 'aalen'
)
#attr(pec[["models"]], "names") <- c("KM", models)
plot(pec)

#Saving data for cross validation
# merged_ZD6474 <- merged
ped_data_long_ZD6474 <- ped_data1_long
ped_data1_ZD6474 <- ped_data1upd
# write.csv(merged_ZD6474, 'DerivedData/merged_ZD6474.csv')
write.csv(ped_data_long_ZD6474, 'DerivedData/ped_long_ZD6474.csv')
write.csv(ped_data1_ZD6474, 'DerivedData/ped_ZD6474.csv')
