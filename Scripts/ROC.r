n_train   <- 300
data1_tte <- data1_tte[data1_tte$USUBJID %in% unique(ped_data1upd$USUBJID), ]
train_idx <- sample(data1_tte$USUBJID, n_train)
train_ped <- ped_data1[ped_data1$USUBJID %in% train_idx, ]
train_ped <- as_ped(
  data    = list(data1_tte[data1_tte$USUBJID %in% train_idx, ], data1hr[data1hr$USUBJID %in% train_idx, ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID")
data1_test_tte <- data1_tte
data1_test_tte$EVENT_TIME <- 994
test_ped <- as_ped(
  data    = list(data1_test_tte[!(data1_test_tte$USUBJID %in% train_idx), ], data1hr[!(data1hr$USUBJID %in% train_idx), ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID")

pam1 <- gam(ped_status ~ s(tend), #+ SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam2 <- pamm(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + ALP + ALT + SLD + s(CREAT),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')

test_ped$intlen <- test_ped$tend - test_ped$tstart

test_ped <- test_ped %>% add_hazard(pam1, overwrite = T)

roc_curve <- roc(test_ped$ped_status, test_ped$hazard, auc.polygon=TRUE)


test_ped <- test_ped %>% add_hazard(pam2, overwrite = T)

roc_curve1 <-roc(test_ped$ped_status, test_ped$hazard, auc.polygon=TRUE)
plot(roc_curve)
plot(roc_curve1)
for_plotting <- data.frame(Model = 'Hazard~s(tend) + SMKSTAT + WHOSTATN',
                           Sensitivity = roc_curve[["sensitivities"]], 
                           Specificity = roc_curve[["specificities"]])
for_plotting1 <- data.frame(Model = 'Hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALP + ALT + SLD +s(CREAT)',
                            Sensitivity = roc_curve1[["sensitivities"]], 
                            Specificity = roc_curve1[["specificities"]])
for_plotting <- rbind(for_plotting, for_plotting1)
for_plotting1 <- data.frame(Model = 'Random Classifier',
                            Sensitivity = c(1, 0), 
                            Specificity = c(0, 1))
for_plotting <- rbind(for_plotting, for_plotting1)
ggplot(aes(x=Specificity, y = Sensitivity, color = Model), data = for_plotting) + geom_line()
ggsave('Results/ROC_all.png')
sbjs <- test_ped$USUBJID
test_ped_new <- data.frame()
for (sbj in sbjs) {
  test_ped_i <- test_ped[test_ped$USUBJID == sbj, ]
  if (nrow(test_ped_i) > 0) {
  test_ped_i <- test_ped_i %>% add_surv_prob(pam1)
  test_ped_new <- rbind(test_ped_new, test_ped_i[,c("USUBJID","tend","surv_prob", "ped_status")])
  }
}
roc_curve <- roc(test_ped_new$ped_status, test_ped_new$surv_prob, auc.polygon=TRUE)
test_ped_new <- data.frame()
for (sbj in sbjs) {
  test_ped_i <- test_ped[test_ped$USUBJID == sbj, ]
  if (nrow(test_ped_i) > 0) {
    test_ped_i <- test_ped_i %>% add_surv_prob(pam2)
    test_ped_new <- rbind(test_ped_new, test_ped_i)
  }
}
roc_curve1 <- roc(test_ped_new$ped_status, test_ped_new$surv_prob, auc.polygon=TRUE)
# for (time_point in seq(25, 1000, 25)) {
time_point <- 200
# test_ped <- test_ped %>% add_hazard(pam1, overwrite = T)
test_ped_timepoint <- test_ped_new[test_ped_new$tend < time_point,]
for_auc <- test_ped_timepoint %>% 
    group_by(USUBJID) %>%
    summarize(status = ifelse(sum(ped_status) > 1, 0, 1),
              predicts = last(surv_prob))

roc_curve1 <- roc(for_auc$status, for_auc$predicts, auc.polygon=TRUE)
test_ped <- test_ped %>% add_hazard(pam2, overwrite = T)
test_ped_timepoint <- test_ped[test_ped$tend < time_point,]
for_auc <- test_ped_timepoint %>% 
  group_by(USUBJID) %>%
  summarize(status = sum(ped_status),
            predicts = sum(hazard))
roc_curve1 <- roc(for_auc$status, for_auc$predicts, auc.polygon=TRUE)
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
title <- paste0('Time: ', time_point)
filename <- paste0('Results/AUC_time_', time_point, '.png')
ggplot(aes(x=Specificity, y = Sensitivity, color = Model), data = for_plotting) + geom_line() + ggtitle(title)
ggsave(filename)
# }

base <- coxph(Surv(EVENT_TIME, EVENT)~1, data1_tte)
head(base)
tr <- timeROC(data1_tte$EVENT_TIME[!(data1_tte$USUBJID %in% train_idx)], 
              data1_tte$EVENT[!(data1_tte$USUBJID %in% train_idx)], 
              data1_tte$WHOSTATN[!(data1_tte$USUBJID %in% train_idx)],
              cause = 1, times = c(200, 400, 600, 800))
pam1 <- gam(ped_status ~ log(tend) + WHOSTATN, #+ SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
test_ped$intlen <- test_ped$tend - test_ped$tstart

test_ped_new <- data.frame()
for (sbj in unique(test_ped$USUBJID)) {
  test_ped_i <- test_ped[test_ped$USUBJID == sbj, ]
  if (nrow(test_ped_i) > 0) {
    test_ped_i <- test_ped_i %>% add_surv_prob(pam1)
    test_ped_new <- rbind(test_ped_new, test_ped_i)
  }
}
simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  return(sum(TPR * dFPR) + sum(dTPR * dFPR)/2)
}


for (time_point in c(200, 400, 600, 800)) {
test_ped_timepoint <- test_ped_new[test_ped_new$tend < time_point,]
test_ped_timepoint<- test_ped_timepoint[
  with(test_ped_timepoint, order(USUBJID, tend)),
]
for_auc <- test_ped_timepoint %>% 
  group_by(USUBJID) %>%
  summarize(status = ifelse(sum(ped_status) > 1, 1, 0),
            predicts = last(surv_prob))
senss <- c(0)
specs <- c(0)
for (tre in sort(for_auc$predicts)) {
  for_auc$res <- ifelse(for_auc$predicts <= tre, 1, 0)
  TP <- sum(for_auc$res * for_auc$status)
  TN <- sum((!for_auc$res) * (!for_auc$status))
  FP <- sum(for_auc$res) - TP
  FN <- sum(!for_auc$res) - TN 
  senss <- c(senss, TP/(TP+FN))
  specs <- c(specs, 1 - TN/(TN + FP))
}
AUC <- simple_auc(senss, specs)
# AUC <- 0
# denom <- 0
# for (i in c(1:nrow(for_auc))) {
#   for (j in c(1:nrow(for_auc))) {
#     AUC <- AUC + ind_1(for_auc$predicts[i], for_auc$predicts[j]) * ind_2(for_auc$status[i], for_auc$status[j])
#     denom <- denom + ind_1(for_auc$predicts[i], for_auc$predicts[j])
#   }
# }
# AUC <- AUC / denom

main_lab <- sprintf("AUC: %.4f Time: %d",AUC, time_point)

filename = paste0('Results/ROC_curves/ROC_', time_point, '_', 'WHOSTATN.png')
png(file=filename)
plot(specs, senss, type = 'l', xlab = '1 - Speceficity',
     ylab = 'Sensetivity', main = main_lab)
dev.off()
}
plot_roc <- function(time_point, dataset, color, x, y, model) {
  test_ped_timepoint <- dataset[dataset$tend < time_point,]
  test_ped_timepoint<- test_ped_timepoint[
    with(test_ped_timepoint, order(USUBJID, tend)),
  ]
  for_auc <- test_ped_timepoint %>% 
    group_by(USUBJID) %>%
    summarize(status = ifelse(sum(ped_status) > 1, 1, 0),
              predicts = last(surv_prob))
  senss <- c(0)
  specs <- c(0)
  for (tre in sort(for_auc$predicts)) {
    for_auc$res <- ifelse(for_auc$predicts <= tre, 1, 0)
    TP <- sum(for_auc$res * for_auc$status)
    TN <- sum((!for_auc$res) * (!for_auc$status))
    FP <- sum(for_auc$res) - TP
    FN <- sum(!for_auc$res) - TN 
    senss <- c(senss, TP/(TP+FN))
    specs <- c(specs, 1 - TN/(TN + FP))
  }
  AUC <- simple_auc(senss, specs)
  print(color)
  print(x)
  print(y)
  print(model)
  Lab <- sprintf("%s AUC: %.4f", model, AUC)
  lines(specs, senss, col = color)
  text(x=x, y=y, label = Lab, col = color)
  return(AUC)
}
predict_surv <- function(test_ped, model) {
  test_ped_new <- data.frame()
  for (sbj in unique(test_ped$USUBJID)) {
    test_ped_i <- test_ped[test_ped$USUBJID == sbj, ]
    if (nrow(test_ped_i) > 0) {
      test_ped_i <- test_ped_i %>% add_surv_prob(model)
      test_ped_new <- rbind(test_ped_new, test_ped_i)
    }
  }
  return(test_ped_new)
}
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
times <- c(200, 400, 600, 800)
models <- c('hazard~s(tend) + SMKSTAT + WHOSTATN\n',
            'hazard~s(tend) + SMKSTAT + WHOSTATN \n+ ALT',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD + ALP',
            'hazard~s(tend) + SMKSTAT + WHOSTATN +\n ALT + SLD + ALP + s(CREAT)'
            )
colors <- brewer.pal(n = 5, name = "Set1")
for (timepoint in times){
  mainLab <- sprintf("Time %d", timepoint)
  filename <- paste0('Results/ROC_curves/INTEREST_ZODIAC', timepoint, '5.png')
  png(file=filename)
  plot(x = c(0, 1), y = c(0, 1), type = 'n', xlab = '1 - Sensetivity', ylab = 'Specificity', 
       main=mainLab)
  
  for (i in c(1:5)) {
    ds_name <- paste0('test_ped_new', i)
    AUCs <- c(AUCs, plot_roc(timepoint, get(ds_name), x = 0.6, y = 0.55 - (i - 1) * 0.05, color = colors[i], 
                             model = models[i]))
  }
  dev.off()
}

ROC_AUCs <- data.frame(AUC = AUCs, Time = rep(times, each = 5), models = rep(models, 4))
ROC_AUCs$models <- as.factor(ROC_AUCs$models)
ggplot(data =ROC_AUCs) + geom_line(aes(x = Time, y = AUC,color = models), size = 1)
ggsave('Results/ROC_curves/all_in_time_INTEREST_ZODIAC.png')


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
pam1 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam2 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam3 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP + SLD,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam4 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN + ALP + ALT + SLD,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam5 <- pamm(ped_status ~ log(tend) + SMKSTAT + WHOSTATN+ ALP + ALT + SLD + s(CREAT),
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pec <- pec(
  list(pam1 = pam1, pam2 = pam2, pam3 = pam3, pam4 = pam4, pam5 = pam5),
  Surv(EVENT_TIME, EVENT) ~ 1, # formula for IPCW
  data = merged[!(merged$id %in% train_idx), ],
  times = seq(.01, 400, by = 10),
  start = .01,
  exact = FALSE,
  reference = T,
  cens.model = 'aalen'
)
plot(pec)
ggsave('Results/PEC.png')
