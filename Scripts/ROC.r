n_train   <- 300
data1_tte <- data1_tte[data1_tte$USUBJID %in% unique(ped_data1upd$USUBJID), ]
train_idx <- sample(data1_tte$USUBJID, n_train)

train_ped <- as_ped(
  data    = list(data1_tte[data1_tte$USUBJID %in% train_idx, ], data1hr[data1hr$USUBJID %in% train_idx, ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID")
test_ped <- as_ped(
  data    = list(data1_tte[!(data1_tte$USUBJID %in% train_idx), ], data1hr[!(data1hr$USUBJID %in% train_idx), ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, ALT, tz_var = "TIME"),
  id      = "USUBJID")

pam1 <- gam(ped_status ~ s(tend) + SMKSTAT + WHOSTATN,
            data = train_ped,
            family = poisson(), offset = offset, 
            method = 'REML')
pam2 <- gam(ped_status ~ s(tend) + WHOSTATN + SMKSTAT + ALP + ALT + SLD + s(CREAT),
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
}
