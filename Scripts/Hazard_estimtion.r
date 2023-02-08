source('Scripts/dependences.r')
data <- read.csv('SourceData/progression_subj_163.csv')
ped <- data %>% as_ped(Surv(PRGDTC, PRGEVENT)~., id = 'USUBJID')
data_prog <- data
data <-read.csv('SourceData/megred_table_JM.csv')
data_gam<- data3[, c('USUBJID', 'ALP', 'AST', 'CREAT', 'LDH', 'WBC', 'NEUT', 'NLR', 'TIME')]
data_gam <- data_gam[data_gam$USUBJID %in% data_prog$USUBJID, ]
#Baseline hazard estimation with gam 

#TDC with gam model
sbj <- data_gam$USUBJID[data_gam$TIME == 0]
sbj <- unlist(data_gam$USUBJID[data_gam$TIME == 0])
data_gam <- data_gam[data_gam$USUBJID %in% sbj, ]
data_prog <- data_prog[data_prog$USUBJID %in% sbj, ]
ped_pbc <- as_ped(
data    = list(data_prog, data_gam),
formula = Surv(PRGDTC, PRGEVENT) ~ . + concurrent(LDH, tz_var = "TIME") +
  + concurrent(NLR, tz_var = "TIME"),
id      = "USUBJID")
pam <- gam(ped_status ~ s(tend) + s(LDH) , data = ped_pbc,
family = poisson(), offset = offset)
ped_pbc <- ped_pbc %>% add_hazard(pam)

ggplot(ped_pbc, aes(x = tend)) +
geom_stephazard(aes(y = hazard)) +
geom_stepribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
ylab(expression(hat(lambda)(t))) + xlab(expression(t))  + ggtitle('Hazrd Function')



#baseline hazard 
base_df <- basehaz(coxph(Surv(PRGDTC, PRGEVENT)~1, data = data_prog)) %>%
  rename(nelson_aalen = hazard)
ggplot(base_df, aes(x = time, y = nelson_aalen)) +
  geom_stephazard() +
  ylab(expression(hat(Lambda)(t))) + xlab("t") +
  ggtitle("Nelson-Aalen estimate of the cumulative hazard")
pam <- gam(ped_status ~ s(tend), data = ped, offset = offset, family = poisson())
summary(pam)
int_df <- make_newdata(ped, tend = unique(tend)) %>%
  add_cumu_hazard(pam)

gg_baseline <- ggplot(int_df, aes(x = tend)) +
  geom_line(aes(y = cumu_hazard, col = "PAM")) +
  geom_stephazard(data = base_df, aes(x=time, y = nelson_aalen, col = "Nelson-Aalen")) +
  scale_color_manual(
    name   = "Method",
    values = c("PAM" = "black", "Nelson-Aalen" = Set1[1])) +
  theme(legend.position = "bottom") +
  ylab(expression(hat(Lambda)(t))) + xlab("t") +
  ggtitle(paste0("Comparison of cumulative hazards estimated by\n",
                 "Cox-PH (Nelson-Aalen) vs. PAM"))
#TDC 

#TDC with gam
time <- c(1:400)
subj <- unique(unlist(data_gam$USUBJID))


data_time <- data.frame(TIME = rep(c(1:400), times = length(subj), USUBJID = rep(subj, each = 400)))
data_time$USUBJID = rep(subj, each = 400)

ct_ALP <- gam(ALP ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$ALT <- predict(ct1, data_time)
ct_AST <- gam(AST ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
ct_AST <- gam(AST ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$AST <- predict(ct_AST, data_time)
ct_CREAT <- gam(CREAT ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$CREAT <- predict(ct_CREAT, data_time)
ct_LDH <- gam(LDH ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$LDH <- predict(ct_LDH, data_time)
ct_NLR <- gam(NLR ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$NLR <- predict(ct_NLR, data_time)
ct_NEUT <- gam(NEUT ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$NEUT <- predict(ct_NEUT, data_time)
ct_WBC <- gam(WBC ~ s(TIME) + s(USUBJID, bs = 're'), data = data_gam)
data_time$WBC <- predict(ct_WBC, data_time)
colnames(data_time)[1] <- 'tend'
data_final <- merge(data_time, ped, by = c('USUBJID', 'tend'))
pam2 <- gam(ped_status ~ s(tend) + s(LDH) + s(NLR),
            data = data_final, offset = offset, family = cox.ph())
data_final$hazard <- predict(pam2, data_final)
data_final$hazard <- exp(data_final$hazard)
data_final$cum_hazard <- cumsum(data_final$hazard) / sum(data_final$hazard)
ggplot(aes(x=tend, y=hazard), data=data_final) + geom_line() + theme_minimal() + 
  xlab('Time, days') + ylab('Hazard')
ggplot(aes(x=tend, y=cum_hazard), data=data_final) + geom_scatter()