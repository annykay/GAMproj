# 'exclude' column in final dataset:
# 0 - don't eclude
# 1 - exclude for sure
# 2 - exclude in the end to
# 3 - duplicates from time adding

source('Scripts/dependences.r')
source('Scripts/dataTransform.r')
source('Scripts/plotFunc.r')
source('Scripts/statistics.r')

data_all <- read.csv('SourceData/megred_table_JM.csv')
data_sets <- c()
studies <- sort(unique(data_all$STUDY))

study <- studies[1]

data <- data_all[data_all$STUDY == study, ]
data <- dfFactorize(data, c('SMKSTAT', 'WHOSTATN',  'EGFRMUTN'))
data <- dfOut(data)
data1 <- dfNormalize(data)
markers <- sort(unique(data1$YTYPE_DESC))[-8]
data1 <- chooseTransform(data1, 'DV', rep('DV', length(markers)), markers)

# filename <- paste0('DerivedData/long', study, '.xlsx')
# write_xlsx(data1, filename)
# 
# markers <- sort(unique(data1$YTYPE_DESC))
# statistic <- basic_biomarker_statistics(data1[!(data1$exclude %in% c(1, 2, 3)), ], markers)
# filename <- paste0('./DerivedData/MarkerStat', study, '_WO_OUTLIERS.xlsx')
# write_xlsx(statistic, filename, format_headers = TRUE)
# 
# statistic <- basic_biomarker_statistics(data1, markers)
# filename <- paste0('./DerivedData/MarkerStat', study, '.xlsx')
# write_xlsx(statistic, filename, format_headers = TRUE)
# 
# 
# factor_statistic <- factor_stat(c('SMKSTAT', 'WHOSTATN', 'EGFRMUTN'), data1[!(data1$exclude %in% c(1, 3)), ])
# filename <- paste0('./DerivedData/FactorStat', study, '.xlsx')
# write_xlsx(factor_statistic, filename, format_headers = TRUE)
  
data2 <- data1[, c('USUBJID', 'SMKSTAT', 'TIME', 'WHOSTATN',  'EGFRMUTN', 'YTYPE_DESC', 'target', 'exclude', 'base')]
# subjects <- unique(data$USUBJID[!is.na(data$SLD)])
# data2 <- data2[data2$USUBJID %in% subjects, ]
data2 <- baseline(data2, 'SLD')
colnames(data2)[which(colnames(data2) == 'base')] <- 'SLDb'

data2 <- spread(data2, key=YTYPE_DESC, value=target)
  
rows <- nrow(data2)
  
#   
# for (i in c(1:rows)) {
#   if (!is.na(data2$SLD[i])) {
#     subj <- data2$USUBJID[i]
#     time <- data2$TIME[i]
#     data2$SLD[data2$USUBJID == subj & abs(data2$TIME - time) < 4 & is.na(data2$SLD)] <- data2$SLD[i]
#   }
# }
subjects <- unique(data2$USUBJID)
data2 <- data2[
  with(data2, order(USUBJID, TIME)),
]
for (i in c(1:(nrow(data2) - 1))) {
  if (is.na(data2$SLD[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$SLD[i+1] <- data2$SLD[i]
  }
  if (is.na(data2$AST[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$AST[i+1] <- data2$AST[i]
  }
  if (is.na(data2$ALP[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$ALP[i+1] <- data2$ALP[i]
  }
  if (is.na(data2$ALT[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$ALT[i+1] <- data2$ALT[i]
  }
  if (is.na(data2$CREAT[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$CREAT[i+1] <- data2$CREAT[i]
  }
  if (is.na(data2$NEUT[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$NEUT[i+1] <- data2$NEUT[i]
  }
  if (is.na(data2$WBC[i+1]) & (data2$USUBJID[i] == data2$USUBJID[i+1])) {
    data2$WBC[i+1] <- data2$WBC[i]
  }
  
}
# data2$exclude[is.na(data2$AST) & is.na(data2$CREAT) & is.na(data2$NLR) & is.na(data2$LDH)] <- 3
# data2$SLD <- na.locf(na.locf(data2$SLD),fromLast=TRUE)
# data2$AST <- na.locf(na.locf(data2$AST),fromLast=TRUE)
# data2$ALP <- na.locf(na.locf(data2$ALP),fromLast=TRUE)
# data2$CREAT <- na.locf(na.locf(data2$CREAT),fromLast=TRUE)
# data2$NEUT <- na.locf(na.locf(data2$NEUT),fromLast=TRUE)
# data2$WBC <- na.locf(na.locf(data2$WBC),fromLast=TRUE)
# data3 <- data2[!(data2$exclude %in% c(1, 2, 3)), ]
data3$USUBJID <- as.factor(data3$USUBJID)

filename <- paste0('DerivedData/wide', study, '.xlsx')
write_xlsx(data2, filename)

surv <- unique(data_all[data_all$STUDY == study, c('USUBJID', 'EVENT', 'EVENT_TIME', 
                                                   'SLDb', 'WHOSTATN', 'EGFRMUTN', 
                                                   'SMKSTAT')])

get_zero_time <- function(biom) {
  sbj <- biom$USUBJID[biom$TIME == 0]
  sbj <- unlist(biom$USUBJID[biom$TIME == 0])
  sbj1 <- unlist(biom$USUBJID[biom$TIME == 1 & !(biom$USUBJID %in% sbj)])
  biom$TIME[biom$USUBJID %in% sbj1 & biom$TIME == 1] <- 0
  sbj <- biom$USUBJID[biom$TIME == 0]
  biom <- biom[biom$USUBJID %in% sbj, ]
  return(biom)
}
biom <- data2
biom <- biom[biom$USUBJID %in% surv$USUBJID, ]
# biom <- get_zero_time(biom)

surv <- surv[surv$USUBJID %in% biom$USUBJID, ]
surv$EGFRMUTN <- as.factor(surv$EGFRMUTN)
surv$WHOSTATN <- as.factor(surv$WHOSTATN)
surv$SMKSTAT <- as.factor(surv$SMKSTAT)
# plot_hazard(surv, biom, study)

data_event <- as_ped(
  data = list(surv, biom),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(AST, ALT, ALP, CREAT, NEUT, WBC, 
                                                     SLD, tz_var = "TIME") + 
    SLDb + SMKSTAT + WHOSTATN + EGFRMUTN,
  id = "USUBJID")
# data_event$SLD <- na.locf(na.locf(data_event$SLD),fromLast=TRUE)
# data_event$AST <- na.locf(na.locf(data_event$AST),fromLast=TRUE)
# data_event$ALP <- na.locf(na.locf(data_event$ALP),fromLast=TRUE)
# data_event$CREAT <- na.locf(na.locf(data_event$CREAT),fromLast=TRUE)
# data_event$NEUT <- na.locf(na.locf(data_event$NEUT),fromLast=TRUE)
# data_event$WBC <- na.locf(na.locf(data_event$WBC),fromLast=TRUE)


pam <- gam(ped_status~s(tend) +SLDb + WHOSTATN + SMKSTAT + s(ALP) + s(SLD) + s(CREAT), 
           select = T, method = 'REML', family=poisson(), offset=offset, data = data_event)
summary(pam)
AIC(pam)
draw(pam)
n_train   <- 150
train_idx <- sample(surv$USUBJID, n_train)
train_ped <- as_ped(
  data    = list(surv[surv$USUBJID %in% train_idx, ], biom[biom$USUBJID %in% train_idx, ]),
  formula = Surv(EVENT_TIME, EVENT) ~ . + concurrent(SLD, CREAT, ALP, tz_var = "TIME"),
  id      = "USUBJID")
pam1 <- pamm(
  formula = ped_status~s(tend) +SLDb + WHOSTATN + SMKSTAT + s(ALP) + s(SLD) + s(CREAT),
  data = train_ped)
pam2 <- pamm(
  formula = ped_status ~ s(tend) + SMKSTAT + WHOSTATN + SLDb,
  data = train_ped)

merged <- merge(surv, biom, by = c('USUBJID', 'SMKSTAT', 'WHOSTATN', 
                                   'EGFRMUTN', 'SLDb'), all =  F)
colnames(merged)[1] <- 'id'
merged$SLD <- na.locf(na.locf(merged$SLD),fromLast=TRUE)
merged$WBC <- na.locf(na.locf(merged$WBC),fromLast=TRUE)
merged$ALP <- na.locf(na.locf(merged$ALP),fromLast=TRUE)

merged <- merged[!is.na(merged$SLD), ]
merged <- merged[!is.na(merged$ALP), ]
merged <- merged[!is.na(merged$CREAT), ]
pec <- pec(
  list(pam1 = pam1, pam2 = pam2),
  Surv(EVENT_TIME, EVENT) ~ 1, # formula for IPCW
  data = merged[!(merged$id %in% train_idx), ],
  times = seq(.01, 400, by = 10),
  start = .01,
  exact = FALSE,
  reference = T,
  cens.model = 'aalen'
)
plot(pec)