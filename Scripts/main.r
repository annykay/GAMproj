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
studies <- unique(data_all$STUDY)

for (study in studies) {
  assign('data', data_all[data_all$STUDY == study, ])
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
  subjects <- unique(data$USUBJID[!is.na(data$SLD)])
  data2 <- data2[data2$USUBJID %in% subjects, ]
  data2 <- baseline(data2, 'SLD')
  colnames(data2)[which(colnames(data2) == 'base')] <- 'SLDb'
  
  data2 <- spread(data2, key=YTYPE_DESC, value=target)
  
  rows <- nrow(data2)
  
  
  for (i in c(1:rows)) {
    if (!is.na(data2$SLD[i])) {
      subj <- data2$USUBJID[i]
      time <- data2$TIME[i]
      data2$SLD[data2$USUBJID == subj & abs(data2$TIME - time) < 4 & is.na(data2$SLD)] <- data2$SLD[i]
    }
  }
  data2$exclude[is.na(data2$AST) & is.na(data2$CREAT) & is.na(data2$NLR) & is.na(data2$LDH)] <- 3
  data3 <- data2[!(data2$exclude %in% c(1, 2, 3)), ]
  data3$USUBJID <- as.factor(data3$USUBJID)
  
  filename <- paste0('DerivedData/wide', study, '.xlsx')
  write_xlsx(data2, filename)
  df_name <- paste0('data2', study)
  assign(df_name, data2)
  df_name <- paste0('data3', study)
  assign(df_name, data3)
}

for (study in studies) {
  df_name <- paste0('surv', study)
  assign(df_name, unique(data_all[data_all$STUDY == study, c('USUBJID', 'EVENT', 'EVENT_TIME', 'SLDb', 'WHOSTATN', 'EGFRMUTN', 'SMKSTAT')]))
}



# set.seed(42)
# 
# for_strat <- data3[, c('USUBJID', 'SMKSTAT', 'WHOSTATN', 'EGFRMUTN')]
# for_strat <- unique(for_strat)
# res1 <- stratified(for_strat, c('SMKSTAT', 'WHOSTATN', 'EGFRMUTN'), 0.7, bothSets= T)
# train <- data3[data3$USUBJID %in% res1$SAMP1$USUBJID, ]
# test <- data3[data3$USUBJID %in% res1$SAMP2$USUBJID, ]

#ct1 <- gam(SLD~
           #ti(ALT, TIME, bs = 'tp', k = 10) + #
           #ti(ALP, TIME, bs = 'tp', k = 10) +
           #s(TIME, k = 15)  +
           #ti(AST, TIME, bs = 'tp', k = 10) +
           #ti(CREAT, bs = 'tp', k = 10) +
           #ti(LDH, TIME, bs = 'tp', k = 10) +
           #ti(NLR, TIME, bs = 'tp', k = 10) +
           #ti(NEUT, TIME, bs = 'tp', k = 10) +
           #ti(WBC, TIME, bs = 'tp', k = 10) +
           #s(ALT, bs = 'tp', k = 10) +
           #s(ALP, bs = 'tp', k = 10) +
           #s(AST, bs = 'tp', k = 10) +
           #s(CREAT, bs = 'tp', k = 10) +
           #s(LDH, bs = 'tp', k = 10) +
           #s(NLR, bs = 'tp', k = 10) +
           #s(NEUT, bs = 'tp', k = 10) +
           #s(WBC, bs = 'tp', k = 10) +
           #SMKSTAT,
           #data = train, select=TRUE,
           #method="REML", family.mgcv = "gaulss")
#AIC(ct1)
