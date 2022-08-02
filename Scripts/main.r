source('./Scripts/dependences.r')
source('./Scripts/dataTransform.r')
source('./Scripts/plotFunc.r')

data <- read_excel('./SourceData/Erlotinib_dataset.xlsx')
data <- dfFactorize(data, c('SMKSTAT', 'WHOSTATN',  'EGFRMUTN'))
data <- dfOut(data)

data1 <- dfNormalize(data)
markers <- sort(unique(data1$YTYPE_DESC))[-8]
data1 <- chooseTransform(data1, 'norm', rep('norm', length(markers)), markers)

data2 <- data1[, c('USUBJID', 'SMKSTAT', 'TIME', 'WHOSTATN',  'EGFRMUTN', 'YTYPE_DESC', 'target', 'exclude')]
data2 <- spread(data2, key=YTYPE_DESC, value=target)

rows <- nrow(data2)

for (i in c(1:rows)) {
    if (!is.na(data2$SLD[i])) {
        subj <- data2$USUBJID[i]
        time <- data2$TIME[i]
        data2$SLD[data2$USUBJID == subj & abs(data2$TIME - time) < 4 & is.na(data2$SLD)] <- data2$SLD[i]
    }
}

print(sum(is.na(data2$AST) & is.na(data2$CREAT) & is.na(data2$NLR) & is.na(data2$LDH)))
data2$exclude[is.na(data2$AST) & is.na(data2$CREAT) & is.na(data2$NLR) & is.na(data2$LDH)] <- 3

set.seed(42)
for_strat <- data2[!(data2$exclude %in% c(1,3)), c('USUBJID', 'SMKSTAT', 'WHOSTATN', 'EGFRMUTN')]
for_strat <- unique(for_strat)
res <- stratified(for_strat, c('SMKSTAT', 'WHOSTATN', 'EGFRMUTN'), 0.7, bothSets= T)
train <- data2[!(data2$exclude %in% c(1,3)) & data2$USUBJID %in% res$SAMP1$USUBJID, ]
test <- data2[!(data2$exclude %in% c(1,3)) & data2$USUBJID %in% res$SAMP2$USUBJID, ]


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
