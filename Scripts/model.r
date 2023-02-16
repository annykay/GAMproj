# Best for norm
#ct1 <- gam(SLD~
#              s(TIME, k = 15)  +
#              s(USUBJID, bs = 're', k = 400) +
#              # te(ALT, TIME, bs = 'tp', k = 10) + #
#              # te(ALP, TIME, bs = 'tp', k = 10) +
#              # ti(AST, TIME, bs = 'tp', k = 10) +
#              # ti(CREAT, TIME, bs = 'tp', k = 10) +
#              # ti(LDH, TIME, bs = 'tp', k = 10) +
#              ti(NLR, TIME, bs = 'tp', k = 10) +
#              # ti(NEUT, TIME, bs = 'tp', k = 10) +
#              # te(WBC, TIME, bs = 'tp', k = 10) +
#              # s(ALT, bs = 'tp', k = 10) +
#              # s(ALP, bs = 'tp', k = 10) +
#              s(AST, bs = 'tp', k = 10) +
#              s(CREAT, bs = 'tp', k = 10) +
#              s(LDH, bs = 'tp', k = 10) +
#              # s(NLR, bs = 'tp', k = 10) +
#              s(NEUT, bs = 'tp', k = 10) +
#              # s(WBC, bs = 'tp', k = 10) +
#             SMKSTAT,
#            data = train, select=TRUE,
#            method="REML", family.mgcv = "gaulss")

ct1 <- gam(SLD~
              s(TIME)  +
              # NEUT +
              SLDb + 
              # s(USUBJID, bs = 're', k = 400) +
              te(ALT, TIME, bs = 'tp', k = 10) + #IPASS #INTEREST
              te(ALP, TIME, bs = 'tp', k = 10) +#IPASS
              te(AST, TIME, bs = 'tp', k = 10) + #IPASS #ZODIAC
              te(CREAT, TIME, bs = 'tp', k = 10) + #ZODIAC
              te(LDH, TIME, bs = 'tp', k = 10) +
              te(NLR, TIME, bs = 'tp', k = 10) +
              te(NEUT, TIME, bs = 'tp', k = 10) +#IPASS
              te(WBC, TIME, bs = 'tp', k = 10)+
              # s(ALT, bs = 'tp', k = 10) +
              # s(ALP, bs = 'tp', k = 10) +
              # s(AST, bs = 'tp', k = 10) +
              # s(CREAT, bs = 'tp', k = 10) +
              # s(LDH, bs = 'tp', k = 10) +
              # s(NLR, bs = 'tp', k = 10) +
              # s(NEUT, bs = 'tp', k = 10) +
              # s(WBC, bs = 'tp', k = 10) +
            SMKSTAT + EGFRMUTN +WHOSTATN,
            data = train, select=TRUE,
            method="REML", family.mgcv = "gaulss")
# saveRDS(ct1, 'Models/BestNormWoID.rds')
saveRDS(pam, 'Models/IPASS_PAMM.rds')
#print(k.check(ct1))
#print(anova(ct1))
#print(AIC(ct1))
# gamsld <- gam(y~s(TIME,  by = WHOSTATN, bs = 'ad', k = 10, id = 1) +
#              s(TIME, by = EGFRMUTN, bs = 'ad', k = 10, id = 1) +
#              s(TIME, by = SMKSTAT, bs = 'ad', k = 10, id = 1) +
#              SMKSTAT + EGFRMUTN + WHOSTATN, data = data1[datq1$exclude != 1 & data1$YTYPE_DESC == 'SLD', ], select=TRUE,
#            method="REML", family.mgcv = "gauss")
# gam.check(ct1)
# plot(ct1, residuals=TRUE, pch=19, ylim = c(-5,5), shade=TRUE, seWithMean=TRUE, scale=0)
#
#   ct2 <- gam(SLD ~ s(TIME, bs = 'ad', k = 10) +
#                ti(ALP, TIME, bs = 'tp', k = 10) +
#                te(AST, TIME, bs = 'cs', k = 10) +
#                s(LDH, k = 10) +
#                s(NEUT, k = 10) +
#               SMKSTAT,
#             data = train,
#             select = TRUE,
#             method = "REML",
#             family.mgcv = 'gaulss')
#
#
