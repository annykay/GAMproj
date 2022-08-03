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
              s(TIME, k = 15)  +
              s(USUBJID, bs = 're', k = 400) +
              # te(ALT, TIME, bs = 'tp', k = 10) + #
              # te(ALP, TIME, bs = 'tp', k = 10) +
              ti(AST, TIME, bs = 'tp', k = 10) +
              te(CREAT, TIME, bs = 'tp', k = 10) +
              te(LDH, TIME, bs = 'tp', k = 10) +
              # ti(NLR, TIME, bs = 'tp', k = 10) +
              ti(NEUT, TIME, bs = 'tp', k = 10) +
              # te(WBC, TIME, bs = 'tp', k = 10) +
              # s(ALT, bs = 'tp', k = 10) +
              # s(ALP, bs = 'tp', k = 10) +
              # s(AST, bs = 'tp', k = 10) +
              # s(CREAT, bs = 'tp', k = 10) +
              # s(LDH, bs = 'tp', k = 10) +
              # s(NLR, bs = 'tp', k = 10) +
              s(NEUT, bs = 'tp', k = 10) +
              # s(WBC, bs = 'tp', k = 10) +
             SMKSTAT,
            data = train, select=TRUE,
            method="REML", family.mgcv = "gaulss")
# saveRDS(ct1, 'Models/BestNormWoID.rds')
# saveRDS(ct1, 'Models/curr.rds')
print(k.check(ct1))
print(anova(ct1))
print(AIC(ct1))