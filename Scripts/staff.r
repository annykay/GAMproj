res <- predict(ct1, train, se.fit = TRUE)
train$SLD_pred <- res$fit
train$SLD_pred_se <- res$se.fit

res <- predict(ct1, test, se.fit = TRUE)
test$SLD_pred <- res$fit
test$SLD_pred_se <- res$se.fit
ggplot(data = test[test$USUBJID %in% unique(test$USUBJID)[1],]) + 
  geom_ribbon(data = test[test$USUBJID %in% unique(test$USUBJID)[1],],
              aes(y = SLD_pred,
                  x = TIME,
                  ymin = SLD_pred - 2 * SLD_pred_se,
                  ymax = SLD_pred + 2 * SLD_pred_se,
                  group = USUBJID,
                  color = USUBJID),
              alpha = 0.25) +
  geom_point(aes(x = TIME,
                 y = SLD,
                 group = USUBJID,
                 color = USUBJID), 
             alpha = 0.7) + 
  geom_line(aes(x = TIME,
                 y = SLD,
                 group = USUBJID,
                 color = USUBJID), 
             alpha = 0.7) + 
  ylim(-2,2) +
  ylab('SLD') + 
  theme_minimal(base_size = 18)


i2 <- table(data1$USUBJID[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
i2 <- rownames(i2)[i2 > 4]
data1[data1$USUBJID == i2[2] & data1$YTYPE_DESC == 'CREAT',]
table(data1$TIME[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 71 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 88 & data1$YTYPE_DESC == 'CREAT'])
