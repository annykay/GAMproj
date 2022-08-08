res <- predict(ct1, train, se.fit = TRUE)
train$SLD_pred <- res$fit
train$SLD_pred_se <- res$se.fit

res <- predict(ct1, test, se.fit = TRUE)
test$SLD_pred <- res$fit
test$SLD_pred_se <- res$se.fit

# ggplot(data = test[test$USUBJID %in% unique(test$USUBJID),]) + 
#   geom_ribbon(data = test[test$USUBJID %in% unique(test$USUBJID)[1],],
#               aes(y = SLD_pred,
#                   x = TIME,
#                   ymin = SLD_pred - 2 * SLD_pred_se,
#                   ymax = SLD_pred + 2 * SLD_pred_se,
#                   group = USUBJID,
#                   color = USUBJID),
#               alpha = 0.25) +
#   geom_point(aes(x = TIME,
#                  y = SLD,
#                  group = USUBJID,
#                  color = USUBJID), 
#              alpha = 0.7) + 
#   geom_line(aes(x = TIME,
#                  y = SLD_pred,
#                  group = USUBJID,
#                  color = USUBJID), 
#              alpha = 0.7) + 
#   ylab('SLD') + 
#   theme_minimal(base_size = 18)


plot_one <- function(data, subj) {
  p <- ggplot(data = test[test$USUBJID == subj,]) + 
        geom_ribbon(data = test[test$USUBJID == subj,],
                    aes(y = SLD_pred,
                        x = TIME,
                        ymin = SLD_pred - 2 * SLD_pred_se,
                        ymax = SLD_pred + 2 * SLD_pred_se,
                        group = USUBJID,
                        color = USUBJID,
                    alpha = 0.25) +
        geom_point(aes(x = TIME,
                       y = SLD,
                       group = USUBJID,
                       color = USUBJID,
                       show.legend = FALSE), 
                   alpha = 0.7) + 
        geom_line(aes(x = TIME,
                      y = SLD_pred,
                      group = USUBJID,
                      color = USUBJID,
                      show.legend = FALSE), 
                  alpha = 0.7) + 
        ylab('SLD') + 
        xlab(subj) +
        theme_minimal(base_size = 18)
        # theme(legend.position="none") 
  return(p)
}

plot_several <- function(data, subjects, test) {

  plots <- lapply(subjects, plot_one, data = data)
  p <- ggarrange(plotlist=plots,
                 ncol = 3, nrow = 3)
  filename <- paste0('Results/', 'Predictions', test, '.png')
  ggsave(filename, p, device = 'png')
  return(p)
}
test_sub <- c('013S6S0LY8', '0F90KK8AGP', '0UTEJBN9VD', '0YLHHEVP8A', 
              '0ZY7F2OMR0', '1E20WTA6X7', '1JLPO0XOJW', '1OV6MXIPM2', 
              '1SRF90PUPP')
i2 <- table(data1$USUBJID[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
i2 <- rownames(i2)[i2 > 4]
data1[data1$USUBJID == i2[2] & data1$YTYPE_DESC == 'CREAT',]
table(data1$TIME[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 71 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 88 & data1$YTYPE_DESC == 'CREAT'])
