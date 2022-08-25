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
  p <- ggplot(data = data[data$USUBJID == subj,]) + 
        geom_ribbon(data = data[data$USUBJID == subj,],
                    aes(y = SLD_pred,
                        x = TIME,
                        ymin = SLD_pred - 2 * SLD_pred_se,
                        ymax = SLD_pred + 2 * SLD_pred_se,
                        group = USUBJID,
                        color = USUBJID, 
                        show.legend = FALSE),
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
        ylab('') + 
        ggtitle(subj) +
        xlab('') +
        theme_minimal(base_size = 14) +
        theme(legend.position="none", plot.title = element_text(size=12))
  return(p)
}

plot_several <- function(data, subjects, test) {
  
  plots <- lapply(subjects, plot_one, data = data)
  p <- ggarrange(plotlist=plots,
                 ncol = 3, nrow = 1)
  a <- annotate_figure(p, 
                       left = textGrob("SLD", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                       bottom = textGrob("Time, days", gp = gpar(cex = 1), y = unit(0.95, "npc")))
  filename <- paste0('./Results/', 'Predictions', test, '.png')
  ggsave(filename, a, device = 'png', 
         dpi = 600, width = 15.11, height = 6.595, units = "cm")
  return(a)
}
test_sub <- c('L9BMJ3OI51', 'HILDJL520L', '2Q09QQJTEQ')#2F15P8EWSF')#, '2Q09QQJTEQ')#, 
              # 'PB72VT2B1L', 'LNWOUSCAG0', 'Z4KRKCB0PK', '917ITST7MQ',
              # 'BK20V0KGC2')
train_sub <- c('01D3N99Y8C', '6R763O4S7S', '6LMRNV5X98')#, 'CPBSGAPA2X')#, 
              # '21KMB4ZBVD', '2JJ4O1B58P', 'GW5IWV0QS2', '3QE0QBSPXF',
              # '4WTY6J3F5W')

plot_several(train, train_sub, 'train')
plot_several(test, test_sub, 'Test')
line45 <- data.frame(x = c(0,max(train$SLD_pred, na.rm = T)), y = c(0, max(train$SLD, na.rm = T)))
ggplot() + geom_point(data = train, aes(x = SLD_pred, y = SLD, color = USUBJID)) + 
  theme_minimal(base_size = 12) + theme(legend.position="none") +
  xlab('Predicted') + ylab('Observed') + 
  geom_line(data = line45, aes(x = x, y = y))
ggsave('Results/ForPresentation/train_pred.png', device = 'png', 
       dpi = 600, width = 6.595, height = 6.595, units = "cm")

line45 <- data.frame(x = c(0,max(test$SLD_pred, na.rm = T)), y = c(0, max(test$SLD, na.rm = T)))
ggplot() + geom_point(data = test, aes(x = SLD_pred, y = SLD, color = USUBJID)) + 
  theme_minimal(base_size = 12) + theme(legend.position="none") +
  xlab('Predicted') + ylab('Observed') + 
  geom_line(data = line45, aes(x = x, y = y))
ggsave('Results/ForPresentation/test_pred.png', device = 'png', 
       dpi = 600, width = 6.595, height = 6.595, units = "cm")


i2 <- table(data1$USUBJID[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
i2 <- rownames(i2)[i2 > 4]
data1[data1$USUBJID == i2[2] & data1$YTYPE_DESC == 'CREAT',]
table(data1$TIME[data1$DV == 62 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 71 & data1$YTYPE_DESC == 'CREAT'])
table(data1$TIME[data1$DV == 88 & data1$YTYPE_DESC == 'CREAT'])

ggplot(data = data1[data1$YTYPE_DESC == 'ALP' & data1$USUBJID %in% c('2KL2236V79', 'AR4ZK3XPH2','Q38IUGVW9C', 'I7ESS1Y9IC'),], 
       aes(x = TIME, y = y, color = USUBJID)) + 
  theme_minimal(base_size = 18) +
  geom_point() + 
  geom_line() + 
  xlab('Time, Days') + 
  ylab('Change from baseline') +
  ggtitle('ALP')
ggsave('Results/ForPresentation/ALPdiff.png', device = 'png')

ggplot(data = data1[data1$YTYPE_DESC == 'ALP' & data1$USUBJID %in% c('2KL2236V79', 'AR4ZK3XPH2','Q38IUGVW9C', 'I7ESS1Y9IC'),], 
       aes(x = TIME, y = norm, color = USUBJID)) + 
  theme_minimal(base_size = 18) +
  geom_point() + 
  geom_line() + 
  xlab('Time, Days') + 
  ylab('Change from baseline') +
  ggtitle('ALP')
ggsave('Results/ForPresentation/ALPDiffNorm.png', device = 'png')


ggplot(data = data1[data1$exclude != 1, ], aes(x=YTYPE_DESC, y=DV, fill=YTYPE_DESC, color=YTYPE_DESC)) +
  geom_violin(width=2.1, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal(base_size = 12) +
  theme(legend.position="none") +
  ylim(0, 500) +
  coord_flip() + 
  xlab("Biomarkers") +
  ylab("Units")
ggsave('Results/ForPresentation/violinMarkers.png', device = 'png', 
       dpi = 600, width = 8.365, height = 6.595, units = "cm")

ggplot(data = data1[data1$exclude != 1, ], aes(x=YTYPE_DESC, y=norm, fill=YTYPE_DESC, color=YTYPE_DESC)) +
  geom_violin(width=2.1, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_minimal(base_size = 12) +
  theme(legend.position="none") +
  ylim(0, 5) +
  coord_flip() + 
  xlab("Biomarkers") +
  ylab("Normalized Units")
ggsave('Results/ForPresentation/violinMarkersNorm.png', device = 'png', 
       dpi = 600, width = 8.365, height = 6.595, units = "cm")

plot_outlier(data1, 'EZHAD9EVZL', 30)

corrplot.mixed(cor_matrix, lower.col = COL1("Blues"))
ggsave('Results/ForPresentation/corrplot.png', device = 'png', 
       dpi = 600, width = 8.365, height = 6.595, units = "cm")

draw(ct1, select = 2)
ggsave('Results/cov1.png', 
       dpi = 600, width = 7.24, height = 6.595, units = "cm")

draw(ct1, select = 3)
ggsave('Results/cov2.png', 
       dpi = 600, width = 7.24, height = 6.595, units = "cm")

draw(ct1, select = 4)
ggsave('Results/cov3.png', 
       dpi = 600, width = 7.24, height = 6.595, units = "cm")

draw(ct1, select = 5)
ggsave('Results/cov4.png', 
       dpi = 600, width = 7.24, height = 6.595, units = "cm")