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
                 ncol = 3, nrow = 3)
  a <- annotate_figure(p,
                       left = textGrob("SLD", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                       bottom = textGrob("Time, days", gp = gpar(cex = 1), y = unit(0.95, "npc")))
  filename <- paste0('./Results/', 'Predictions_lin_', test, '.png')
  ggsave(filename, a, device = 'png')
  return(a)
}
test_sub <- c('L9BMJ3OI51', 'HILDJL520L', '2Q09QQJTEQ', '2F15P8EWSF', '2Q09QQJTEQ',
              'PB72VT2B1L', 'LNWOUSCAG0', 'Z4KRKCB0PK', '917ITST7MQ')
train_sub <- c('01D3N99Y8C', '6R763O4S7S', '6LMRNV5X98', 'CPBSGAPA2X',
              '21KMB4ZBVD', '2JJ4O1B58P', 'GW5IWV0QS2', '3QE0QBSPXF',
              '4WTY6J3F5W')

plot_several(train, train_sub, 'train')
plot_several(test, test_sub, 'Test')
line45 <- data.frame(x = c(0,max(train$SLD, na.rm = T)), y = c(0, max(train$SLD, na.rm = T)))
ggplot() + geom_point(data = train, aes(x = SLD_pred, y = SLD, color = USUBJID)) +
  theme_minimal(base_size = 12) + theme(legend.position="none") +
  xlab('Predicted') + ylab('Observed') +
  geom_line(data = line45, aes(x = x, y = y))
ggsave('Results/ForPresentation/train_lin_pred.png', device = 'png',
       dpi = 600, width = 6.595, height = 6.595, units = "cm")

line45 <- data.frame(x = c(0,max(test$SLD, na.rm = T)), y = c(0, max(test$SLD, na.rm = T)))
ggplot() + geom_point(data = test, aes(x = SLD_pred, y = SLD, color = USUBJID)) +
  theme_minimal(base_size = 12) + theme(legend.position="none") +
  xlab('Predicted') + ylab('Observed') +
  geom_line(data = line45, aes(x = x, y = y))
ggsave('Results/ForPresentation/test_lin_pred.png', device = 'png',
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

       fviz_pca_ind(res.pca,
                    col.ind = pca_data$q_SLD, # color by groups
                    palette = c("#00AFBB",  "#FC4E07"),
                    addEllipses = TRUE, # Concentration ellipses
                    ellipse.type = "confidence",
                    legend.title = "Groups",
                    repel = TRUE
                    )

library(ggbiplot)
pca_data <- data3[data3$TIME == 0, c('ALT', 'ALP', 'AST', 'CREAT', 'LDH', 'NLR', 'NEUT', 'WBC', 'SLD')]
 #pca_data <- data3[data3$TIME == 0, c('AST', 'CREAT', 'LDH', 'NLR',  'SLD')]
row.names(pca_data) <- data3$USUBJID[data3$TIME  == 0]
 pca_data <-na.omit(pca_data)
res.pca <- prcomp(pca_data[, c('ALT', 'ALP', 'AST', 'CREAT', 'LDH', 'NLR', 'NEUT', 'WBC')], scale = T, center = T)
 #res.pca <- prcomp(pca_data[, c('AST', 'CREAT', 'LDH', 'NLR')], scale = T, center = T)

q <- quantile(pca_data$SLD, c(0.25, 0.5, 0.75))
pca_data$q_SLD <- '0-25%'
pca_data$q_SLD[pca_data$SLD > q[[1]]] <- '25-50%'
pca_data$q_SLD[pca_data$SLD > q[[2]]] <- '50-75%'
pca_data$q_SLD[pca_data$SLD > q[[3]]] <- '75-100%'
pca_data$q_SLD <- as.factor(pca_data$q_SLD)
ggbiplot(res.pca, ellipse = T, scale = T, center = T, groups = pca_data$q_SLD, var.axes=FALSE) +theme_minimal()

## permutatios
data_perm <- data3
data_perm <- data_perm[, -6]

for (col in colnames(data_perm)) {
    print(col)
    if (col %in% c("ALP", "ALT", "AST", "CREAT", "LDH", "NEUT", "NLR", "WBC")) {
      data_perm[, col] <- sample(unlist(data3[,col]))
    }
}
shuffle <-function(i) {
  res1 <- stratified(for_strat, c('SMKSTAT', 'WHOSTATN', 'EGFRMUTN'), 0.7, bothSets= T)
  train <- data_perm[data_perm$USUBJID %in% res1$SAMP1$USUBJID, ]
  test <- data_perm[data_perm$USUBJID %in% res1$SAMP2$USUBJID, ]
  source('Scripts/model.r')
  res <- predict(ct1, train)
  train$SLD_pred <- res

  res <- predict(ct1, test)
  test$SLD_pred <- res

  r2_test <- 1 -  sum((test$SLD_pred - test$SLD) ** 2, na.rm = T) / sum((test$SLD - mean(test$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_train <- 1 -  sum((train$SLD_pred - train$SLD) ** 2, na.rm = T) / sum((train$SLD - mean(train$SLD, na.rm = T)) ** 2, na.rm = T)
  print(i)
  return(c(r2_test, r2_train))
}
r2_test_train <- lapply(c(1:10), shuffle)
r2_test_train_shuffle <- c()
r2_train_train_shuffle <- c()
r2_test_not_shuffle <- c()
r2_train_not_shuffle <- c()

for_strat <- data_perm[, c('USUBJID', 'SMKSTAT', 'WHOSTATN', 'EGFRMUTN')]
for_strat <- unique(for_strat)
for (i in c(1:10)) {
  res1 <- stratified(for_strat, c('SMKSTAT', 'WHOSTATN', 'EGFRMUTN'), 0.7, bothSets= T)
  train <- data_perm[data_perm$USUBJID %in% res1$SAMP1$USUBJID, ]
  test <- data3[data3$USUBJID %in% res1$SAMP2$USUBJID, ]
  source('Scripts/model.r')
  res <- predict(ct1, train)
  train$SLD_pred <- res

  res <- predict(ct1, test)
  test$SLD_pred <- res

  r2_test <- 1 -  sum((test$SLD_pred - test$SLD) ** 2, na.rm = T) / sum((test$SLD - mean(test$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_train <- 1 -  sum((train$SLD_pred - train$SLD) ** 2, na.rm = T) / sum((train$SLD - mean(train$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_test_train_shuffle <- c(r2_test_train_shuffle, r2_test)
  r2_train_train_shuffle <- c(r2_train_train_shuffle, r2_train)
  train <- data3[data3$USUBJID %in% res1$SAMP1$USUBJID, ]
  test <- data3[data3$USUBJID %in% res1$SAMP2$USUBJID, ]
  source('Scripts/model.r')
  res <- predict(ct1, train)
  train$SLD_pred <- res
  
  res <- predict(ct1, test)
  test$SLD_pred <- res
  
  r2_test <- 1 -  sum((test$SLD_pred - test$SLD) ** 2, na.rm = T) / sum((test$SLD - mean(test$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_train <- 1 -  sum((train$SLD_pred - train$SLD) ** 2, na.rm = T) / sum((train$SLD - mean(train$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_test_not_shuffle <- c(r2_test_not_shuffle, r2_test)
  r2_train_not_shuffle <- c(r2_train_not_shuffle, r2_train)
  
  print(i)
}
data_r2 <- data.frame(test = r2_test_not_shuffle[c(1:11)], 
                      train = r2_train_not_shuffle[c(1:11)],
                      test_shuffle = r2_test_train_shuffle[c(1:11)],
                      train_shuffle = r2_train_train_shuffle[c(1:11)])
data_meld <- melt(data_r2)
write.csv(data_r2, 'DerivedData/IPASSBLR2.csv')
my_xlab <- paste(levels(data_meld$variable),"\n(N=",table(data_meld$variable),")",sep="")

ggplot(data_meld, aes(x=variable, y=value, fill=variable)) +
      geom_boxplot(varwidth = TRUE, alpha=0.2)  +theme_minimal() +
      theme(legend.position="none") +
      scale_x_discrete(labels=my_xlab) + ylab('r2') + xlab("")



slope_mut <- rnorm(50, mean = 20, sd = 3)
slope_norm <- rnorm(50, mean = 2, sd = 3)
slope <- c(slope_mut, slope_norm)

mut <- rep(1, 50)
norm <- rep(0, 50)
mut <-c(mut, norm)
baseline <- rnorm(100, mean = 2, sd = 1)
time <- c(0:10)
data_mut <- slope %*% t(time) + baseline

data_mut <- slope %*% t(time)
data_mut <- data.frame(data_mut)
data_mut$USUBLID <- paste0("ID", c(1:100))
data_mut$mut <- mut
data_mut_long <- reshape(data_mut, varying  = list(1:10), direction  = 'long')
synthetic_data <- data_mut_long[, c('USUBLID', 'mut', 'time', 'X1')]
colnames(synthetic_data) = c('USUBJID', 'MUTATION', 'TIME', 'SLD')
for_strat <- synthetic_data[, c('USUBJID', 'MUTATION')]
for_strat <- unique(for_strat)

r2_test_train_shuffle <- c()
r2_train_train_shuffle <- c()

r2_test_not_shuffle <- c()
r2_train_not_shuffle <- c()

for (i in c(1:40)) {
  res1 <- stratified(for_strat, c('MUTATION'), 0.7, bothSets= T)
  train <- synthetic_data[synthetic_data$USUBJID %in% res1$SAMP1$USUBJID, ]
  test <- synthetic_data[synthetic_data$USUBJID %in% res1$SAMP2$USUBJID, ]
  #source('Scripts/model.r')
  ct_lin <- glmer(SLD ~ TIME|MUTATION, data = train)
  res <- predict(ct_lin, train)
  train$SLD_pred <- res

  res <- predict(ct_lin, test)
  test$SLD_pred <- res

  r2_test <-  1 -  sum((test$SLD_pred - test$SLD) ** 2, na.rm = T) / sum((test$SLD - mean(test$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_train <- 1 -  sum((train$SLD_pred - train$SLD) ** 2, na.rm = T) / sum((train$SLD - mean(train$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_test_not_shuffle <- c(r2_test_not_shuffle, r2_test)
  r2_train_not_shuffle <- c(r2_train_not_shuffle, r2_train)

  train[, 'SLD'] <- sample(unlist(train[, 'SLD']))
  ct_lin <- glmer(SLD ~ TIME|MUTATION, data = train)
  res <- predict(ct_lin, train)
  train$SLD_pred <- res

  res <- predict(ct_lin, test)
  test$SLD_pred <- res

  r2_test <- 1 -  sum((test$SLD_pred - test$SLD) ** 2, na.rm = T) / sum((test$SLD - mean(test$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_train <- 1 -  sum((train$SLD_pred - train$SLD) ** 2, na.rm = T) / sum((train$SLD - mean(train$SLD, na.rm = T)) ** 2, na.rm = T)
  r2_test_train_shuffle <- c(r2_test_train_shuffle, r2_test)
  r2_train_train_shuffle <- c(r2_train_train_shuffle, r2_train)

  print(i)
}
data_r2 <- data.frame(test = r2_test_not_shuffle, train = r2_train_not_shuffle,
                      test_shuffle = r2_test_train_shuffle, train_shuffle = r2_train_train_shuffle)
data_meld <- melt(data_r2)

my_xlab <- paste(levels(data_meld$variable),"\n(N=",table(data_meld$variable),")",sep="")

ggplot(data_meld, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2)  +theme_minimal() +
  theme(legend.position="none") +
  scale_x_discrete(labels=my_xlab) + ylab('r2') + xlab("")

