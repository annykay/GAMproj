plotMarkers <- function(data, title_end, filename_end, dirname, strat='') {
    descriptions <- sort(unique(data$YTYPE_DESC))
    for (i in c(1:9)) {
        baseline <- paste0(descriptions[i], 'b')
        if (strat != '')
            df <- data[data$YTYPE_DESC == descriptions[i], c('USUBJID', 'TIME', 'DV', baseline, strat)]
        else
            df <- data[data$YTYPE_DESC == descriptions[i], c('USUBJID', 'TIME', 'DV', baseline)]
        colnames(df)[4] <- 'base'
        df$curr <- df$DV - df$base

        title <- paste0(descriptions[i], title_end)
        filename <-  paste0(dirname, descriptions[i], filename_end, '.png')
        if (strat != '')
            ggplot(df, aes(TIME, curr, color = get(strat))) + scale_color_discrete(name=strat) + geom_point() +
                    geom_smooth() + xlab('Time, Days') + ylab(descriptions[i]) + ggtitle(title) + theme_minimal()
        else
            ggplot(df, aes(TIME, curr)) + geom_point() + geom_smooth() + xlab('Time, Days') +
                   ylab(descriptions[i]) + ggtitle(title) + theme_minimal()
        ggsave(filename)


        subjlist <- unique(df[df$TIME %in% c(0), 'USUBJID'])
        title <- paste0(descriptions[i], title_end, ' for patients came at first two days')
        filename <-  paste0(dirname, descriptions[i], filename_end, 'First', '.png')
        if (strat != '')
            ggplot(df[df$USUBJID %in% unlist(subjlist), ], aes(TIME, curr, color = get(strat))) + geom_point() +
                geom_smooth() + xlab('Time, Days') + ylab(descriptions[i]) + ggtitle(title) +
                scale_color_discrete(name=strat)
        else
            ggplot(df[df$USUBJID %in% unlist(subjlist), ], aes(TIME, curr)) + geom_point() +
                geom_smooth() + xlab('Time, Days') + ylab(descriptions[i]) + ggtitle(title)
        ggsave(filename)

        df$norm <- (df$curr - mean(df$curr)) / sd(df$curr)
        subjlist <- unique(df[df$TIME %in% c(1, 0), 'USUBJID'])
        title <- paste0(descriptions[i], title_end, ' normalized')
        filename <-  paste0(dirname, descriptions[i], filename_end,'Norm', '.png')
        if (strat != '')
            ggplot(df, aes(TIME, norm, color = get(strat))) + geom_point() + geom_smooth() +
                xlab('Time, Days') + ylab(descriptions[i]) + ggtitle(title) + ylim(-10, 10) +
                scale_color_discrete(name=strat)
        else
            ggplot(df, aes(TIME, norm)) + geom_point() + geom_smooth() + xlab('Time, Days') +
                ylab(descriptions[i]) + ggtitle(title) + ylim(-10, 10) + theme_minimal()
        ggsave(filename)

        df$normSu <- 0

        for (j in unlist(unique(df$USUBJID))) {
            curr_vec<- df[df$USUBJID == j, 'curr']
            df[df$USUBJID == j, 'normSu'] <- (curr_vec - mean(curr_vec)) / sd(curr_vec)
        }
        title <- paste0(descriptions[i], title_end, '\nStandart normalization for each patient')
        filename <-  paste0(dirname, descriptions[i], filename_end, 'NormBySubj', '.png')
        if (strat != '')
            ggplot(df, aes(TIME, normSu, color = get(strat))) + geom_point() + geom_smooth() +
                xlab('Time, Days') + ylab(descriptions[i]) + ggtitle(title) + ylim(-5, 5) +
                scale_color_discrete(name=strat) + theme_minimal()
        else
            ggplot(df, aes(TIME, normSu)) + geom_point() + geom_smooth() + xlab('Time, Days') +
                ylab(descriptions[i]) + ggtitle(title) + ylim(-5, 5)+ theme_minimal()

        ggsave(filename)
    }
}

plotBar <- function(data, strat, filename, title, ymax) {

    data_count = data.frame(table(unique(data[, c('USUBJID', strat)])[, strat]))
    colnames(data_count) <- c('names', 'heigth')
    png(filename=filename)


    my_bar <- barplot(data_count$heigth , border=F ,
                      las=2 ,
                      col=c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6) , rgb(0.3,0.9,0.4,0.6)) ,
                      ylim = c(0, ymax),
                      cex.main=2,
                      main=title)

    text(my_bar, data_count$heigth + 40, paste("n: ", data_count$heigth, sep="") ,cex=2)

    #Legende
    legend("topright", legend = data_count$names,
         col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6) , rgb(0.3,0.9,0.4,0.6) ) ,
         bty = "n", pch=20 , pt.cex = 2, cex = 2, horiz = FALSE, inset = c(0.05, 0.05))

        dev.off()
}

plot_spagetti <- function(data, marker, dir, maxy){
    row <- ggplot(data, aes(x = TIME, y = DV, group = USUBJID, color = USUBJID, show.legend = FALSE)) +
        geom_line() + geom_point() + theme(legend.position="none") + ylab(marker) +
            theme(axis.title=element_text(size=18,face="bold"))#+ ylim(0, maxy)

    basetitle <- paste0(marker, ' change frome baseline')
    base <- ggplot(data, aes(x = TIME, y = y, group = USUBJID, color = USUBJID, show.legend = FALSE)) +
        geom_line() + geom_point()  + theme(legend.position="none") +
        ylab(basetitle) + theme(axis.title=element_text(size=18,face="bold")) #+ ylim(0, maxy)

    p <- ggarrange(row, base,
              labels = c("A", "B"),
              ncol = 2, nrow = 1)

    filename <- paste0(dir, '/s', marker, '.png')
    ggsave(filename, device = 'png', width = 12, height = 8)
}
plot_diff_dist <- function(data){
  to_plot <- df_diff_dist(data)
  xmax <- mean(data$DV) + 3 * sd(data$DV)
  to_plot <- melt(to_plot, id.vars = 1)
  print(xmax)
  p <- ggplot(to_plot, aes(x = Count, y= value, color =variable)) + geom_point() + xlim(0, xmax)
  return(p)
}
plot_outlier <- function(data, subj) {
  markers <- sort(unique(data$YTYPE_DESC))

  new_data <- data[data$USUBJID == subj, ]
  plots <- lapply(markers, plot_marker, data = new_data)
  p <- ggarrange(plotlist=plots,
                 ncol = 3, nrow = 3,
                 common.legend = T)
  filename <- paste0('Results/', 'Out', subj, '.png')
  ggsave(filename, p, device = 'png')
  return(p)
}

plot_marker <- function(data, marker) {
  pl <- ggplot(data = data[data$YTYPE_DESC == marker, ],
               aes(x = TIME, y = DV, color = USUBJID)) +
    geom_point() + 
    geom_line() +
    xlab('TIME') + 
    ylab(marker) + 
    theme_minimal(base_size = 18) + 
    theme(legend.position="none")
    return(pl)
}