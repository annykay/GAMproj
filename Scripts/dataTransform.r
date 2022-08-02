make_df <- function(data, mark) {
  base <- paste0(mark, 'b')
  df <- data[data$YTYPE_DESC == mark, c('USUBJID', 'SMKSTAT', 'TIME', 'WHOSTATN', 'YTYPE_DESC', 'DV',
                                          'EGFRMUTN', base, 'exclude')]
  df[, 'y'] <- df[,'DV'] - df[, base]
  df[, 'logy'] <- log(1 + df[, 'DV']) - log(df[, base])
  df[, 'logDV'] <- log(1 + df[, 'DV'])
  for (j in unlist(unique(df$USUBJID))) {
    curr_vec <- apply(df[df$USUBJID == j, 'y'], 2, as.numeric)
    df[df$USUBJID == j, 'norm'] <- (curr_vec - mean(curr_vec)) / sd(curr_vec)
  }
  return(df)
}

dfFactorize <- function(data, factors) {
  for (factor in factors){
    data[, factor] <- apply(data[,factor], 2, as.factor)
  }
  return(data)
}

dfOut <- function(data) {
  data$exclude <- 0

  data$exclude[data$USUBJID == 'EZHAD9EVZL'] <- 1
  data$exclude[data$USUBJID == '6J3DRPZZKH'] <- 2
  data$exclude[data$USUBJID == 'L5K6B1LQE5'] <- 2
  data$exclude[data$USUBJID == 'ZESXWGGXVU'] <- 2

  time_lim <- sort(unique(data$TIME[data$USUBJID == '2KL2236V79']),  decreasing = T)[3]
  data$exclude[data$USUBJID == '2KL2236V79' & data$TIME > time_lim] <- 2

  data$exclude[data$USUBJID == 'GHKU91P6DO' & data$TIME > 100] <- 2
  i1 <- table(data$USUBJID[data$YTYPE_DESC == 'SLD'])
  i1 <- rownames(i1)[i1 < 2]
  print(length(i1))
  for (i in i1) {
    min_time <- min(data$TIME[data$USUBJID == i])
    data$exclude[data$USUBJID == i & data$TIME > min_time] <- 1
  }

  #markers <- sort(unique(data$YTYPE_DESC))
  #for (i in markers) {
  #  base <- paste0(i, 'b')
  #  print(sum(!(data$USUBJID %in% i1)))
  #  data$exclude[data$USUBJID %in% i1 & data$YTYPE_DESC == i &
  #    (data$DV != data[, base] | data$TIME > 7) ] <- 1
  #}
  return(data)
}

dfNormalize <- function(data) {
  for (i in sort(unique(data$YTYPE_DESC))) {
    df_name <- paste0('df', i)
    assign(df_name, make_df(data, i))
  }
  data1 <- dfALT
  colnames(data1)[8] <- 'base'
  j <- 1
  for (i in sort(unique(data$YTYPE_DESC))) {
    df_name <- paste0('df', i)
    df <- get(df_name)
    colnames(df)[8] <- 'base'
    j <- j + 1
    if (i != 'ALT')
      data1 <- rbind(data1, df)
  }
  return(data1)
}
chooseTransform <- function(data, Ytrans, Xtranses, markers) {
  data$target <- 0
  data$target[data$YTYPE_DESC == 'SLD'] <- unlist(data[data$YTYPE_DESC == 'SLD', Ytrans])
  len <- length(markers)
  for (i in c(1:len)) {
    data$target[data$YTYPE_DESC == markers[i]] <- unlist(data[data$YTYPE_DESC == markers[i], Xtranses[i]])
  }
  return(data)
 }
