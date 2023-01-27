make_df <- function(data, mark) {
  base <- paste0(mark, 'b')
  df <- data[data$YTYPE_DESC == mark, c('USUBJID', 'SMKSTAT', 'TIME', 'WHOSTATN', 'YTYPE_DESC', 'DV',
                                          'EGFRMUTN', base, 'exclude')]
  df[, 'y'] <- df[,'DV'] - df[, base]
  df[, 'logy'] <- log(1 + df[, 'DV']) - log(1 + df[, base])
  df[, 'logDV'] <- log(1 + df[, 'DV'])
  for (j in unlist(unique(df$USUBJID))) {
    curr_vec <- as.numeric(df$y[df$USUBJID == j])
    df$norm[df$USUBJID == j] <- (curr_vec - mean(curr_vec)) / sd(curr_vec)
  }
  return(df)
}

dfFactorize <- function(data, factors) {
  for (factor in factors) {
    data[, factor] <- as.factor(data[, factor])
  }
  return(data)
}

dfOut <- function(data) {
  data$exclude <- 0
  
  data$exclude[data$USUBJID == 'ZUP61WOQNN'] <- 2
  data$exclude[data$USUBJID == '15N1NV7ZXJ'] <- 2
  data$exclude[data$USUBJID == 'DAMYDEUFTA'] <- 2
  data$exclude[data$USUBJID == 'KWNCNYVBAA'] <- 2
  data$exclude[data$USUBJID == '0WJRHUB7FB'] <- 2
  data$exclude[data$USUBJID == '16E1J5IX1I'] <- 2
  
  data$exclude[data$USUBJID == 'ZUP61WOQNN'] <- 2
  data$exclude[data$USUBJID == '15N1NV7ZXJ'] <- 2
  data$exclude[data$USUBJID == 'DAMYDEUFTA'] <- 2
  data$exclude[data$USUBJID == 'KWNCNYVBAA'] <- 2
  
  # data$exclude[data$USUBJID == 'EZHAD9EVZL'] <- 1
  # data$exclude[data$USUBJID == '6J3DRPZZKH'] <- 2
  # data$exclude[data$USUBJID == 'L5K6B1LQE5'] <- 2
  # data$exclude[data$USUBJID == 'ZESXWGGXVU'] <- 2

  time_lim <- sort(unique(data$TIME[data$USUBJID == '2KL2236V79']),  decreasing = T)[3]
  data$exclude[data$USUBJID == '2KL2236V79' & data$TIME > time_lim] <- 2

  data$exclude[data$USUBJID == 'GHKU91P6DO' & data$TIME > 100] <- 2
  
  i1 <- table(data$USUBJID[data$YTYPE_DESC == 'SLD'])
  i1 <- rownames(i1)[i1 < 2]
  data$exclude[data$USUBJID %in% i1] <- 1
  
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
 
get_times <- function(data, time_lim){
   times <- data.frame(times = colnames(t(sort(table(data$TIME), decreasing = T))))
   
   times$times <- apply(times, 1, as.numeric)
   timepoints <- c(times$times[1])
   for (i in times$times) { 
     last <- length(timepoints)
     if ((i - 3)> timepoints[last] && i < time_lim) {
       timepoints <- c(timepoints, i)
     }
   }
   differences <- c()
   len <- length(timepoints)
   for (i in c(1:(len-1)))
     differences <- c(differences, (timepoints[i + 1] + timepoints[i]) / 2)
   return(differences)
 }
 
baseline <- function(data, marker) {
  subjects <- unique(data$USUBJID)
  for (subj in subjects){
    base <- unique(data$base[data$USUBJID == subj & data$YTYPE_DESC == 'SLD'])
    data$base[data$USUBJID == subj] <- base
  }
  return(data)
}
 