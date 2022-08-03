basic_biomarker_statistics <- function(data, markers) {
  npoints <- unlist(lapply(markers, count_points, data = data))
  npatients <- unlist(lapply(markers, count_patients, data = data))
  means <- unlist(lapply(markers, calc_mean, data = data))
  
  medians <- unlist(lapply(markers, calc_median, data = data))
  sds <- unlist(lapply(markers, calc_sd, data = data))
  mads <- unlist(lapply(markers, calc_mad, data = data))
  
  iqrs <- unlist(lapply(markers, calc_iqr, data = data))
  q0 <- unlist(lapply(markers, calc_quart, p = 0, data = data))
  q25 <- unlist(lapply(markers, calc_quart, p = 0.25, data = data))
  
  q50 <- unlist(lapply(markers, calc_quart, p = 0.5, data = data))
  q75 <- unlist(lapply(markers, calc_quart, p = 0.75, data = data))
  q100 <- unlist(lapply(markers, calc_quart, p = 1, data = data))
  
  lower <- means - sds
  upper <- means + sds
  
  basic_statistics <- data.frame(Marker = markers, 
                                 Count = npoints,
                                 Patients = npatients,
                                 Mean = means, 
                                 Sd = sds, Lower = lower, Upper = upper, 
                                 Median = medians, Mad = mads, IQR = iqrs,
                                 q0 = q0, q25 = q25, q50 = q50, 
                                 q75 = q75, q100 = q100)
  colnames(basic_statistics)[c(6, 7, 11, 12, 13, 14, 15)] <- 
    c("CI Lower", "CI Upper", "0%", "25%", "50%", "75%", "100%")
  return(basic_statistics)
}
count_points <- function(marker, data) {
  return(nrow(data[data$YTYPE_DESC == marker,]))
}
count_patients <- function(marker, data) {
  return(length(unique(data$USUBJID[data$YTYPE_DESC == marker])))
}
calc_mean <- function(marker, data) {
  return(round(mean(data$DV[data$YTYPE_DESC == marker], na.rm = T), 2))
}
calc_median <- function(marker, data) {
  return(round(median(data$DV[data$YTYPE_DESC == marker], na.rm = T), 2))
}
calc_sd <- function(marker, data) {
  return(round(sd(data$DV[data$YTYPE_DESC == marker], na.rm = T), 2))
}
calc_mad <- function(marker, data) {
  return(round(mad(data$DV[data$YTYPE_DESC == marker]), 2))
}
calc_iqr <- function(marker, data) {
  return(round(IQR(data$DV[data$YTYPE_DESC == marker], na.rm = T), 2))
}
calc_quart <- function(marker, data, p) {
  return(round(quantile(data$DV[data$YTYPE_DESC == marker], p, na.rm = T), 2))
}

factor_stat <- function(factors, data) {
  zeroes <- unlist(lapply(factors, count_amount, val = 0, data = data))
  ones <- unlist(lapply(factors, count_amount, val = 1, data = data))
  twoes <- unlist(lapply(factors, count_amount, val = 2, data = data))
  res <- data.frame(Factor = factors, zeroes = zeroes, one = ones, two = twoes)
  colnames(res) <- c("Factor", "0", "1", "2")
  return(res)
}

count_amount <- function(factor, val, data){
  return(length(unique(data$USUBJID[data[, factor] == val])))
}