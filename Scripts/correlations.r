calc_corr <- function(data, first, second, subjid) {
    df1 = data[data$USUBJID == subjid, ]
    df1$curr <- 0
    first_base <- paste0(first, 'b')
    second_base <- paste0(second, 'b')
    df1[df1$YTYPE_DESC == first, 'curr'] <- df1[df1$YTYPE_DESC == first, 'DV'] - df1[df1$YTYPE_DESC == first,
                                                                                     first_base]
    df1[df1$YTYPE_DESC == second, 'curr'] <- df1[df1$YTYPE_DESC == second, 'DV'] - df1[df1$YTYPE_DESC == second,
                                                                                     second_base]
    df1 <- df1[df1$YTYPE_DESC %in% c(first, second), c('curr', 'TIME', 'YTYPE_DESC')]
    spread(df1, key = TIME, value = curr)
    df2 = t(spread(df1, key = TIME, value = curr))
    df2 <- data.frame(df2)
    if (ncol(df2) == 2){

        colnames(df2) <- c(first, second)
        df2 <- df2[-1, ]
        df2 <- df2[!is.na(df2[, first]), ]
        df2 <- df2[!is.na(df2[, second]), ]
        df2 <- apply(df2, 2, as.numeric)

        if (is.matrix(df2) && nrow(df2) > 2)
            return (cor(df2[, first], df2[, second], method = 'kendall'))
        else
            return (0)
    } else {
        return(0)
    }
}

calc_cor_matrix <- function(data, subjid) {
    cor_matrix = matrix(0, nrow = 8, ncol = 8)
    markers = sort(unique(data$YTYPE_DESC))
    markers = markers[-8]
    for (i in c(1:7)){
        for (j in c((i+1):8)){
            correlation <-  calc_corr(data, markers[i], markers[j], subjid)
            if (!is.na(correlation))
                cor_matrix[i, j] <- correlation[1]
                cor_matrix[j, i] <- correlation[1]
        }
    }
    return (cor_matrix)
}
calc_av <- function(data) {
    av_cor_matrix <- matrix(0, nrow = 8, ncol = 8)
    usubjlist <- unlist(unique(data$USUBJID))
    temp_matrix <- matrix(0, nrow = 8, ncol = 8)
    for (i in usubjlist) {
        temp_matrix <- calc_cor_matrix(data, i)
        av_cor_matrix <- av_cor_matrix + temp_matrix
    }

    av_cor_matrix <- av_cor_matrix /length(usubjlist)
    for (i in (1:8))
        av_cor_matrix[i, i] <- 1
    return (av_cor_matrix)
}
