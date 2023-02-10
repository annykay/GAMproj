library(haven)
library(dplyr)
data <- read_sas('./SourceData/NCT00457392/lab_safe.sas7bdat')
data <- data[!is.na(data$LABVALUE), ]
data$LABVALUE <- as.numeric(data$LABVALUE)
data$MULTIPLIER <- -1
data$VALUE <- -1
data$UNIT <- "N/A"
turn_to_standart_enzym <- function(data, biomarker) {
    checker <- c(data$LBTEST == biomarker)
    checker_temp <- c(data$LABUNITR == "IU/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 1
    checker_temp <- c(data$LABUNITR == "UKT/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60
    checker_temp <- c(data$LABUNITR == "umol/l")
    data[checker & checker_temp, "MULTIPLIER"] <- 60
    checker_temp <- c(data$LABUNITR == "U/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 1
    checker_temp <- c(data$LABUNITR == "UMOL/H/ML")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/SEC*L")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** (-3) * 60
    checker_temp <- c(data$LABUNITR == "MMOL/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** (-3)
    checker_temp <- c(data$LABUNITR == "MMOL/ML/H")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 6 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/L*HR")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "UMOL/HXML")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "UMOL/H.ML")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/HR X L")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MKMOL/HR*ML")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/HXL")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MCMOL/MLXHOUR")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/H.L")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MMOL/H-L")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "UMOL/ML/HR")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "MKMOL/HXML")
    data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 3 / 60
    checker_temp <- c(data$LABUNITR == "UMOL/C*L") #???????
    data[checker & checker_temp, "MULTIPLIER"] <- 1 / 60
    checker_temp <- c(data$LABUNITR == "NMOL/SEC*L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 / 10 ** 3
    checker_temp <- c(data$LABUNITR == "NMOL/SEK-L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 / 10 ** 3
    checker_temp <- c(data$LABUNITR == "NMOLS/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 / 10 ** 3
    checker_temp <- c(data$LABUNITR == "MMOLS/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 * 10 ** 3
    checker_temp <- c(data$LABUNITR == "NMOL/SEK.L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 / 10 ** 3
    checker_temp <- c(data$LABUNITR == "NMOL S/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60 / 10 ** 3
    checker_temp <- c(data$LABUNITR == "UMOL/SEC.L")
    data[checker & checker_temp, "MULTIPLIER"] <- 60
    checker_temp <- c(data$LABUNITR == "UKAT/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 1000 / 60
    checker_temp <- c(data$LABUNITR == "UKAL/L")
    data[checker & checker_temp, "MULTIPLIER"] <- 1000 / 60
    checker_temp <- c(data$LABUNITR == "M/E") #????
    data[checker & checker_temp, "MULTIPLIER"] <- 1
    checker_temp <- c(data$LABUNITR == "E/L") #????
    data[checker & checker_temp, "MULTIPLIER"] <- 1
    checker_temp <- c(data$LABUNITR == "UNIT/L") #????
    data[checker & checker_temp, "MULTIPLIER"] <- 1
    data$VALUE <- data$MULTIPLIER * data$LABVALUE
    data$UNIT <- "UMOL/L"
    return(data)
}
turn_to_standart_conc <- function(data, biomarker, M) {
  checker <- c(data$LBTEST == biomarker)
  checker_temp <- c(data$LABUNITR == "UMOL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "MG/DL")
  data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 4 / M
  checker_temp <- c(data$LABUNITR == "g/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 4 / M
  checker_temp <- c(data$LABUNITR == "mg/dL")
  data[checker & checker_temp, "MULTIPLIER"] <- 10 ** 4 / M
  checker_temp <- c(data$LABUNITR == "MKMOL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "MICROMOL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "MMOL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1 ** 3
  checker_temp <- c(data$LABUNITR == "MCMOL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "MKMIOL/L") # ????
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "UMOL/SL") # ????
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "NMOL/L") # ???? value in same range wo multiplier
  data[checker & checker_temp, "MULTIPLIER"] <- 1000
  checker_temp <- c(data$LABUNITR == "UIM/L") # ????
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  data$VALUE <- data$MULTIPLIER * data$LABVALUE
  data$UNIT <- "X10E3/UL"
  return(data)
}
turn_to_standart_cells <- function(data, biomarker) {
  checker <- c(data$LBTEST == biomarker)
  checker_temp <- c(data$LABUNITR == "X10E3/UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "10**3/mm**3")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "X10^3/UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "K/UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "THOUSAND/MM^3")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "X10^9/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "10*9/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == '10**3/mm**3" "')
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "%") ## As I need NLR, can be the same
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "10^3/UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "CELLS/UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 10 ** (-3)
  checker_temp <- c(data$LABUNITR == "X10^3") ## dont understand, mb exclud
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "X10^3 UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "#") ## again dont understand
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "K/CMM")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "BIL/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "X10 9TH/L")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "X10^3UL")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  checker_temp <- c(data$LABUNITR == "K/CUMM")
  data[checker & checker_temp, "MULTIPLIER"] <- 1
  data$VALUE <- data$MULTIPLIER * data$LABVALUE
  data$UNIT <- "X10E3/UL"
  return(data)
  }

same <- function(data, marker) {
  data$VALUE[data$LBTEST == marker] <- data$LABVALUE[data$LBTEST == marker]
  data$UNIT[data$LBTEST == marker] <- data$LABUNITR[data$LBTEST == marker]
  return(data)
}
#ALT, ALP no need to change, only one unit
data <- turn_to_standart_enzym(data, "ALANINE AMINOTRANSFERASE (ALT)")
data <- turn_to_standart_enzym(data, "ALKALINE PHOSPHATASE")
data <- data$LABUNITR[data$LBTEST == 'ALBUMIN' != data$LABUNITR == "%"]
data < same(data, 'ALBUMIN')
data <- same(data, 'BASOPHILS (%)')
#BASOPHILS (%) seems to have same Units
data <- turn_to_standart_cells(data, "BASOPHILS (ABSOLUTE)")
data <- turn_to_standart_conc(data,  'BILIRUBIN (TOTAL)', 584.66)
data <- turn_to_standart_conc(data,  'CALCIUM', 40.078)
data < same(data, 'CHLORIDE')
data < same(data, 'EOSINOPHILS (%)')
data <- turn_to_standart_cells(data, "EOSINOPHILS (ABSOLUTE)")
data < same(data, 'HEMOGLOBIN')
data < same(data, 'LYMPHOCYTES (%)')
data < same(data, 'MAGNESIUM')
data < same(data, 'MONOCYTES (%)')
data < same(data, 'MONOCYTES (ABSOLUTE)')
data < same(data, 'PHOSPHATE')
data < same(data, 'PLATELETS')
data < same(data, 'POTASSIUM')
data < same(data, 'PROTEIN (TOTAL)') #???
data < same(data, 'MONOCYTES (%)')


data <- turn_to_standart_enzym(data, "ASPARTATE AMINOTRANSFERASE (AST)")
data <- turn_to_standart_enzym(data, "LACTATE DEHYDROGENASE")
data <- turn_to_standart_conc(data, "CREATININE", 113)
data <- turn_to_standart_cells(data, "NEUTROPHILS (ABSOLUTE)")
data <- turn_to_standart_cells(data, "LYMPHOCYTES (ABSOLUTE)")
data <- data[data$LBTEST %in% c("ASPARTATE AMINOTRANSFERASE (AST)",
                                 "LACTATE DEHYDROGENASE", "CREATININE",
                                  "NEUTROPHILS (ABSOLUTE)", "LYMPHOCYTES (ABSOLUTE)"), ]
data_new <- data[, c('PID_A', 'COLLDAY', 'LBTEST', 'VALUE')]
colnames(data_new) <- c('USUBJID', 'TIME', 'YTYPE_DESC', 'VALUE')
#wide_data <- cast(data_new, PID_A + COLLDAY ~ LBTEST)
data_new$VALUE[c(6642, 6643, 6648, 6649)] <- mean(data_new$VALUE[c(6642, 6643, 6648, 6649)])
data_new <- distinct(data_new)
wide_data <- spread(data_new, key=YTYPE_DESC, value=VALUE)
 colnames(wide_data) <- c('USUBJID', 'TIME', 'AST', 'CREAT', 'LDH', 'LYMPH', 'NEUT')
wide_data$NLR = wide_data$NEUT / wide_data$LYMPH
SLD_data <- read_sas('./SourceData/NCT00457392/tmm_p.sas7bdat')
SLD_data <- SLD_data[SLD_data$TMMDIS == "LUNG", ]
SLD_agg <- SLD_data[SLD_data$LESTYPE == 'Target', ] %>%
           group_by(PID_A, EFDAY) %>%
           summarise(TMMDIA = sum(TMMDIA, na.rm = T), .groups = 'drop')
colnames(SLD_agg) <- c('USUBJID', 'TIME', 'SLD')
SLD_agg$TIME[SLD_agg$TIME < 0] <- 0
SLD_b <- unique(SLD_agg[SLD_agg$TIME == 0, c('USUBJID', 'SLD')])
colnames(SLD_b) <- c('USUBJID', 'SLDb')
wide_data$TIME[wide_data$TIME < 0] <- 0
final <- merge(wide_data, SLD_agg, by = c("USUBJID", "TIME"), all = T)
rows <- nrow(final)
for (i in c(1:rows)) {
    if (!is.na(final$SLD[i])) {
        subj <- final$USUBJID[i]
        time <- final$TIME[i]
        final$SLD[final$USUBJID == subj & abs(final$TIME - time) < 7 & is.na(final$SLD)] <- final$SLD[i]
    }
}
data <- merge(final, SLD_b, by = "USUBJID", all = T)

data <- data[is.finite(data$NLR),]

wide_data[wide_data$TIME == 0, ] %>% summarise(AST = mean(AST, na.rm = T),
                                     LDH = median(LDH, na.rm = T),
                                     CREAT = median(CREAT, na.rm = T),
                                     NLR = median(NLR, na.rm = T),
                                       .groups = 'drop')

df_new_study <- data
data_base <- df_new_study[df_new_study$TIME == 0, ] %>% group_by(USUBJID) %>%
            summarise(ASTb = mean(AST),
                      LDHb = mean(LDH),
                      CREATb = mean(CREAT),
                      NLRb = mean(NLR),
                      .groups = 'drop')

merged <- merge(df_new_study, data_base, by = 'USUBJID')  %>%
          group_by(USUBJID) %>%
          summarize(TIME = TIME,
                    AST = AST - ASTb,
                    LDH = LDH - LDHb,
                    CREAT = CREAT - CREATb,
                    NLR = NLR - NLRb,
                    SLD = SLD,
                    SLDb = SLDb,
                    .groups = 'drop')
final <- merged %>% group_by(USUBJID) %>%
              summarise(AST = (AST - mean(AST, na.rm = T)) / sd(AST, na.rm = T),
                        LDH = (LDH - mean(LDH, na.rm = T)) / sd(LDH, na.rm = T),
                        CREAT = (CREAT - mean(CREAT, na.rm = T)) / sd(CREAT, na.rm = T),
                        NLR = (NLR - mean(NLR, na.rm = T)) / sd(NLR, na.rm = T),
                        TIME = TIME,
                        SLD = SLD,
                        SLDb= SLDb,
                        .groups = 'drop')

res <- predict(ct1, final, se.fit = TRUE)
final$SLD_pred <- res$fit
line45 <- data.frame(x = c(0,max(final$SLD_pred, na.rm = T)), y = c(0, max(final$SLD_pred, na.rm = T)))
ggplot() + geom_point(data = final[final$TIME != 0,], aes(x = SLD_pred, y = SLD, color = USUBJID)) +
  theme_minimal(base_size = 12) + theme(legend.position="none") +
  xlab('Predicted') + ylab('Observed') +
  geom_line(data = line45, aes(x = x, y = y))
ggsave('Results/ForPresentation/new_ds_wo_base_pred.png', device = 'png',
       dpi = 600, width = 6.595, height = 6.595, units = "cm")
