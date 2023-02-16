packages = c('readxl', 'ggplot2', 'ggpubr',
             'stats', 'reshape2', 'ggrepel',
             'corrplot', 'tidyr', 'splitstackshape',
             'mgcv', 'writexl', 'gratia',
             'grid', 'viridis', 'pammtools',
             'pec')

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

options(warn=-1)
