#' @import foreach
#' @import parallel
#' @import stats
#' @import survivalROC
#' @import doParallel
#' @import glmnet
#' @import survivalsvm
#' @import ggplot2
#' @import minfi
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import survival
#' @import utils


utils::globalVariables(c('%dopar%','FPF','TPF','aes','chr_id','coef','coord_cartesian','doParallel','data',
                         'element_blank','element_line','foreach','geom_abline','geom_roc','ggplot','rnorm',
                         'ggplot2','ggsave','groups','median','new','parallel','pheatmap','plot','demo_data',
                         'prediction','prognosticROC','qnorm','style_roc','survivalROC','theme',
                         'theme_grey','theme_light'))

