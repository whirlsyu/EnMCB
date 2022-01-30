#' @import stats
#' @import survivalROC
#' @import glmnet
#' @import survivalsvm
#' @import ggplot2
#' @import survival
#' @import utils
#' @import BiocFileCache
#' @importFrom methods as
#' @importFrom mboost CoxPH
#' @importFrom mboost boost_control
#' @importFrom mboost glmboost
#' @importFrom mboost predict.glmboost



utils::globalVariables(c('FPF','TPF','aes','chr_id','coef','coord_cartesian','doParallel','data',
                         'element_blank','element_line','foreach','geom_abline','geom_roc','ggplot','rnorm',
                         'ggplot2','ggsave','groups','median','new','parallel','pheatmap','plot','demo_data',
                         'prediction','prognosticROC','qnorm','style_roc','survivalROC','theme','boot','boot.ci',
                         'theme_grey','theme_light'))

