# Data-driven analytical algorithm for basophil activation data: Derivation of basophil activation measures
# Written by: Sarita Patil, MD. Function published 6/8/2017

#' @import rio
#' @import magrittr
#' @import xtable
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import DescTools
#' @import drc


# Modeling functions
mean_drm <- function(df){
  drm(
    CD63hi~agconc, data=df,  robust='mean',
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
    control = drmc(maxIt = 1000)
  )
}

median_drm <- function(df){
  drm(
    CD63hi~agconc, data=df,  robust='median',
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
    control = drmc(maxIt = 1000)
  )
}

lms_drm  <- function(df){
  drm(
    CD63hi~agconc, data=df,  robust='lms',
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
    control = drmc(maxIt = 1000)
  )
}

cons_drm  <- function(df){
  drm(
    CD63hi~agconc, data=df,  robust='mean',
    fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
    lowerl=c(NA, 0, 0, NA), upperl=c(NA, 25, 100, NA),
    control = drmc(maxIt = 1000)
  )
}

#' BATmeasures: Measures of basophil activation testing
#'
#' @return This function analyzes basophil activation dose-response curve to produce measures of basophil activation
#'
#' @author Sarita Patil, \email{sarita.patil@@mgh.harvard.edu}
#' @param data A dataframe that contains several variables:
#'        groupvariables which is a grouping variable that identifies all the conditions in a dose-response curve, agconc which is the antigen concentration used in basophil stimulation, and CD63hi which is the percentage of CD63hi basophils from flow analysis in each stimulation condition)
#' @param groupvariables The variable name of the grouping variable from the data
#' @param output.file The name of the output file
#'
#' @return A csv file containing calculated AUC and ED50 values based on input data
#'
#' @keywords models, nonlinear
#' @examples
#' BATmeasures(data, patientid, output.file="BAT_AUC_ED50.csv")
#' @export
BATmeasures <- function(data, groupvariables, output.file){
#  data$agconc <- as.numeric(data$agconc)
#  data$CD63hi <- as.numeric(data$CD63hi)

  auc <- data %>%
    group_by(groupvariables) %>%
    nest(agconc,CD63hi) %>%
    dplyr::mutate(n   = data %>% map_dbl(~nrow(.)),
           sd  = data %>% map_dbl(~sd(.$CD63hi)),
           AUC = data %>% map_dbl(~AUC(log10(.$agconc), .$CD63hi, method = "trapezoid")))

  intermed <- auc %>%
    filter(sd>0) %>%
    dplyr::mutate( mean_drm   = map(data, possibly(mean_drm,NULL)),
            median_drm = map(data, possibly(median_drm,NULL)),
            lms_drm    = map(data, possibly(lms_drm,NULL)),
            cons_drm   = map(data, possibly(cons_drm,NULL)))

  ed50 <- intermed %>%
    gather(model,value,mean_drm,median_drm,lms_drm,cons_drm) %>%
    dplyr::select(groupvariables, data,n,sd,AUC,model,value) %>%
    mutate(is_null = map_lgl(value, is.null)) %>%
    filter(!is_null) %>%
    dplyr::select(-is_null) %>%
    mutate(MX   = map_dbl(data,~ max(.$agconc)),
           MN   = map_dbl(data,~ min(.$agconc)),
           ED50 = map_dbl(value,~ED(.,50,display=F)[1]),
           SE   = map_dbl(value,~ED(.,50,display=F)[2]),
           UL   = map_dbl(value,~coefficients(.)[3]),
           LL   = map_dbl(value,~coefficients(.)[2]),
           SL   = map_dbl(value,~coefficients(.)[1]),
           AIC  = map(value,~mselect(.,icfct=BIC)) %>% map_dbl("IC")) %>%
    separate(model,c("model","ex")) %>% dplyr::select(-ex)

  ed50 <- ed50 %>%
    filter(between(UL,0,100), between(LL,0,100), SL < 0, ED50<MX, ED50>MN) %>%
    arrange(groupvariables) %>%
    dplyr::select(groupvariables, ED50,UL,LL,SL) %>%
    group_by(groupvariables) %>%
    summarize_each(funs(median(.,na.rm=T)), ED50, UL, LL, SL)

  all <- left_join(auc %>% dplyr::select(-data),ed50)

  write.csv(all, output.file)
}




