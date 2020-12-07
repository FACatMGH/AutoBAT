# for extracting and formatting to use by Sarita's BATmeasures function
rm(list=ls())
library(readxl)
library(tidyverse)
library(caTools)
library(glue)
library(ggpubr)

workingdir <- "/Volumes/ShreffLabRemote/Vedanta BAT"
xlfiles <- list.files(workingdir, pattern = "*.xls$", full.names = TRUE, recursive = TRUE)
# manually removing the first 18 done prior to correcting low Ca issue
xlfiles <- xlfiles[!grepl("Suboptimal", xlfiles)]

exp_names <- map(strsplit(xlfiles, split = "/"), 5)
exp_names <- rep(exp_names, each = 23)

rawdata <- map_df(xlfiles, read_excel) %>%
  mutate(exp_name = as.character(exp_names))
# name and clean
names(rawdata) <- c("stim", "count", "CD203mfi", "CD63hi", "exp_name")
rawdata$CD63hi <- as.numeric(strsplit(rawdata$CD63hi, " %"))
rawdata <- filter(rawdata, stim != "Mean" & stim != "SD")
# use lookup table to clean stim conditions and allow easy edit if i have the Ag concentrations wrong

lut_agconc <- c(
  "1" = "NA",
  "2" = "NA",
  "3" = 2,
  "4" = 0.4,
  "5" = 0.08,
  "6" = 2,
  "7" = 0.4,
  "8" = 0.08,
  "9" = 0.016,
  "10" = 0.0032,
  "11" = 0.00064,
  "12" = 0.000128,
  "13" = 0.000026,
  "14" = 2,
  "15" = 0.4,
  "16" = 0.08,
  "17" = 0.016,
  "18" = 0.0032,
  "19" = 0.00064,
  "20" = 0.000128,
  "21" = 0.000026
)

rawdata$agconc <- NA
rawdata$stim <- as.numeric(map(strsplit(rawdata$stim, "\\.", perl = T), 1))
rawdata$agconc <- log10(as.numeric(lut_agconc[rawdata$stim]))

lut_stim <- c(
  "1" = "RPMI",
  "2" = "IL3",
  "3" = "a-IgE",
  "4" = "a-IgE",
  "5" = "a-IgE",
  "6" = "PN",
  "7" = "PN",
  "8" = "PN",
  "9" = "PN",
  "10" = "PN",
  "11" = "PN",
  "12" = "PN",
  "13" = "PN",
  "14" = "ARAH2",
  "15" = "ARAH2",
  "16" = "ARAH2",
  "17" = "ARAH2",
  "18" = "ARAH2",
  "19" = "ARAH2",
  "20" = "ARAH2",
  "21" = "ARAH2"
)
# now write back over stim with the corresponding lut_stim value
rawdata$stim <- lut_stim[rawdata$stim]
rawdata <- rawdata %>%
  mutate(patient = map(strsplit(rawdata$exp_name, "-"), 1)) %>%
  mutate(visit = map(strsplit(rawdata$exp_name, "-"), 2))


nested_data <- rawdata %>%
  group_by(exp_name, stim) %>%
  nest()


# integration for AUC using caTools -- note that AUC function Sarita used has been deprecated
#note that trapz requires data be ordered by x variable (for us, agconc)
mod_trapz <- function(df)
  trapz(df$agconc[order(df$agconc)], df$CD63hi[order(df$agconc)])

#functions for getting a measure of sensitivity: x value for peak response
# ended up using these in place of ED50 after extensively exploring the curve fitting models in drm package, which Sarita used, and finding variable and often poor fitting with any of them across all data in comparison to stat_smooth
# first function returns the x value, second function is for adding it along with other formatting to plots


cdmax.x <- function(gg.object) {
  gb <- ggplot_build(gg.object)
  max.x.index <- which(diff(sign(diff(gb$data[[1]]$y)))==-2)+1
  if(length(max.x.index)>1) {max.x <- gb$data[[1]]$x[which(gb$data[[1]]$y == max(gb$data[[1]]$y[max.x.index]))]} else {max.x <- gb$data[[1]]$x[max.x.index]}
  if(is_empty(max.x)) {max.x <- gb$data[[1]]$x[which(gb$data[[1]]$y == max(gb$data[[1]]$y))]}
  return(max.x)
}

final.plot <- function(gg.object) {
  gb <- ggplot_build(gg.object)
  max.x.index <- which(diff(sign(diff(gb$data[[1]]$y)))==-2)+1
  if(length(max.x.index)>1) {max.x <- gb$data[[1]]$x[which(gb$data[[1]]$y == max(gb$data[[1]]$y[max.x.index]))]} else {max.x <- gb$data[[1]]$x[max.x.index]}
  if(is_empty(max.x)) {max.x <- gb$data[[1]]$x[which(gb$data[[1]]$y == max(gb$data[[1]]$y))]}
  finalPlot <- gg.object + geom_vline(xintercept=max.x)
  return(finalPlot)
}


auc <- nested_data %>%
  dplyr::mutate(n   = data %>% map_dbl(~nrow(.)),
                sd  = data %>% map_dbl(~sd(.$CD63hi)),
                AUC = map(data, mod_trapz))


n_plots <- auc %>%
  filter(n > 1) %>%
  mutate(rawplot = map2(data, stim, ~ggplot(data = ., aes(x = agconc, y = CD63hi)) +
                      theme_light() +
                      geom_point() +
                      geom_smooth(span = 0.75) +
                      scale_y_continuous(limits = c(0, 100)) +
                      labs(title = glue({exp_name},": ", {stim}, "\n AUC: ", {round(as.numeric(AUC),2)}))
                      )
         )

edmax_with_plots <- auc %>%
  filter(n > 1 & AUC >10) %>%
#  mutate(drm_model = map(data, crs_drm)) %>%
  mutate(rawplot = map2(data, stim, ~ggplot(data = ., aes(x = agconc, y = CD63hi)) +
                          geom_smooth(span = 0.75)),
         max.x = map_dbl(rawplot, possibly(cdmax.x, NULL)),
         finalplot = map2(data, stim, ~ggplot(data = ., aes(x = agconc, y = CD63hi)) +
                          theme_light() +
                          geom_smooth(span = 0.75) +
                          geom_point() +
                          scale_y_continuous(limits = c(0, 100)) +
                          geom_vline(xintercept = max.x) +
                          labs(title = glue({as.character(exp_name)},": ", {as.character(stim)}, "\nAUC: ", {round(as.numeric(AUC),1)}," CDMax.x: ", {round(max.x,3)})))
  )


# export data

final_auc <- dplyr::select(edmax_with_plots, -data, -rawplot, -finalplot)
final_auc$AUC <- unlist(final_auc$AUC)
final_auc$max.x <- unlist(final_auc$max.x)
final_auc$stim <- unlist(final_auc$stim)
write.csv(as.data.frame(final_auc), "~/Dropbox (Partners HealthCare)/Projects/Vedanta-1453/analysis/BAT_analysis/summary_auc_data.csv")

rawdata$patient <- unlist(rawdata$patient)
rawdata$exp_name <- unlist(rawdata$exp_name)
rawdata$stim <- unlist(rawdata$stim)
rawdata$visit <- unlist(rawdata$visit)
write.csv(as.data.frame(rawdata), "~/Dropbox (Partners HealthCare)/Projects/Vedanta-1453/analysis/BAT_analysis/raw_data.csv")


# export plots
setwd("~/Dropbox (Partners HealthCare)/Projects/Vedanta-1453/analysis/BAT_analysis/")
nested_plots <- edmax_with_plots %>%
  dplyr::group_by(exp_name) %>%
  dplyr::select(exp_name, finalplot) %>%
  nest()



map2(paste0(nested_plots$exp_name, ".pdf", sep=""), map(nested_plots$data, ~cowplot::plot_grid(plotlist = .$finalplot)), ggsave)


