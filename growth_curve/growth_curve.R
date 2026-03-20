library(ggplot2)
library(tidyverse)

#############################################################################
### file reading
setwd("")
lacidophilus_48h_od600 <- read.csv("lacidophilus_growth_curve_48h.csv", header = T)
ecoli_24h_od600 <- read.csv("ecoli_37d_1nacl_LB_od600.csv", header = T)

ecoli_24h_od600_summary <- ecoli_24h_od600 %>%
  group_by(time) %>%
  summarise(mean_od600 = mean(od_600),
            mean_ln600 = mean(ln600),
            mean_log10_cfu=mean(log10_cfu),
            mean_cfu=mean(cfu_ml),
            sd_od600 = sd(od_600),
            sd_ln600 = sd(ln600), 
            sd_log10_cfu=sd(log10_cfu), 
            sd_cfu=sd(cfu_ml), .groups = "drop")

lacidophilus_48h_od600_summary <- lacidophilus_48h_od600 %>%
  group_by(time) %>%
  summarise(mean_od600 = mean(od_600),
            mean_ln600 = mean(ln600),
            mean_log10_cfu=mean(log10_cfu),
            mean_cfu=mean(cfu_ml),
            sd_od600 = sd(od_600),
            sd_ln600 = sd(ln600), 
            sd_log10_cfu=sd(log10_cfu), 
            sd_cfu=sd(cfu_ml), .groups = "drop")


#############################################################################
### functions of Gompertz model for data imputation
fit_gompertz <- function(data_table){
  pre_data_form <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("time", "pre_od600", "model", "high_time"))))
  data=data_table
  i <- which.max(diff(data$mean_od600)) ## find out the highest growth rate in the data frame
  
  ## calculate the values for the formula
  a=max(data$mean_od600)
  mu=max(diff(data$mean_od600))/(data[i,"time"]-data[i-1, "time"])
  lambda=(mu*data[i, "time"]-data[i, "mean_od600"])/mu
  ## calculate the values for the formula
  starting_values <- list(a=a, ## the asymptote
                          mu=mu, ## the highest growth rate
                          lambda=lambda) ## the lag time
  print(starting_values$lambda)
  print(c)
  formula_gompertz <- "mean_od600 ~ a*exp(-exp(mu*exp(1)/a*(lambda-time)+1))"
  fitted_gompertz_model <- nlsLM(formula_gompertz, data, starting_values)
  fitted_a <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["a"]]
  fitted_mu_time <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["mu"]][["time"]]
  fitted_lambada <- environment(fitted_gompertz_model[["m"]][["predict"]])[["env"]][["lambda"]][["time"]]
  
  pre_time=seq(0,48,by=0.1)
  y <- fitted_a*exp(-exp(fitted_mu_time*exp(1)/fitted_a*(fitted_lambada-pre_time)+1))
  mod <- "Gompertz"
  inter_dataframe <- data.frame(time=pre_time, pre_od600=y)
  high_num <- which.max(diff(inter_dataframe$pre_od600))
  high_time <- inter_dataframe[high_num, "time"]
  inter_dataframe$model <- mod
  inter_dataframe$high_time <- high_time
  pre_data_form <- pre_data_form %>%
    rbind(inter_dataframe)
  return(pre_data_form)
}

#############################################################################
### execute functions for data imputation
lacidophilus_prediction_gompertz <- fit_gompertz(lacidophilus_48h_od600_summary)
ecoli_prediction_gompertz <- fit_gompertz(ecoli_24h_od600_summary)

#############################################################################
### visualize the observed and imputated growth curve
ggplot() +
  geom_line(data=lacidophilus_prediction_gompertz, aes(x=time, y=pre_od600, color="pre_od600"), size=1, alpha=0.9) +
  geom_line(data=lacidophilus_48h_od600_summary, aes(x=time, y=mean_od600, color="mean_od600"), size=1, alpha=0.9) +
  scale_color_manual(name="Data type",
                     values=c("pre_od600"="#fd990d", "mean_od600"="#644b77"),
                     labels=c("Observed OD600", "Predicted OD600")) +
  theme_bw() +
  labs(x="Time (h)", y="OD600") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 22),
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank())


ggplot() +
  geom_line(data=ecoli_prediction_gompertz, aes(x=time, y=pre_od600, color="pre_od600"), size=1, alpha=0.9) +
  geom_line(data=ecoli_24h_od600_summary, aes(x=time, y=mean_od600, color="mean_od600"), size=1, alpha=0.9) +
  scale_color_manual(name="Data type",
                     values=c("pre_od600"="#fd990d", "mean_od600"="#644b77"),
                     labels=c("Observed OD600", "Predicted OD600")) +
  theme_bw() +
  labs(x="Time (h)", y="OD600") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 22),
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank())