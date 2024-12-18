library(tidyverse)

rm(list = ls())
source("project/proj_utils.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# load data and light tidying #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
wd <- read_csv("project/wd_gundeaths.csv") 
lower_names <- str_to_lower(names(wd))
wd <- wd |> `colnames<-`(lower_names)
names(wd)[str_detect(names(wd), "percent.+")] <- 
  str_replace(names(wd)[str_detect(names(wd), "percent.+")], "percent[s]{0,1}","pct")
names(wd) <- str_replace_all(names(wd),"\\.","_")
wd <- wd |> filter(state == "Texas") |> filter(year >= 1999)


# declare outcome variable
# ------------------------
total_firearm_deaths <- ts(1e4 * wd$total_firearm_deaths/wd$population)

# declare covariates
# ------------------
xcovariates <- cbind(wd$permit, wd$poverty_rate)


# CV each ARIMA model
# -------------------
n <- 0
orders <- vector(mode = "list") #, length = 32)
ARs <- c(0,1,2,3,4)
MAs <- c(0,1,2,3,4)
Is  <- c(0)
for(i in Is) {
  for(j in ARs) {
    for(k in MAs){
        n <- n + 1
        orders[[n]] <- c(j,i,k)
    }
  }
}
names(orders) <- as.character(orders)
orders

# obtain test MAE and RMSE using sliding window CV
# ------------------------------------------------
results_ar_sw <- sapply(orders, 
                   ts_cv, 
                   ts = total_firearm_deaths,
                   type = "linear",
                   xcovars = xcovariates,
                   include_mean = TRUE,
                   min_window = 6,
                   sliding_trainset = TRUE,
                   forward = 1,
                   verbose = FALSE)

plot_df_ar_sw <-  results_ar_sw |> 
                  t() |> 
                  as.data.frame() |> 
                  rownames_to_column("params") |>
                  mutate(mse = as.numeric(mse),
                         MAEs = as.numeric(MAEs))

myARMAplot(plot_df_ar_sw, "MSE by ARMA degree\n(sliding window CV)")


# obtain test MAE and RMSE using rolling window CV
# ------------------------------------------------
results_ar_fc <- sapply(orders, 
                   ts_cv, 
                   ts = total_firearm_deaths,
                   type = "linear",
                   xcovars = xcovariates,
                   include_mean = TRUE,
                   min_window = 6,
                   sliding_trainset = FALSE,
                   forward = 1,
                   verbose = FALSE)

plot_df_ar_fc <-  results_ar_fc |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column("params") |>
  mutate(mse = as.numeric(mse),
         MAEs = as.numeric(MAEs))

myARMAplot(plot_df_ar_fc, "MSE by ARMA degree\n(forward chaining CV)")


# chart of AR(0-4)
# chart of MA(0-4)
# chart of AR(0-4) with MA(0)
# chart of AR(0-4) with MA(1)
# chart of AR(0-4) with MA(2)
# chart of AR(0-4) with MA(3)
# chart of AR(0-4) with MA(4)




