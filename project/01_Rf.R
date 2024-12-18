library(tidyverse)
library(randomForest)

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
wd <- wd |> filter(state == "Texas") |> filter(year >= 1999 & year < 2015)
names(wd)[names(wd) == "total_hunting_licenses_tags_permits_and_stamps"] <- "hunt_licenses"


# declare outcome variable
# ------------------------
total_firearm_deaths <- ts(1e4*wd$total_firearm_deaths/wd$population)

# declare covariates
# ------------------
wd$pct_risk_age <- wd$pct_45_to_54_years + 
                   wd$pct_55_to_59_years + 
                   wd$pct_60_to_64_years

xcovariates <- cbind(wd$permit, 
                     wd$poverty_rate, 
                     wd$hunt_licenses, 
                     wd$unemployment_rate,
                     wd$cpi,
                     wd$pct_risk_age
                    )


# CV each Boosted Trees model
# ---------------------------
n <- 0
params_set <- vector(mode = "list") #, length = 32)
ntrees <- c(5, 25, 100) 
mtry <- c(1, 2, 3, 4)
for(j in ntrees) {
  for(i in mtry) {
      n <- n + 1
      params_set[[n]] <- c(i,j)
  }
}
names(params_set) <- as.character(params_set)
params_set


# obtain test MAE and RMSE with sliding window CV
# -----------------------------------------------
results_rf_sw <- sapply(params_set, 
                   ts_cv, 
                   ts = total_firearm_deaths,
                   type = "rf",
                   xcovars = xcovariates,
                   include_mean = TRUE,
                   min_window = 6,
                   sliding_trainset = TRUE,
                   forward = 1,
                   verbose = FALSE)

plot_df_rf_sw <- results_rf_sw |>
                 t() |>
                 as.data.frame() |> 
                 rownames_to_column("params") |>
                 mutate(mse = as.numeric(mse),
                        MAEs = as.numeric(MAEs))

myRFplot(plot_df_rf_sw, 
         "MSE by Random Forest parameter set\n(sliding window CV)")


# obtain test MAE and RMSE with rolling window CV
# -----------------------------------------------
results_rf_fc <- sapply(params_set, 
                   ts_cv, 
                   ts = total_firearm_deaths,
                   type = "rf",
                   xcovars = xcovariates,
                   include_mean = TRUE,
                   min_window = 6,
                   sliding_trainset = FALSE,
                   forward = 1,
                   verbose = FALSE)

plot_df_rf_fc <- results_rf_fc |>
  t() |>
  as.data.frame() |> 
  rownames_to_column("params") |>
  mutate(mse = as.numeric(mse),
         MAEs = as.numeric(MAEs))

myRFplot(plot_df_rf_fc, 
         "MSE by Random Forest parameter set\n(forward chaining CV)")


# Testing RF function
# -------------------

from <- 1
to <- 6
test <- from + to

x_train <- xcovariates_boost[from:to, ] |> as.data.frame()
y_train <- total_firearm_deaths[from:to]
x_test <- xcovariates_boost[test, ] |> t() |> as.data.frame()

fit_tree <- rpart::rpart(y_train~., 
                         data = cbind(y_train, x_train), 
                         control = list(maxdepth=3))
yhat_tree <- predict(fit_tree, newdata = x_test)

fit_rf <- randomForest(x = x_train, y = y_train, ntrees = 50, mtry = 1)
yhat_rf <- predict(fit_rf, newdata = x_test)

yhat_tree
yhat_rf[[1]]
total_firearm_deaths[test]


