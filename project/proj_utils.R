ts_cv <- function(
    ts, # ts object
    type = "linear", # type of prediction function c("linear, "boost_trees)
    params = NA, # params of ARMA model when running linear reg or set of parameters when running boosted trees
    xcovars, # matrix of covariates, use cbind to construct. eg, cbind(df$x1, df$x2, ..., df$xk)
    include_mean, # include intercept in ARIMA model TRUE/FALSE
    min_window = 5, # minimum length of training set K
    forward = 1, # number of forward periods when forecasting
    sliding_trainset = TRUE, # TRUE = sliding window; FALSE = expanding window scheme
    seasonal = NA,
    verbose = FALSE
    ) {
    stopifnot(class(ts) == "ts")
    ts_start <- tsp(ts)[1] # the starting point of ts t
    ts_end <- tsp(ts)[2] # the end point of ts
    ts_freq <- tsp(ts)[3] # the frequency of ts
    ts_length <- length(ts) # The total length of time series T

    # CV 
    n_iters <- ts_length - min_window
    mae <- matrix(NA, nrow = n_iters, ncol = forward) # initialize a matrix to store performance metrics
    mse <- matrix(NA, nrow = n_iters, ncol = forward) 
    for (i in 1:n_iters) {
      # compute the starting and end points of training and test sets
      if (sliding_trainset == TRUE) {
        train_start <- ts_start + (i - 1) / ts_freq
      } else {
        train_start <- ts_start
      }
      train_end <- ts_start + (i - 1 + min_window - 1) / ts_freq
      test_start <- ts_start + (min_window + i - 1) / ts_freq
      test_end <- min(ts_start + (min_window + i - 1 + forward - 1) / ts_freq, ts_end)
      
      # construct training and test sets
      y_train <- window(ts, start = train_start, end = train_end)
      y_test <- window(ts, start = test_start, end = test_end)
      
      # also truncate covariates time series to the same length
      x_train <- window(xcovars, start = train_start, end = train_end)
      x_test <- window(xcovars, start = test_start, end = test_end)
      
      # run model on training set
      if(type == "linear") {
          train_fit <- forecast::Arima(y_train, 
                                       order=params, 
                                       xreg=x_train, 
                                       include.mean=include_mean, 
                                       seasonal=seasonal,
                                       method = "ML")
     
          pred <- forecast::forecast(train_fit, h = forward, xreg = x_test)
          pred <- pred$mean
      }
      
      if(type == "rf") {
          ntrees <- params[1]
          mtry <- params[2]
          #maxdepth <- params[3]
          
          x_train_df  <- as.data.frame(x_train) |> `colnames<-`(letters[1:ncol(x_train)])
          x_test_df   <- as.data.frame(x_test) |> `colnames<-`(letters[1:ncol(x_test)])
          y_train_df <- y_train |> as.vector()
          
          fit <- randomForest(x = x_train_df, y = y_train_df, 
                              ntrees = ntrees, 
                              mtry = mtry)
          pred <- predict(fit, newdata = x_test_df)
      }
      
      # compute performance metrics on test set
      mae[i , 1:length(y_test)] <-   abs(pred - y_test) #test set for i-fold is of size=1
      mse[i, 1:length(y_test)] <-   (pred - y_test)^2 #sqrt(pred - y_test)^2 
    }

    # average performance metrics across iterations
    if(verbose == TRUE) {
      print(mae)
      print(mse)
    }
    mae <- colMeans(mae, na.rm = TRUE)
    mse <- colMeans(mse, na.rm = TRUE)
    results <- list("MAEs" = mae, "mse" = mse)
    return(results)
}


make_table <- function(model, names, model_name) {
  mse   <- sqrt(model$residuals^2) |> mean()
  aic    <- model$aic
  coeffs <- model$coef |> unname() 
  SEs    <- diag(model$var.coef) |> unname() |> sqrt() 
  N      <- length(model$residuals)
  names  <- names
  
  mytable <- data.frame(Names = names,
                        Coefficients = round(coeffs, 1) |> as.character(), 
                        `Std.Errors` = paste0("(", as.character(round(SEs, 2)), ")")
  ) |> tidyr::pivot_longer(c("Coefficients","Std.Errors"))
  
  row1 <- data.frame(Names = "AIC" , name = "", value = as.character(round(aic, 1)))
  row2 <- data.frame(Names = "mse", name = "", value = as.character(round(mse, 1)))
  row3 <- data.frame(Names = "N"   , name = "", value = as.character(N))
  
  mytable <- mytable |> dplyr::bind_rows(row1, row2, row3)
  names(mytable) <- c("Parameters", "Type", model_name)
  mytable
}

myARMAplot <- function(data, mytitle = "") {
  
  # reshape data
  pd <- data[1:5, ] |> 
    rbind(data[c(1,6,11,16,21), ]) |>
    rbind(data[6:25, ]) |>
    `row.names<-`(NULL) |>
    mutate(i = rep(0:4,6),
           model = c(rep("AR(i)",5), 
                     rep("MA(i)",5), 
                     rep("ARMA(1,i)",5),
                     rep("ARMA(2,i)",5),
                     rep("ARMA(3,i)",5),
                     rep("ARMA(4,i)",5))
    )
  
  # build plot
  ggplot(pd, aes(x=i, y=mse, color=model)) + 
    geom_point() + 
    geom_line() + 
    ggtitle(mytitle) +
    theme(#panel.grid = element_blank(),
          #panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    geom_text(aes(y = ifelse(i==4, mse, NA), 
                  x = 4, 
                  label = model, 
                  vjust = -0.5), size = 2.5)
  
}


myRFplot <- function(data, mytitle = "") {
  
  # reshape data
  pd <- data |>
        mutate(p = rep(1:4, 3),
               model = c(rep("5 trees",4), 
                         rep("25 trees",4), 
                         rep("100 trees",4))
  )
  
  # build plot
  ggplot(pd, aes(x=p, y=mse, color=model)) + 
    geom_point() + 
    geom_line() + 
    ggtitle(mytitle) +
    theme(#panel.grid = element_blank(),
          #panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    geom_text(aes(y = ifelse(p==4, mse, NA), 
                  x = 4, 
                  label = model, 
                  vjust = -0.5), size = 2.5)
  
}









