x <- c(1,2,3,4,5)
y <- c(9,1,7,2,3)
h <- 3

###  local_polyB <- function(formula, x, y, kernel = "epinechnikov", h, x0) { 
###    
###    # extract elements from the formula
###    RHS <- attr(terms(formula), which = "term.labels") |> stringr::str_remove_all(" ")
###    deg <- stringr::str_extract(RHS, "\\,[0-9]+\\)") |> stringr::str_extract("[0-9]+") |> as.numeric()
###    
###    # Kernel smoother
###    if(kernel == "epinechnikov") { 
###        
###        yhat <- vector("numeric", length = length(x0))
###        for(i in 1:length(x0)) {
###          t <- abs(x-xi)/h
###          w <- ifelse(abs(t)<=1, 3/4*(1-t^2), 0)
###          new_x <- x[w!=0]
###          new_y <- y[w!=0]
###          
###          # not with lm
###          Z <- model.matrix(~., poly(new_x, degree = deg, raw = TRUE))
###          if(length(w[w!=0])<2) { 
###            W <- matrix(w[w != 0], ncol = 1, nrow = 1)
###          } else{ 
###            W <- matrix(diag(w), ncol = length(w)) 
###          }
###          bhat <- solve(t(Z)%*%W%*%Z) %*% t(Z)%*%W%*%new_y
###          z0 <- model.matrix(~., poly(xi, degree = deg, raw = TRUE)) 
###          yhat[i] <- z0%*%bhat
###          
###        }
###    }else {
###      
###      yhat <- vector("numeric", length = length(x0))
###      for(i in 1:length(x0)) {
###        xi <- x0[i]
###        w <- iflese((xi-x)<=h, 1, 0)
###        new_x <- x[w!=0]
###        new_y <- y[w!=0]
###        
###        # in the case when h only allows in one observation, pad with 2 zeros
###        if(length(w[w!=0])<2) { 
###          w <- c(0, w[w!=0], 0) 
###          new_x <- c(0, new_x, 0)
###          new_y <- c(0, new_y, 0)
###        }      
###        
###        # with no lm
###        Z <- model.matrix(~., poly(new_x, degree = deg, raw = TRUE))
###        W <- matrix(diag(w[w != 0]), ncol = length(w[w!=0]))
###        bhat <- solve(t(Z)%*%W%*%Z) %*% t(Z)%*%W%*%new_y
###        z0 <- model.matrix(~., poly(xi, degree = deg, raw = TRUE)) 
###        yhat[i] <- z0%*%bhat     
###        
###        # with lm
###        lm <- lm(new_y ~ poly(new_x), degree= deg, raw = TRUE)
###        yhat[i] <- predict(lm, new_data = xi)
###      }
###    }
###    return(yhat)
###  }

local_polyB <- function(formula, x, y, kernel = "epinechnikov", h, x0) {
  
  # Kernel smoother
  if(kernel == "epinechnikov") { 
    t <- abs(x-x0)/h
    w <- ifelse(abs(t)<=1, 3/4*(1-t^2), 0)
  }else {
    w <- ifelse((x0-x)<=h, 1, 0)
  }
  
  # extract elements from the formula
  RHS <- attr(terms(formula), which = "term.labels") |> stringr::str_remove_all(" ")
  deg <- stringr::str_extract(RHS, "\\,[0-9]+\\)") |> stringr::str_extract("[0-9]+") |> as.numeric()
  
  # remove values outside support
  new_x <- x[w!=0]
  new_y <- y[w!=0]
  if(length(new_x) < 2) { print(paste0("too few observations inside bandwidth = ", h, ". Try a larger h."))}
  stopifnot(length(new_x)>1) 
  
  # build Z, W, and bhat
  Z <- model.matrix(~., poly(new_x, degree = deg, raw = TRUE))
  W <- matrix(diag(w[w!=0]), ncol = length(w[w!=0]))
  bhat <- solve(t(Z)%*%W%*%Z) %*% t(Z)%*%W%*%new_y
  
  # prediction
  z0 <- model.matrix(~., poly(x0, degree = deg, raw = TRUE)) 
  yhat <- z0%*%bhat
  return(yhat[,1])
}



epinechnikov <- function(x0, x, h) {
  t <- abs(x-x0)/h
  d <- ifelse(abs(t)<=1, 3/4*(1-t^2), 0)
  #d[d != 0] #<- so d is the same length as 'ys_of_xs_within_h'
}



epinechnikov <- function(x0, x, h) {
  
  d <- vector("numeric", length = length(x))
  for(i in 1:length(x)) {
    
    t <- abs(x[i]-x0)/h
    if(abs(t)<=1) {
      d[i] <- 3/4*(1-t^2) 
    } else d[i] <- 0
  }
  d[d != 0]
}




nad_wat <- function(x, y, h, kernel = "epinechnikov") {
  
  yhat <- vector("numeric", length = length(x))
  for(i in 1:length(x)){
    k <- epinechnikov(x[i], x, h)
    ys_of_xs_within_h <- ifelse(x %in% x[(x<x[i]+h & x[i]-h<x)], y, 0)
    stopifnot(length(y) == length(k))
    yhat[i] <- sum(k*ys_of_xs_within_h)/sum(k)
  }
  o <- data.frame("yhat" = yhat, "x" = x)
  o
}


local_poly <- function(formula, X, y, kernel = "gauss", h, x0) {
  
  # Kernel Smoother
  if(kernel == "gauss"){
    w <- exp(-((x0-x)^2)/(2*h^2) )
  }else if(kernel == "epinechnikov") { 
    t <- abs(x-x0)/h
    w <- ifelse(abs(t)<=1, 3/4*(1-t^2), 0)
  }else {
    w <- iflese((x0-x)<=h, 1, 0)
  }
  
  # extract elements from the formula
  RHS <- attr(terms(formula), which = "term.labels") |> stringr::str_remove_all(" ")
  deg <- stringr::str_extract(RHS, "\\,[0-9]+\\)") |> stringr::str_extract("[0-9]+") |> as.numeric()

  # build Z, W, and bhat
  Z <- model.matrix(~., poly(x, degree = deg))
  W <- matrix(diag(w), ncol = length(w))
  bhat <- solve(t(Z)%*%W%*%Z) %*% t(Z)%*%W%*%y

  # prediction
 ### if(length(x0)<=deg) {
 ###   z0 <- matrix(NA, nrow = length(x0), ncol = deg + 1)
 ###   for(i in 0:deg) {
 ###     z0[, i+1] <- x0^i
 ###   }
 ### } else {
 ###   z0 <- model.matrix(~., poly(x0, degree = deg)) 
 ### }
  lenx0 <- length(x0)
  if(lenx0 <= deg) { x0 <- c(x0, x0+.01, x0-.01)}
  z0 <- model.matrix(~., poly(x0, degree = deg)) 
  yhat <- z0%*%bhat
  print(yhat)
  return(yhat[1:lenx0,1])

}


