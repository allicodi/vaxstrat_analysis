
SL.earth.mod <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
          ncross = 1, minspan = 0, endspan = 0, ...) 
{
  #.SL.require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                              ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if (family$family == "binomial") {
    
    # Check if all(Y) == 1 or == 0
    if(all(Y == 1)){
      fit.earth <- list(msg = "All Y == 1") 
      class(fit.earth) <- c("SL.earth.mod")
      pred <- predict(fit.earth, newdata = newX, type = "response")
    } else if (all(Y == 0)){
      fit.earth <- list(msg = "All Y == 0")   
      class(fit.earth) <- c("SL.earth.mod")
      pred <- predict(fit.earth, newdata = newX, type = "response")
    } else{
      fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
                                nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                ncross = ncross, minspan = minspan, endspan = endspan, 
                                glm = list(family = binomial))
      pred <- predict(fit.earth, newdata = newX, type = "response")
    }
  }
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth.mod")
  return(out)
}
predict.SL.earth.mod <- function (object, newdata, ...) 
{
  #.SL.require("earth")
  if(is.null(object$msg)){
    pred <- predict(object$object, newdata = newdata, type = "response")
  } else if(object$msg == "All Y == 1"){
    pred <- rep(1, nrow(newdata))
  } else if(object$msg == "All Y == 0"){
    pred <- rep(0, nrow(newdata))
  } else{
    pred <- NA
  }
  
  return(pred)
}