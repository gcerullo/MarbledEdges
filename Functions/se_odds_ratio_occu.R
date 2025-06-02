

#' Compute Standard Error of odds ratio
#'
#' Calculates the difference in predicted probabilities (occupancy or detection)
#' between two observations in a new dataset, along with the associated standard error,
#' from a fitted `unmarked` occupancy model.
#'
#' @param newdat A data frame containing the covariate values for prediction. Must contain the two rows we want to contrast.
#' @param fit A fitted model object of class `unmarkedFitOccu` (from the `unmarked` package).
#' @param type A character string indicating which component to analyze: `"state"` for occupancy or `"det"` for detection. Defaults to `"state"`.
#' @param order A numeric vector of length 2 specifying the row indices in `newdat` to compare. 
#' The function returns the predicted probability of the index listed in the first element in `order` 
#' minus the second. Defaults to `c(2, 1)`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`estim`}{The estimated difference in predicted probabilities.}
#'   \item{`se`}{The standard error of the estimated difference.}
#' }
#'
#' @details
#' The function constructs the appropriate model matrix for either the occupancy or detection model,
#' computes the predicted probabilities using the inverse-logit transformation, and then calculates the difference
#' between the two specified observations. The standard error is computed using the delta method via numerical
#' differentiation (`numDeriv::grad`).
#'
#' @importFrom numDeriv grad
#' @importFrom stats coef vcov model.matrix as.formula
#'
#' @examples
#' \dontrun{
#' library(unmarked)
#' data(frogs)
#' umf <- unmarkedFrameOccu(y = frogs[, 1:3], siteCovs = frogs[, 4:6])
#' fit <- occu(~ detect1 + detect2 ~ site1 + site2, data = umf)
#' newdat <- frogs[1:2, 4:6]
#' se_diff_occ(newdat, fit, type = "state", order = c(2, 1))
#' }
#'
#' @export
#' 
se_odds_ratio_occ <- function(newdat, fit, type = "state", order = c(2,1)){
  
  # model matrix for detection
  f_p <- deparse(fit@formula[[2]])
  if(length(f_p) > 1){
    f_p <- paste(f_p, collapse = "")
  }
  X_p <- model.matrix(
    as.formula(f_p), 
    newdat
  )
  
  # model matrix for occupancy
  f_psi <- deparse(fit@formula[[3]])
  if(length(f_psi) > 1){
    f_psi <- paste(f_psi, collapse = "")
  }
  X_psi <- model.matrix(
    as.formula(paste("~", f_psi)),
    newdat
  )
  
  # extract parameter estimates
  estims <- coef(fit)
  
  # construct difference in probability function
  or_occu <- function(x, X, order = order){
    preds <- as.double(X %*% x)
    #return(preds[order[1]] - preds[order[2]])
    # Calculate odds ratio
    or <- exp(preds[order[1]] - preds[order[2]])
    return(or)
  }
  
  # now get the estimates and model matrix depending on type
  if(type == "state"){
    indexes <- 1:ncol(X_psi)
    pred <- or_occu(estims[indexes], X_psi, order)
    gradient <- numDeriv::grad(or_occu, estims[indexes], X = X_psi, order = order)
    V <- vcov(fit)[indexes, indexes]
  } else if(type == "det"){
    indexes <- (ncol(X_psi) + 1):(ncol(X_psi) + ncol(X_p))
    pred <- or_occu(estims[indexes], X_p, order)
    gradient <- numDeriv::grad(or_occu, estims[indexes], X = X_p, order = order)
    V <- vcov(fit)[indexes, indexes]
  } else {
    stop("Unrecognized type argument. Expecting `state` or `det`.")
  }
  
  # define list of returned values
  return(
    list(
      estim = pred,
      se = sqrt(gradient %*% V %*% gradient) |>
        as.double()
    )
  )
  
  
}
  
  
