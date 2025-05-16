#model calculator 

model <- readRDS("Models/final_model_5thMay2025.rds")
source("scripts/02_OrganiseMurreletData.R")

#read in means and SDs of covariates from original model to enable coorrect scalings
meansAndSds



predict_fun <- function(h100, h2000, e100, e2000, PC1_t1, dist) {
  # Scale the input parameters using the means and standard deviations
  scaleHabAmount100 <- (h100 - meansAndSds$meanHabAmountDich100) / meansAndSds$sdHabAmountDich100
  scaleHabAmount2000 <- (h2000 - meansAndSds$meanHabAmountDich2000) / meansAndSds$sdHabAmountDich2000
  scaleEdgeDens100 <- (e100 - meansAndSds$meanEdgeDens100) / meansAndSds$sdEdgeRook100
  scaleEdgeDens2000 <- (e2000 - meansAndSds$meanEdgeDens2000) / meansAndSds$sdEdgeRook2000
  PC1_t1 <- PC1_t1
  scaleCoastDist  <- (dist - meansAndSds$meanCoastDist ) / meansAndSds$sdCoastDist 
  
  
  # Create a data frame for predictions
  predict_df <- data.frame(
    scaleHabAmount100 = scaleHabAmount100,
    scaleHabAmount2000 = scaleHabAmount2000,
    scaleEdgeDens100 = scaleEdgeDens100,
    scaleEdgeDens2000 = scaleEdgeDens2000,
    PC1_t1 = PC1_t1,
    scaleCoastDist = scaleCoastDist
  )
  
  # Generate predictions using the model
  predictions <- predict(
    model,
    newdata = predict_df,
    type = "state",  # Replace "state" with the correct type for your model if needed
    se.fit = TRUE    # Obtain standard errors for predictions
  )
  
  # Add predictions and confidence intervals to the data frame
  predict_df$Occupancy <- predictions$Predicted
  # predict_df$SE <- predictions$se.fit
  # predict_df$LowerCI <- predict_df$Occupancy - 1.96 * predict_df$SE
  # predict_df$UpperCI <- predict_df$Occupancy + 1.96 * predict_df$SE
  
  return(predict_df$Occupancy)
}

# Test occupancy 
define_params <- list(
  h100 = 1,
  h2000 = 0.73135,
  e100 = 0.3323,
  e2000 = 0.2513,
  PC1_t1 = 1.04,
  dist =31700
)

PC1_quantiles
define_params <- list(
  h100 = 1,
  h2000 = 0.91743,
  e100 = 0.1680,
  e2000 = 0.09609,
  PC1_t1 = 0.973  ,
  dist =28693
)

results <- do.call(predict_fun, define_params)
results
