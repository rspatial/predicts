

RMSE <- function(obs, prd, na.rm=FALSE) {
	sqrt(mean((obs - prd)^2, na.rm=na.rm))
}

RMSE_null <- function(obs, pred, na.rm=FALSE) {
	r <- RMSE(obs, pred, na.rm=na.rm)
	null <- RMSE(obs, mean(obs))
	(null - r) / null
}

