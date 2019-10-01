# Author: Robert Hijmans
# January 2010
# License GPL3

varImportance <- function(mod, dat, vars, n=10) {
	RMSE <- matrix(nrow=n, ncol=length(vars))
	colnames(RMSE) <- vars
	for (i in 1:length(vars)) {
	rd <- dat
	v <- vars[i]
		for (j in 1:n) {
			rd[[v]] <- sample(rd[[v]])
			p <- predict(mod, rd)
			RMSE[j,i] <- predicts::RMSE(rd$SOC, p)
		}
	}
	return(RMSE)
}


partialResponse <- function(mod, dat, var, nsteps=25) {
	if (is.factor(data[[var]])) {
		steps <- levels(data[[var]]])
	} else {
		if (is.null(rng)) { 
			rng <- range(dat[[var]]) 
		}
		increment <- (rng[2] - rng[1])/(nsteps-2)
		steps <- seq(rng[1]-increment, rng[2]+increment, increment)
	}
	res <- rep(NA, length(steps))
	for (i in 1:length(steps)) {
		dat[[var]] <- steps[i]
		p <- predict(mod, dat)
		res[i] <- mean(p)
	}
	data.frame(var=steps, p=res)
}



