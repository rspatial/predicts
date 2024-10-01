# Author: Robert Hijmans
# January 2010
# License GPL3

.get_model_data <- function(m) {
	if (inherits(m, "lm") || inherits(m, "glm")) {
		m$model
	} else if (inherits(m, "SDM")) {
		rbind(m@presence, m@absence)
	} else {
		NULL
	}
}


varImportance <- function(model, data, vars=colnames(data), n=10, ...) {
	RMSE <- matrix(nrow=n, ncol=length(vars))
	colnames(RMSE) <- vars

	if (missing(data)) {
		data <- .get_model_data(model)
		if (is.null(data)) {
			stop("data argument cannot be missing when using this model type")
		}
	}

	P <- predict(model, data, ...)
	for (i in 1:length(vars)) {
		rd <- data
		v <- vars[i]
		for (j in 1:n) {
			rd[[v]] <- sample(rd[[v]])
			p <- predict(model, rd)
			RMSE[j,i] <- predicts::RMSE(P, p)
		}
	}
	colMeans(RMSE) 
}


partialResponse <- function(model, data, var=1, rng=NULL, nsteps=25, ...) {

	if (missing(data)) {
		data <- .get_model_data(model)
		if (is.null(data)) {
			stop("data argument cannot be missing when using this model type")
		}
	}
	
	if (is.numeric(var)) {
		stopifnot(var > 0 & var <= ncol(data))
		var <- names(data)[var]
	} else {
		stopifnot(all(var %in% names(data)))
	}
	
	if (is.factor(data[[var]])) {
		steps <- levels(data[[var]])
	} else {
		if (is.null(rng)) { 
			rng <- range(data[[var]]) 
		}
		increment <- (rng[2] - rng[1])/(nsteps-2)
		steps <- seq(rng[1]-increment, rng[2]+increment, increment)
	}
	res <- rep(NA, length(steps))
	for (i in 1:length(steps)) {
		# to handle factors (#16)
		data[data[[var]] != steps[i], var] <- steps[i]

		p <- predict(model, data, ...)
		res[i] <- mean(p)
	}
	x <- data.frame(steps, res)
	names(x) <- c(var, "p")
	x
}


partialResponse2 <- function(model, data, var, var2, var2levels, rng=NULL, nsteps=25, ...) {
	if (is.factor(data[[var]])) {
		steps <- levels(data[[var]])
	} else {
		if (is.null(rng)) { 
			rng <- range(data[[var]]) 
		}
		increment <- (rng[2] - rng[1])/(nsteps-2)
		steps <- seq(rng[1]-increment, rng[2]+increment, increment)
	}
	res <- rep(NA, length(steps))
	out <- data.frame(var=steps)	
	for (v in var2levels) {
		data[[var2]] <- v
		for (i in 1:length(steps)) {
			# to handle factors (#16)
			data[data[[var]] != steps[i], var] <- steps[i]
			##data[[var]] <- steps[i]
			p <- stats::predict(model, data, ...)
			res[i] <- mean(p)
		}
		out[[paste(var2, v, sep="_")]] <- res
	}
	out
}


