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


varImportance <- function(model, y, x, n=10, stat, ...) {

#	vars <- vars[vars %in% colnames(x)]
#	if (length(vars) < 1) {
#		stop("no valid names in vars")
#	}
	vars <- colnames(x)
	eva <- matrix(nrow=n, ncol=length(vars))
	colnames(eva) <- vars

	if (missing(x)) {
		x <- .get_model_data(model)
		if (is.null(x)) {
			stop("data argument cannot be missing when using this model type")
		}
	}

	P <- predict(model, x, ...)
	if (is.factor(P)) {
		if (missing(stat)) stat = "overall"
		stopifnot(stat %in% c("overall", "kappa"))
		efun <- function(y, x) {
			tab <- table(x, y)
			1 - cm_evaluate(tab, stat)
		}
	} else {
		if (missing(stat)) stat = "RMSE"
		stopifnot(stat %in% c("RMSE", "AUC", "cor"))
		if (stat == "AUC") {
			efun <- function(y, x) {
				i <- y == 1
				1 - pa_evaluate(x[i], x[!i])@stats$auc
			}
		} else if (stat == "cor"){
			efun <- function(y, x) {
				i <- y == 1
				1 - pa_evaluate(x[i], x[!i])@stats$cor
			}
		} else {
			efun <- predicts::RMSE
		}
	}

	base <- efun(y, P)

	for (i in 1:length(vars)) {
		rd <- x
		v <- vars[i]
		for (j in 1:n) {
			rd[[v]] <- sample(rd[[v]])
			p <- predict(model, rd, ...)
			eva[j,i] <- efun(y, p)
		}
	}
	colMeans(eva) - base 
}


.pnrnc <- function(nr, nc, nl) {
	if (missing(nc)) {
		nc <- ceiling(sqrt(nl))
	} else {
		nc <- max(1, min(nl, round(nc)))
	}
	if (missing(nr)) {
		nr <- ceiling(nl / nc)
	} else {
		nr <- max(1, min(nl, round(nr)))
		nc <- ceiling(nl / nr)
	}
	c(nr, nc)
}


partialResponse <- function(model, data, var=NULL, rng=NULL, nsteps=25, plot=TRUE, nr, nc, ...) {

	if (missing(data)) {
		data <- .get_model_data(model)
		if (is.null(data)) {
			stop("data argument cannot be missing when using this model type")
		}
	}
	
	if (is.null(var)) {
		var <- names(data)
	} else if (is.numeric(var)) {
		var <- round(var)
		stopifnot(all(var > 0 & var <= ncol(data)))
		var <- names(data)[var]
	} else {
		stopifnot(all(var %in% names(data)))
	}
	
	out <- lapply(var, function(v) {
		if (is.factor(data[[v]])) {
			steps <- levels(data[[v]])
		} else {
			if (is.null(rng)) { 
				rng <- range(data[[v]]) 
			}
			increment <- (rng[2] - rng[1])/(nsteps-2)
			steps <- seq(rng[1]-increment, rng[2]+increment, increment)
		}
		res <- rep(NA, length(steps))
		for (i in 1:length(steps)) {
			d <- data
			# to handle factors (#16)
			d[d[[v]] != steps[i], v] <- steps[i]
			p <- predict(model, d, ...)
			res[i] <- mean(p)
		}
		x <- data.frame(steps, res)
		names(x) <- c(v, "p")
		x
	})
	names(out) <- var
	if (plot) {
		nrnc <- .pnrnc(nr, nc, length(var))
		old.par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old.par))
		graphics::par(mfrow=nrnc)
		for (i in 1:length(out)) {
			plot(out[[i]], type="l", las=1)
		}
		invisible(out)
	} else {
		out
	}
}


partialResponse2 <- function(model, data, var1, var2, var2levels, rng=NULL, nsteps=25, ...) {
	if (is.factor(data[[var1]])) {
		steps <- levels(data[[var1]])
	} else {
		if (is.null(rng)) { 
			rng <- range(data[[var1]]) 
		}
		increment <- (rng[2] - rng[1])/(nsteps-2)
		steps <- seq(rng[1]-increment, rng[2]+increment, increment)
	}
	res <- rep(NA, length(steps))
	out <- data.frame(var1=steps)	
	for (v in var2levels) {
		data[[var2]] <- v
		for (i in 1:length(steps)) {
			# to handle factors (#16)
			data[data[[var1]] != steps[i], var1] <- steps[i]
			##data[[var1]] <- steps[i]
			p <- stats::predict(model, data, ...)
			res[i] <- mean(p)
		}
		out[[paste(var2, v, sep="_")]] <- res
	}
	out
}


