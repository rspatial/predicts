# Author: Robert J. Hijmans
# Date :  December 2009
# Version 0.1
# Licence GPL v3



RMSE <- function(obs, prd, na.rm=FALSE) {
	sqrt(mean((obs - prd)^2, na.rm=na.rm))
}

RMSE_null <- function(obs, pred, na.rm=FALSE) {
	r <- RMSE(obs, pred, na.rm=na.rm)
	null <- RMSE(obs, mean(obs))
	(null - r) / null
}



.auctest <- function(p, a) {
	w <- stats::wilcox.test(p, a)
	pauc <- w$p.value
	auc <- as.vector(w$statistic) / (length(a) * length(p))
	cbind(auc, pauc)
}



setClass("paModelEvaluation",
	representation (
		presence = "vector",
		absence = "vector",
		confusion = "matrix",
		stats = "data.frame",
		tr_stats = "data.frame",
		thresholds = "data.frame"
	)
)



pa_evaluate <- function(p, a, tr) {
	p <- stats::na.omit(p)
	a <- stats::na.omit(a)
	np <- length(p)
	na <- length(a)
	if (na == 0 | np == 0) {
		stop("cannot evaluate a model without absence and presence data that are not NA")
	}

	if (missing(tr)) {
		if (length(p) > 1000) {
			tr <- as.vector(stats::quantile(p, 0:1000/1000))
		} else {
			tr <- p
		}
		if (length(a) > 1000) {
			tr <- c(tr, as.vector(stats::quantile(a, 0:1000/1000)))
		} else {
			tr <- c(tr, a)
		}
		tr <- sort(unique( round(tr, 8)))
		tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
	} else {
		tr <- sort(as.vector(tr))
	}
	
	N <- na + np

	xc <- methods::new("paModelEvaluation")
	xc@presence = p
	xc@absence = a
		
	R <- sum(rank(c(p, a))[1:np]) - (np*(np+1)/2)
	auc <- R / (as.numeric(na) * as.numeric(np))
	
	cr <- try( stats::cor.test(c(p,a), c(rep(1, length(p)), rep(0, length(a))) ), silent=TRUE )
	corc <- pcor <- NA
	if (class(cr) != "try-error") {
		corc <- cr$estimate
		pcor <- cr$p.value
	} 
	
	res <- matrix(ncol=4, nrow=length(tr))
	colnames(res) <- c("tp", "fp", "fn", "tn")
	for (i in 1:length(tr)) {
		res[i,1] <- length(p[p>=tr[i]])  # a  true positives
		res[i,2] <- length(a[a>=tr[i]])  # b  false positives
		res[i,3] <- length(p[p<tr[i]])    # c  false negatives
		res[i,4] <- length(a[a<tr[i]])    # d  true negatives
	}
	xc@confusion = res
	a = res[,1]
	b = res[,2]
	c = res[,3]
	d = res[,4]
# after Fielding and Bell	
	np <- as.integer(np)
	na <- as.integer(na)
	prevalence = (a[1] + c[1]) / N
 # overall diagnostic power
	ODP = (b[1] + d[1]) / N
	xc@stats <- data.frame(np, na, prevalence, auc, cor=corc, pcor, ODP)
	rownames(xc@stats) <- NULL
 # correct classification rate
	CCR = (a + d) / N
 # sensitivity, or true positive rate
	TPR = a / (a + c)
 # specificity, or true negative rate
	TNR = d / (b + d)
 # False positive rate
	FPR = b / (b + d)
 # False negative rate
	FNR = c/(a + c)
	PPP = a/(a + b)
	NPP = d/(c + d)
 # misclassification rate
	MCR = (b + c)/N
 # odds ratio
	OR = (a*d)/(c*b)

	prA = (a+d)/N
	prY = (a+b)/N * (a+c)/N
	prN = (c+d)/N * (b+d)/N
	prE = prY + prN
	kappa = (prA - prE) / (1-prE)
	xc@tr_stats <- data.frame(treshold=tr, kappa, CCR, TPR, TNR, FPR, FNR, PPP, NPP, MCR, OR)



	max_kappa <- tr[which.max(kappa)]
	# maximum sum of the sensitivity (true positive rate) and specificity (true negative rate)
	max_spec_sens <- tr[which.max(TPR + TNR)]
	# no omission
	no_omission <- tr[max(which(res[, "fn"] == 0))]
	# Suggestions by Diego Nieto-Lugilde		
	# equal prevalence
	equal_prevalence = tr[which.min(abs(tr - prevalence))] 
	# equal sensitivity and specificity
	equal_sens_spec <- tr[which.min(abs(TPR - TNR))]
		
	xc@thresholds <- data.frame(max_kappa, max_spec_sens, no_omission,  equal_prevalence, equal_sens_spec)

	return(xc)
}




setMethod ("show" , "paModelEvaluation", 
	function(object) {
		cat("@stats\n")
		print(round(object@stats, 3))
		cat("\n")
		cat("@thresholds\n")
		print(round(object@thresholds, 3))
		cat("\n")
		cat("@tr_stats\n")
		x <- rbind(head(object@tr_stats, 4), tail(object@tr_stats, 3))
		x <- round(x, 2)
		x[4, ] <- "..."
		print(x)		
	}
)	



setMethod("plot", signature(x="paModelEvaluation"), 
	function(x, y="ROC", col="red", ...) {
		if (y == "boxplot") {
			prediction <- c(x@presence, x@absence)
			group <- c(rep("presence", length(x@presence)), rep("absence", length(x@absence)) )
			boxplot(prediction~group, data=data.frame(prediction, group), ...)
		} else if (y == "ROC") {
			txt = paste("AUC=", round(x@stats$auc,3))
			plot(x@tr_stats$FPR, x@tr_stats$TPR, xlim=c(0,1), ylim=c(0,1), xlab="False postive rate", ylab="True positive rate", col=col, main=txt, ...)
			lines(x@tr_stats$FPR, x@tr_stats$TPR, col=col)
			lines(rbind(c(0,0), c(1,1)), lwd=2, col="grey")
		} else if (y == "density") {
			pr <- density(x@presence)
			ab <- density(x@absence, bw=pr$bw )
			yl = c(min(ab$y, pr$y), max(ab$y, pr$y))
			xl = c(min(x@tr_stats$treshold), max(x@tr_stats$treshold))
			plot(ab, main="", ylab=paste("Density. Bandwidth=",round(pr$bw,5),paste=""), xlab="predicted value", xlim=xl, ylim=yl, lwd=2, lty=2, col="blue", ...)
			lines(pr, col="red", lwd=2)
		} else if (y %in% colnames(x@tr_stats)) {
			dat <- x@tr_stats[, y]
			#mx = x@tr_stats$treshold[which.max(dat)]
			#mn = x@tr_stats$treshold[which.min(dat)]
			plot(x@tr_stats$treshold, dat, xlab="threshold", ylab=y, ...)
			lines(x@tr_stats$treshold, dat, col=col)
		} else {
			stop("invalid option")
		}
	}	
)
