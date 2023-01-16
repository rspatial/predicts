
setClass("SDM",
	contains = "VIRTUAL",
	representation (
		presence = "matrix",
		absence = "matrix",
		hasabsence = "logical"
	),	
	prototype (	
		presence = matrix(nrow=0, ncol=2),
		absence = matrix(nrow=0, ncol=2),
		hasabsence = FALSE
	),
	validity = function(object)	{
		(ncol(object@presence) == 2) && (ncol(object@absence) == 2)
	}
)	



setMethod ("show" , "SDM", 
	function(object) {
		cat("class    :" , class(object), "\n\n")
		cat("variables:", colnames(object@presence), "\n\n")
		pp <- nrow(object@presence)
		cat("\npresence points:", pp, "\n")
		if (pp < 10) {
			print(object@presence)
		} else {
			print(object@presence[1:10,])
			cat("  (... ...  ...)\n\n")
		}
		if (object@hasabsence) {
			pp <- nrow(object@absence)
			cat("\nabsence points:", pp, "\n")
			if (pp < 10) {
				print(object@absence)
			} else {
				print(object@absence[1:10,])
				cat("  (... ...  ...)\n\n")
			}
		}
	}
)	

