
setClass("SDM",
	contains = "VIRTUAL",
	representation (
		presence = "data.frame",
		absence = "data.frame",
		hasabsence = "logical"
	),	
	prototype (	
		hasabsence = FALSE
	),
	validity = function(object)	{
		(ncol(object@presence) == 2) && (ncol(object@absence) == 2)
	}
)	



setMethod ("show" , "SDM", 
	function(object) {
		cat("class    :" , class(object), "\n")
		cat("variables:", colnames(object@presence), "\n")
		pp <- nrow(object@presence)
		cat("presence points:", pp, "\n")
		if (pp < 10) {
			print(object@presence)
		} else {
			print(object@presence[1:6,])
			cat("      ... \n")
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

