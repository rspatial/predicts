
if (!isGeneric("hullModel")) {setGeneric("hullModel", function(p, ...)	standardGeneric("hullModel"))}
if (!isGeneric("geometry")) { setGeneric("geometry", function(x, ...) standardGeneric("geometry"))}

if (!isGeneric("MaxEnt")) { setGeneric("MaxEnt", function(x, p, ...) standardGeneric("MaxEnt"))}	
if (!isGeneric("mess")) { setGeneric("mess", function(x, ...) standardGeneric("mess")) }	

if (!isGeneric("plot")) { setGeneric("plot", function(x, y,...) standardGeneric("plot"))}
if (!isGeneric("predict")) { setGeneric("predict", function(object, ...) standardGeneric("predict"))}	

