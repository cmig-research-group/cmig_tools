args       <- commandArgs(trailingOnly = TRUE)
valvec     <- args[1]
knots      <- args[2]
splineType <- args[3]
intercept  <- as.logical(args[4])
outName    <- args[5]
callPath   <- args[6]

# Note that this script is designed to be called via MATLAB; 
# regression is not performed and all modifications of the basis functions are 
# done in MATLAB

# Read data
valvec <- as.matrix(read.table(file = valvec, header = F))

# Read knots
knots <- as.matrix(read.table(file = knots, header = F))

# Call
source(file.path(callPath, "createBasis.R"))
basis <- createBasis(valvec, knots    = knots, splineType = splineType, 
                             method   = "raw", intercept  = intercept, 
                             doSubset = FALSE, doInterp   = FALSE)
# Write output
write.table(basis, file = outName, sep = "\t", row.names = F, col.names = F, na = "NaN")