args       <- commandArgs(trailingOnly = TRUE)
splineType <- args[1]
xvec       <- args[2]
knots      <- args[3]
intercept  <- as.logical(args[4])
df         <- as.logical(args[5])
outName    <- args[6]
callPath   <- args[7]

# Read data
xvec <- as.matrix(read.table(file = xvec, header = F))

# Read knots
knots <- as.matrix(read.table(file = knots, header = F))

# Decide if knots or df
if (df)
{
  df    <- knots
  knots <- NULL
} else
{
  df <- NULL
}

# Call
source(file.path(callPath, "createBasis.R"))
basis <- createBasis(xvec, splineType = splineType, knots = knots, df = df, 
                     intercept  = intercept,  demean = FALSE, Boundary.knots = range(xvec), 
                     remove_col = NULL, dframe = FALSE)
  
# Write output
write.table(basis, file = outName, sep = "\t", row.names = F, col.names = F, na = "NaN")