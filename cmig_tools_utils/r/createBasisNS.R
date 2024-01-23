createBasisNS <- function(xvec, knots = NULL, df = NULL, intercept = TRUE, demean = TRUE, Boundary.knots = range(xvec), remove_col = NULL) {
  # args      <- commandArgs(trailingOnly = TRUE)
  # splinevar       <- args[1]
  # knots     <- args[2]
  # intercept <- as.logical(args[3])
  # outName   <- args[4]

  # Get the required library
  require(splines)

  # set default values for knots if not specified
  if (is.null(knots) & is.null(df)){
    knots = linspace(min(xvec), max(xvec), 6) 
    knots = knots[2:(length(knots)-1)] # default will be 5 basis functions with evenly spaced knots
  }

  # Create basis functions
  basis <- data.frame(ns(xvec, df, knots, intercept, Boundary.knots), row.names = xvec)
  colnames(basis) = paste0('bf_',c(1:dim(basis)[2]))

  if (demean) {
    # demean basis functions
    basis_demeaned = sapply(basis,function(x)x-mean(x))
    colnames(basis_demeaned) <- paste0("bf_demean_",c(1:dim(basis)[2]))
    row.names(basis_demeaned) = row.names(basis)

    # resulting matrix is rank deficient; remove middle by default
    if (is.null(remove_col)){
      remove_col = round(median(1:dim(basis_demeaned)[2]))
    }
    
    return(basis_demeaned[,-remove_col])
  } else {
    return (basis) 
  }
}

extract_basis <- function(x, basis) {
  basis_values = data.frame(matrix(ncol = dim(basis)[2], nrow = 1))
  if (!is.na(x)){
    rownum = which.min(abs(as.numeric(rownames(basis)) - x))
    basis_values = as.data.frame(basis[rownum,])
  }
}

get_basis_values <- function(data, basis, varname){

  require(data.table)

  # get basis values
  basis_values <- lapply(data[,varname], extract_basis, basis=basis)

  # Combine basis values as columns in dataframe
  basis_values_df <- rbindlist(lapply(basis_values, function(x) if(is.null(x[])) data.frame(matrix(ncol = dim(basis)[2], nrow = 1)) else transpose(x)), use.names=FALSE)

  return (basis_values_df)

}