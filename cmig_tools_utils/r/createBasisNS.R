createBasisNS <- function(data, splinevar, knots = NULL, df = NULL, intercept = TRUE) {
  # args      <- commandArgs(trailingOnly = TRUE)
  # splinevar       <- args[1]
  # knots     <- args[2]
  # intercept <- as.logical(args[3])
  # outName   <- args[4]

  # Get the required library
  require(splines)
  require(pracma)

  # Create xvec
  xvec = linspace(min(data[,splinevar],na.rm = TRUE),max(data[,splinevar],na.rm = TRUE),101)

  # Create basis functions
  basis <- data.frame(ns(x=xvec,knots = knots, df=df,intercept = intercept), row.names = xvec)
  colnames(basis) = paste0('bf_',c(1:df)) 

  return (basis)
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
  basis_values_df <- rbindlist(lapply(basis_values, function(x) if(is.null(x)) data.frame(matrix(ncol = dim(basis)[2], nrow = 1)) else x), use.names=FALSE)

  return (basis_values_df)

}