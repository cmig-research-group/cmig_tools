createBasisSVD <- function(xvec, splineType = "ns", knots  = NULL, df = NULL, 
                              intercept  = TRUE, demean = FALSE, Boundary.knots = range(xvec), 
                              remove_col = NULL, dframe = TRUE, Xvars = NULL, 
							  method = NULL, xSubset = NULL, addConst = FALSE) 
{
  # Function to create basisFunction functions, given a vector of numeric values and
  # optional parameters
  # Inputs: 
  # xvec:           a vector containing the age/time or other values using 
  #                 which basisFunction functions will be created
  #
  # splineType:     can be one of the following:
  #                   * "ns" (default)
  #                   * "bs"
  #                   * "nsk"
  #
  # knots:          a vector of values specifying where knots should be placed;
  #                 if NULL, then knots are evenly placed 
  #
  # df:             instead of providing knots, df can be specified to create 
  #                 basisFunction functions (knots decided internally) 
  #
  # intercept:      logical, specifying if intercept (during creation of 
  #                 splines) should be included or not (default: TRUE) 
  #
  # demean:         logical, specifying if basisFunction functions should be column-wise 
  #                 mean centered (default: TRUE) 
  #
  # Boundary.knots: a vector of values specifying boundary knots (defaults to 
  #                 the range of the xvec) 
  #
  # remove_col:     specifies if one or more basisFunction functions should be dropped; 
  #                 can be one of the following:
  #                   * "first"  (first basisFunction function dropped)
  #                   * "middle" (middle basisFunction function dropped; default)
  #                   * "last"   (last basisFunction function dropped)
  #                   * "none"   (no basisFunction function dropped even if demean = TRUE)
  #                   * NULL     (middle basisFunction function dropped if demean = TRUE)
  #                   * one or more numbers 
  #
  # dframe:         logical, specifying if returned variable should be a data
  #                 frame (default: TRUE); if TRUE, then the column names are 
  #                 either "bf_x" or "bf_demean_x" (x is the basisFunction function 
  #                 number) and rows are named using the values in xvec (i.e., 
  #                 this particular output style assumes that xvec contains 
  #                 unique values)
  # Xvars:            vector or matrix of X variables that should be
  #                   (linearly) regressed out from the basis functions
  #
  # addConst:         logical; if true, a vector of ones is added as the 
  #                   first column of the output basisFunction
  #
  # method:           character; one of the following (see Notes):
  #                       * 'default'
  #                       * 'svd
  #
  # Output:
  # A data frame or a matrix containing the created basisFunction functions; if the
  # output is a data frame, then the column names are "bf_x" or "bf_demean_x"
  # (where x is the basisFunction function number) and row numbers are the values in
  # xvec (if xvec contains non-unique values such as full data, then specify 
  # dframe = FALSE)
  #
  # Requirements:
  # "splines" library for cases when splineType is "bs" or "ns" and "splines2" 
  # library when splineType is "nsk"
  #
  # Notes:
  # There are two types of methods supported:
  # default:          this returns the output as such from R
  #
  # svd:              first, the basis functions are created; then, if X
  #                   variables are specified, they are regressed out of the
  #                   basis functions. Next, the basis functions are
  #                   column-wise mean centered, followed by creating their
  #                   orthonormal basis (using singular value decomposition;
  #                   this might drop one or more basis functions); these are
  #                   then min-max scaled to have values between 0 and 1.
  
  # Get the required library

  if (splineType == "nsk")
  {
    require(splines2)
  } else
  {
    require(splines)
  }
  
  # Set default values for knots if not specified
  if (is.null(knots) & is.null(df))
  {
    knots = seq(from = min(xvec), to = max(xvec), length = 5)
    knots = knots[2:(length(knots) - 1)] # default will be 5 basisFunction functions with evenly spaced knots
  }
  
  # Create basisFunction functions
  if (splineType == "bs")
  {
    basisFunction <- bs(xvec, knots = knots, df = df, intercept = intercept, 
                degree = 3, Boundary.knots = Boundary.knots)
  } else
  {
    if (splineType == "ns")
    {
      basisFunction <- ns(xvec, knots = knots, df = df, intercept = intercept, 
                  Boundary.knots = Boundary.knots)
    } else
    {
      basisFunction <- nsk(xvec, knots = knots, df = df, intercept = intercept, 
                   Boundary.knots = Boundary.knots)
    }
  }
  
  # Save the number of columns for later use
  actCols <- dim(basisFunction)[2]
  
  # Demean the basisFunction functions (if required)
  if (demean)
  {
    # demean basisFunction functions: overwrite basisFunction
    basisFunction <- apply(basisFunction, 2, function(x) x - mean(x))
  }
  
  # Decide which column to drop (if required); defaults to middle; and drop
  if (demean & is.null(remove_col))
  {
    # Can result in fractions; for example median(1:6) = 3.5
    # remove_col <- round(median(1:dim(basis_demeaned)[2]))
    remove_col <- round(median(1:dim(basisFunction)[2]), digits = 0)
    basisFunction      <- basisFunction[, -remove_col]
  } else
  {
    # Handle case where remove_col is character
    if (is.character(remove_col))
    {
      remove_col <- tolower(remove_col)
      if (remove_col == "first")
      {
        remove_col <- 1
        basisFunction <- basisFunction[, -remove_col]
      } else
      {
        if (remove_col == "last")
        {
          remove_col <- ncol(basisFunction)
          basisFunction      <- basisFunction[, -remove_col]
        } else
        {
          if (remove_col == "middle")
          {
            remove_col <- round(median(1:dim(basisFunction)[2]), digits = 0)
            basisFunction      <- basisFunction[, -remove_col]
          } else
          {
            if (remove_col == "none")
            {
              basisFunction <- basisFunction # Do nothing
            } else
            {
              warning("Unknown value for remove_col provided; defaulting to middle")
              remove_col <- round(median(1:dim(basisFunction)[2]), digits = 0)
              basisFunction      <- basisFunction[, -remove_col]
            }
          }
        }
      }
    } else
    {
      if (as.logical(sum(remove_col >= 1 & remove_col <= ncol(basisFunction))))
      {
        basisFunction <- basisFunction[, -remove_col]
      } else
      {
        warning("Invalid/out of range value for remove_col provided; defaulting to middle")
        remove_col <- round(median(1:dim(basisFunction)[2]), digits = 0)
        basisFunction <- basisFunction[, -remove_col]
      }
    }
  }

  # Regress out Xvars if provided 
  if (!is.null(Xvars)) {
    lm_res <- lm(basisFunction ~ Xvars)
    basisFunction <- resid(lm_res)
  }

  # Perform orthonormalization if method is 'svd'
  if (method == 'svd') {
    # Demean the matrix (subtract the mean of each column)
	if (demean) {
      warning("Basis functions are already demeaned; orthogonalization will be performed on demeaned basis functions")
	  basis_demeaned <- basisFunction
    } else {
      basis_demeaned <- scale(basisFunction, scale = FALSE)
    }
    # Perform SVD
    svd_res <- svd(basis_demeaned)

    # Set a tolerance to determine when singular values are "zero" (i.e., dependent)
    tol <- 1e-10
    non_zero_singulars <- svd_res$d > tol  # Keep only non-zero (independent) components

    # Retain only independent columns 
    basis_orth <- svd_res$u[, non_zero_singulars] 

    # Rescale by maximum absolute value
    basis_orth_scaled <- apply(basis_orth, 2, function(x) x / max(abs(x)))

    basisFunction <- basis_orth_scaled
    # # Return orthogonalized and scaled basisFunction
    # return(list(basisFunction = basis_orth_scaled, rank = qr(basisFunction)$rank))

	# Do we need to add a constant?
    if (addConst) {
        basisFunction <-cbind(rep(1, length(xvec)), basisFunction)
    }
  }
  # Return
  if (!dframe)
  {
    return (basisFunction)
  } else
  {
    basisFunction <- as.data.frame(basisFunction, row.names = xvec)
    if (demean)
    {
      colnames(basisFunction) <- gsub('V', 'bf_demean_', colnames(basisFunction))
    } else if (method == 'svd') {
	  colnames(basisFunction) <- gsub('V', 'bf_svd_', colnames(basisFunction))
	} else {
	  colnames(basisFunction) <- gsub('V', 'bf_demean_', colnames(basisFunction))
	}
	
    return(basisFunction)
  }
}


extract_basis <- function(x, basisFunction) {
  basis_values = data.frame(matrix(ncol = dim(basisFunction)[2], nrow = 1))
  if (!is.na(x)){
    rownum = which.min(abs(as.numeric(rownames(basisFunction)) - x))
    basis_values = as.data.frame(basisFunction[rownum,])
  }
}


get_basis_values <- function(data, basisFunction, varname){
# This function takes, as inputs, a dataframe `data` and corresponding column
# name `varname` and the basisFunction function `basisFunction` returned by createBasis and
# returns the corresponding basisFunction function values for the values in
# data$varname; can be used to map basisFunction function values to a full dataset; note
# that this function does not perform any interpolation or extrapoliation
  require(data.table)
  
  # get basisFunction values
  basis_values <- lapply(data[,varname], extract_basis, basisFunction=basisFunction)
  
  # Combine basisFunction values as columns in dataframe
  basis_values_df <- rbindlist(lapply(basis_values, function(x) if(is.null(x[])) data.frame(matrix(ncol = dim(basisFunction)[2], nrow = 1)) else x), use.names=FALSE)
  
  return (basis_values_df)
  
}