createBasis <- function(valvec, knots = NULL, splineType = "nsk", Xpowers = NULL,
                                method = "svd", minMax = NULL, intercept = TRUE, 
                                addConst = FALSE, doSubset = TRUE, doInterp = TRUE) 
{
  # Function to create basis functions, given a set of values, knot values 
  # and other parameters - supports creation of natural cubic splines (ns), 
  # B-splines (bs), and natural cubic splines with unit heights at knots
  # (nsk): ns and bs functions are created using the splines package in R 
  # while nsk is created using the splines2 package in R
  #
  # The approach for creating basis functions:
  # 1) Calculate min and max for the values in valvec
  # 2) Create a linearly spaced vector of 100 values between the min and max
  # 3) Create basis functions
  # 4) Linearly regress out the powers of X variables (if required)
  # 5) Perform SVD on the demeaned residuals / basis functions (if required)
  # 6) Perform interpolation and scale data back to the values in valvec
  #
  # The powers in step 4 refers to the powers of the vector created in step 2
  #
  # Inputs: 
  # valvec:         vector of values for which basis functions need to be
  #                 created (for example, age / time)
  #
  # knots:          vector of values which serve as knots
  #
  # splineType:     can be one of the following:
  #                   * "nsk"
  #                   * "ns"
  #                   * "bs"
  #
  # Xpowers:        vector indicating which powers of the linearly spaced
  #                 variables should be regressed from the basis functions;
  #                 for example, 0:1 would mean regressing out the
  #                 intercept and the linear effect of the variable; 0:2
  #                 would mean regressing out the intercept, the linear
  #                 effect, and the quadratic effect of the variable
  #
  # method:         can be one of the following:
  #                   * "svd"
  #                   * "raw"
  #
  # minMax:         numeric vector: if a single value is provided, this is
  #                 the number of linearly spaced values that will be
  #                 created; if two values are provided, these are the
  #                 minimum and the maximum values which are used for
  #                 creating a vector linearly spaced values; if three
  #                 values are provided, the third value is used as the
  #                 number of linearly spaced values; if left empty, the
  #                 min and max of valvec is used
  #
  # intercept:      logical, specifying if intercept (during creation of 
  #                 splines) should be included or not (default: TRUE) 
  #
  # addConst:       logical; if true, a vector of ones is added as the 
  #                 first column of the output basisFunction
  #
  # doSubset:       logical; if true, a subset vector of linearly spaced values
  #                 is created, which is then used for creating the basis 
  #                 function; set to FALSE if the input valvec is already the 
  #                 subset of values
  #
  # doInterp:       logical; if true, the basis function subset is interpolated
  #                 so that the basis function values are corresponding to the 
  #                 values in valvec; set to FALSE to get the basisSubset as the 
  #                 returned value
  #
  # Output:
  # A  matrix containing the created basis functions; if doInterp = FALSE, then
  # the output is the basis subset for the linearly spaced values
  #
  # Requirements:
  # "splines" library for cases when splineType is "bs" or "ns" and "splines2" 
  # library when splineType is "nsk"
  #
  # Notes:
  # There are two types of methods supported:
  # svd:              first, the basis functions are created; then, the
  #                   powers of the uniformly spaced vector are regressed
  #                   out; next, the residualized basis functions are
  #                   column-wise mean centered, followed by creating their
  #                   orthonormal basis (using singular value decomposition;
  #                   this might drop one or more basis functions); these are
  #                   then scaled to have values between 0 and 1; finally,
  #                   the values are interpolated to generate the basis
  #                   function for the full data
  # 
  # raw:              this returns the output as such; the basis functions
  #                   are computed for the linearly spaced vector, powers are
  #                   regressed out (if the user specifies so), and then the
  #                   values are interpolated to generate the basis functions
  #                   for the full data
  
  # Get the required library
  if (splineType == "nsk")
  {
    require(splines2)
  } else
  {
    require(splines)
  }
  # require(pracma) # pracma is required for interp1
  
  # Make sure valvec is a vector
  if (is.data.frame(valvec))
  {
    valvec <- as.matrix(valvec)
  }
  
  # Genetate a vector of linearly spaced values
  # If called via MATLAB, doSubset is FALSE; the valvec provided is then
  # used as such to generate basis functions; interpolation should be off in that case
  if (doSubset)
  {
    if (is.null(minMax))
    {
      minX    <- min(valvec)
      maxX    <- max(valvec)
      howMany <- 101
    } else
    {
      if (length(minMax) == 1)
      {
        minX    <- min(valvec)
        maxX    <- max(valvec)
        howMany <- minMax
      } else
      {
        if (length(minMax) == 2)
        {
          minX    <- min(minMax)
          maxX    <- max(minMax)
          howMany <- 101
        } else
        {
          if (length(minMax) == 3)
          {
            howMany <- minMax[3]
            minMax  <- minMax[1:2]
            minX    <- min(minMax)
            maxX    <- max(minMax)
          }
        }
      }
    }
    Xvars <- seq(from = minX, to = maxX, length.out = howMany)
  } else
  {
    Xvars <- valvec
  }
  
  # Determine knots
  if (is.null(knots))
  {
    knots <- quantile(valvec, c(0.25, 0.50, 0.75))
  }
  
  # Create basis functions
  if (splineType == "bs")
  {
    basisSubset <- bs(Xvars, knots = knots, intercept = intercept)
  } else
  {
    if (splineType == "ns")
    {
      basisSubset <- ns(Xvars, knots = knots, intercept = intercept)
    } else
    {
      if (splineType == "nsk")
      {
        basisSubset <- nsk(Xvars, knots = knots, intercept = intercept)
      } else
      {
        stop("Invalid value of splineType provided; should be ns, bs, or nsk")
      }
    }
  }
  
  # Regress out Xpowers, if provided 
  if (!is.null(Xpowers))
  {
    # Expand Xvars by raising it to Xpowers
    Xexpanded <- outer(Xvars, Xpowers, FUN = "^")
    
    # Do not explicitly add an intercept as it could be part of Xexpanded
    # if one of the powers is zero
    mdl <- lm(basisSubset ~ -1 + Xexpanded)
    
    # Extract residuals
    basisSubset <- residuals(mdl)
  }
  
  # Perform orthon-ormalization if method is 'svd'
  if (method == 'svd')
  {
    # Demean the matrix (subtract the mean of each column)
    basisSubset <- scale(basisSubset, scale = FALSE)
    
    # Perform SVD
    svd_res <- svd(basisSubset)
    
    # Drop columns
    # tol <- 1e-10
    # non_zero_singulars <- svd_res$d > tol  # Keep only non-zero (independent) components
    if (is.null(Xpowers))
    {
      # The number of columns in Xvars is zero; can subtract one; otherwise need to transpose
      colsToRetain <- 1:(ncol(basisSubset) - 1)
    } else
    {
      colsToRetain <- 1:(ncol(basisSubset) - ncol(Xexpanded))
    }
    basisSubset <- svd_res$u[, colsToRetain]
    
    # Rescale by maximum absolute value
    basisSubset <- apply(basisSubset, 2, function(x) x / max(abs(x)))
  }
  
  # Interpolate to full dataset
  if (doInterp)
  {
    basisFunction <- matrix(data = NA, nrow = nrow(valvec), ncol = ncol(basisSubset))
    for (ii in 1:ncol(basisSubset))
    {
      vals <- approx(x = Xvars, y = basisSubset[,ii], xout = valvec, method = "linear")
      basisFunction[,ii] <- vals$y
    }
    # Solution using pracma
    # for (ii in 1:ncol(basisSubset))
    # {
    #   basisFunction[,ii] <- interp1(Xvars, basisSubset[,ii], valvec, method = "linear")
    # }
  }
  
  # Do we need to add a constant?
  if (addConst)
  {
    basisFunction <-cbind(rep(1, nrow(basisFunction)), basisFunction)
  }
  
  # Return results
  if (doInterp)
  {
    return (basisFunction)  
  } else
  {
    return(basisSubset)
  }
}

#   # Return
#   if (!dframe)
#   {
#     return (basisFunction)
#   } else
#   {
#     basisFunction <- as.data.frame(basisFunction, row.names = xvec)
#     if (demean)
#     {
#       colnames(basisFunction) <- gsub('V', 'bf_demean_', colnames(basisFunction))
#     } else if (method == 'svd') {
#       colnames(basisFunction) <- gsub('V', 'bf_svd_', colnames(basisFunction))
#     } else {
#       colnames(basisFunction) <- gsub('V', 'bf_demean_', colnames(basisFunction))
#     }
#     
#     return(basisFunction)
#   }
# }

# # DEPRECATED
# extract_basis <- function(x, basisFunction) {
#   basis_values = data.frame(matrix(ncol = dim(basisFunction)[2], nrow = 1))
#   if (!is.na(x)){
#     rownum = which.min(abs(as.numeric(rownames(basisFunction)) - x))
#     basis_values = as.data.frame(basisFunction[rownum,])
#   }
# }
# 
# 
# get_basis_values <- function(data, basisFunction, varname){
#   # This function takes, as inputs, a dataframe `data` and corresponding column
#   # name `varname` and the basisFunction function `basisFunction` returned by createBasis and
#   # returns the corresponding basisFunction function values for the values in
#   # data$varname; can be used to map basisFunction function values to a full dataset; note
#   # that this function does not perform any interpolation or extrapoliation
#   require(data.table)
#   
#   # get basisFunction values
#   basis_values <- lapply(data[,varname], extract_basis, basisFunction=basisFunction)
#   
#   # Combine basisFunction values as columns in dataframe
#   basis_values_df <- rbindlist(lapply(basis_values, function(x) if(is.null(x[])) data.frame(matrix(ncol = dim(basisFunction)[2], nrow = 1)) else x), use.names=FALSE)
#   
#   return (basis_values_df)
#   
# }