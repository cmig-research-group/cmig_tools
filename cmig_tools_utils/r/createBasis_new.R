createBasis_new <- function(xvec, knots = NULL, splineType = 'ns', Xvars = NULL, 
                                   dfFlag = FALSE, intercept = TRUE, method = 'svd',
                                   xSubset = NULL, addConst = FALSE, outDir = getwd(),
                                   cleanUp = TRUE, instance = 1, dframe = TRUE) 
  {
  # This function creates basis functions given age/time values and knots.
  # It supports natural cubic splines or B-splines using ns or bs functions from 
  # the splines package, or creating natural cubic splines with unit heights at knots 
  # using the nsk function from the splines2 package
  
  # Inputs:
  # xvec:           a vector containing the age/time or other values using 
  #                 which basis functions will be created
  #
  # knots:          vector of values that serve as knots
  #
  # splineType:     character; one of the following: 
  #                     * 'ns' (default)
  #                     * 'bs'  (B-splines)
  #                     * 'nsk' (natural cubic splines with unit heights at knots)
  #
  # Xvars:          vector or matrix of variables to be regressed out from the basis functions (optional)
  # 
  # dfFlag:         logical; specifies whether degrees of freedom should be used instead of knots; 
  #                 enter a single number corresponding to the degrees of freedom for "knots" and 
  #                 specify "dfFlag" as true
  #
  # intercept:      logical; if TRUE, an intercept will be used during spline creation (default TRUE)
  # 
  # method:         character; one of the following: 
  #                     * 'default' (return basis functions) 
  #                     * 'svd' (perform orthonormalization)
  #
  # xSubset:        a subset of x values to be used for creating basis functions
  #
  # addConst:       logical; if TRUE, a vector of ones is added as the first column of the output basis function
  #
  # outDir:         full path to where the output file should be 
  #                 temporarily saved; if empty, pwd is used
  #
  # cleanUp:        logical; if true, deletes all temporary files that were
  #                 created along the way
  #
  # instance:       numeric; useful if function is being called in parallel
  #                 to ensure independence of each call (and that files are
  #                 not deleted, etc.)
  #
  # Outputs:
  # basisFunction:  a matrix containing the spline basis functions
  #
  # bfRank:         rank of the basis functions (useful as a sanity check)
  #
  # basisSubset:    if xSubset was specified, contains the basis functions for the subset prior to interpolation
  #
  # settings:       a list of settings used for creating the basis functions
  #
  # timing:         timing information for various steps in the function
  #
  # Notes:
  # There are two types of methods supported:
  # default:        this returns the output as such from R
  #
  # svd:            first, the basis functions are created; then, if X
  #                 variables are specified, they are regressed out of the
  #                 basis functions. Next, the basis functions are
  #                 column-wise mean centered, followed by creating their
  #                 orthonormal basis (using singular value decomposition;
  #                 this might drop one or more basis functions); these are
  #                 then min-max scaled to have values between 0 and 1.

  # Requirements:
  # "splines" library for cases when splineType is "bs" or "ns" and "splines2" 
  # library when splineType is "nsk"
  # "data.table" library 
  
  # Get the required library
  if (splineType == "nsk") {
    if (!require(splines2)) {
        install.packages(splines2)  # Install splines2 if not already installed
        library(splines2)             # Load splines2
    }
  } else {
    if (!require(splines)) {
     install.packages(splines)   # Install splines if not already installed
     library(splines)              # Load splines
   }
  }

  if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
  }

  # Initialize timing
  timing <- list()
  tInit <- Sys.time()

  # Input validation
  if (is.null(xvec) || length(xvec) == 0) stop('Please provide a vector of age')
  
  xvec <- as.vector(xvec)

  # Set default values for knots if not specified and dfFlag is FALSE
  if (is.null(knots) && !dfFlag) {
    # Generate evenly spaced knots between the min and max of xvec
    knots <- seq(from = min(xvec), to = max(xvec), length.out = 5)
    # Remove the first and last points to place the knots within the range
    knots <- knots[2:(length(knots) - 1)]  # This creates 3 interior knots
  }

  if (!splineType %in% c('ns', 'bs', 'nsk')) stop("splineType should be 'ns', 'bs', or 'nsk'")

  if (!is.null(Xvars) && nrow(Xvars) != length(xvec)) stop('Mismatch between number of entries in xvec and rows in Xvars')

  if (dfFlag && length(knots) > 1) stop('When dfFlag is TRUE, only one value for knots should be provided.')
  
  # Handle optional xSubset for interpolation
  if (!is.null(xSubset)) {
    xSubset <- na.omit(xSubset)  # Remove NA values from xSubset
    useSubset <- TRUE
  } else {
    useSubset <- FALSE
  }

  # Create basis functions based on spline type using xvec as input
  if (splineType == 'ns') {
    if (dfFlag) {
      basisFunction <- ns(xvec, df = knots, intercept = intercept)
    } else {
      basisFunction <- ns(xvec, knots = knots, intercept = intercept)
    }
  } else if (splineType == 'bs') {
    if (dfFlag) {
      basisFunction <- bs(xvec, df = knots, intercept = intercept)
    } else {
      basisFunction <- bs(xvec, knots = knots, intercept = intercept)
    }
  } else if (splineType == 'nsk') {
    basisFunction <- nsk(xvec, knots = knots, intercept = intercept)
  }

  # Add constant (vector of ones) if addConst is TRUE
  if (addConst) {
    basisFunction <- cbind(1, basisFunction)
  }

  # Interpolation if xSubset is used
  if (useSubset && length(xSubset) < length(xvec)) {
    basisFunction_full <- matrix(NA, nrow = length(xvec), ncol = ncol(basisFunction))
    for (i in 1:ncol(basisFunction)) {
      basisFunction_full[, i] <- approx(xSubset, basisFunction[, i], xout = xvec)$y
    }
    basisFunction <- basisFunction_full
  }

  # Regress out Xvars if provided
  if (!is.null(Xvars)) {
    lm_res <- lm(basisFunction ~ Xvars)
    basisFunction <- resid(lm_res)
  }

  # Perform orthonormalization if method is 'svd'
  #browser()
  if (method == 'svd') {
    # Demean the matrix (subtract the mean of each column)
    basis_demeaned <- scale(basisFunction, scale = FALSE)

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
    # # Return orthogonalized and scaled basis
    # return(list(basisFunction = basis_orth_scaled, rank = qr(basisFunction)$rank))
  }

  # Return basis function as a formatted data frame
  if (!dframe) {
    return(basisFunction)
  } else {
    basisFunction <- as.data.frame(basisFunction, row.names = xvec)
    if (method == 'svd') {
      colnames(basisFunction) <- paste0("bf_svd_", 1:ncol(basisFunction))
    } else {
      colnames(basisFunction) <- paste0("bf_", 1:ncol(basisFunction))
    }
  }

  # # Calculate rank of the basis functions
  # bfRank <- qr(basisFunction)$rank

  # # Check for rank deficiency and issue a warning if necessary
  # if (bfRank < ncol(basisFunction)) {
  #   warning("Rank deficient basis function: bfRank = ", bfRank, 
  #         " is less than the number of basis functions = ", ncol(basisFunction))
  # }   

  # # Record timing for overall execution
  # timing$tOverall <- Sys.time() - tInit

#   # Output settings
#   settings <- list(splineType = splineType, knots = knots, dfFlag = dfFlag, 
#                    intercept = intercept, method = method, addConst = addConst, instance = instance)

#   # Return results
#   list(basisFunction = basisFunction, bfRank = bfRank, settings = settings, timing = timing)
# }

# Return only the basis function (formatted as a data frame)
  return(basisFunction)
}

# Helper function to extract basis function values for a single data point
extract_basis <- function(x, basis) {
  basis_values <- data.frame(matrix(ncol = dim(basis)[2], nrow = 1))
  if (!is.na(x)) {
    rownum <- which.min(abs(as.numeric(rownames(basis)) - x))
    basis_values <- as.data.frame(basis[rownum, ])
  }
  return(basis_values)
}

# Function to map basis function values to a dataset
get_basis_values <- function(data, basis, varname) {
  # This function maps basis function values to a dataset based on the variable in 'varname'
  
  # get basis values for each row in the data
  basis_values <- lapply(data[, varname], extract_basis, basis = basis)
  
  # Combine basis values as columns in dataframe
  basis_values_df <- rbindlist(lapply(basis_values, function(x) 
    if (is.null(x[])) data.frame(matrix(ncol = dim(basis)[2], nrow = 1)) else x), use.names = FALSE)
  
  return(basis_values_df)
}