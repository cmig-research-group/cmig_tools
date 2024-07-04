createBasis <- function(xvec, splineType = "ns", knots  = NULL, df = NULL, 
                              intercept  = TRUE, demean = TRUE, Boundary.knots = range(xvec), 
                              remove_col = NULL, dframe = TRUE)
{
  # Function to create basis functions, given a vector of numeric values and
  # optional parameters
  # Inputs: 
  # xvec:           a vector containing the age/time or other values using 
  #                 which basis functions will be created
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
  #                 basis functions (knots decided internally) 
  #
  # intercept:      logical, specifying if intercept (during creation of 
  #                 splines) should be included or not (default: TRUE) 
  #
  # demean:         logical, specifying if basis functions should be column-wise 
  #                 mean centered (default: TRUE) 
  #
  # Boundary.knots: a vector of values specifying boundary knots (defaults to 
  #                 the range of the xvec) 
  #
  # remove_col:     specifies if one or more basis functions should be dropped; 
  #                 can be one of the following:
  #                   * "first"  (first basis function dropped)
  #                   * "middle" (middle basis function dropped; default)
  #                   * "last"   (last basis function dropped)
  #                   * "none"   (no basis function dropped even if demean = TRUE)
  #                   * NULL     (middle basis function dropped if demean = TRUE)
  #                   * one or more numbers 
  #
  # dframe:         logical, specifying if returned variable should be a data
  #                 frame (default: TRUE); if TRUE, then the column names are 
  #                 either "bf_x" or "bf_demean_x" (x is the basis function 
  #                 number) and rows are named using the values in xvec (i.e., 
  #                 this particular output style assumes that xvec contains 
  #                 unique values)
  #
  # Output:
  # A data frame or a matrix containing the created basis functions; if the
  # output is a data frame, then the column names are "bf_x" or "bf_demean_x"
  # (where x is the basis function number) and row numbers are the values in
  # xvec (if xvec contains non-unique values such as full data, then specify 
  # dframe = FALSE)
  #
  # Requirements:
  # "splines" library for cases when splineType is "bs" or "ns" and "splines2" 
  # library when splineType is "nsk"
  
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
    knots = knots[2:(length(knots) - 1)] # default will be 5 basis functions with evenly spaced knots
  }
  
  # Create basis functions
  if (splineType == "bs")
  {
    basis <- bs(xvec, knots = knots, df = df, intercept = intercept, 
                degree = 3, Boundary.knots = Boundary.knots)
  } else
  {
    if (splineType == "ns")
    {
      basis <- ns(xvec, knots = knots, df = df, intercept = intercept, 
                  Boundary.knots = Boundary.knots)
    } else
    {
      basis <- nsk(xvec, knots = knots, df = df, intercept = intercept, 
                   Boundary.knots = Boundary.knots)
    }
  }
  
  # Save the number of columns for later use
  actCols <- dim(basis)[2]
  
  # Demean the basis functions (if required)
  if (demean)
  {
    # demean basis functions: overwrite basis
    basis <- apply(basis, 2, function(x) x - mean(x))
  }
  
  # Decide which column to drop (if required); defaults to middle; and drop
  if (demean & is.null(remove_col))
  {
    # Can result in fractions; for example median(1:6) = 3.5
    # remove_col <- round(median(1:dim(basis_demeaned)[2]))
    remove_col <- round(median(1:dim(basis)[2]), digits = 0)
    basis      <- basis[, -remove_col]
  } else
  {
    # Handle case where remove_col is character
    if (is.character(remove_col))
    {
      remove_col <- tolower(remove_col)
      if (remove_col == "first")
      {
        remove_col <- 1
        basis <- basis[, -remove_col]
      } else
      {
        if (remove_col == "last")
        {
          remove_col <- ncol(basis)
          basis      <- basis[, -remove_col]
        } else
        {
          if (remove_col == "middle")
          {
            remove_col <- round(median(1:dim(basis)[2]), digits = 0)
            basis      <- basis[, -remove_col]
          } else
          {
            if (remove_col == "none")
            {
              basis <- basis # Do nothing
            } else
            {
              warning("Unknown value for remove_col provided; defaulting to middle")
              remove_col <- round(median(1:dim(basis)[2]), digits = 0)
              basis      <- basis[, -remove_col]
            }
          }
        }
      }
    } else
    {
      if (as.logical(sum(remove_col >= 1 & remove_col <= ncol(basis))))
      {
        basis <- basis[, -remove_col]
      } else
      {
        warning("Invalid/out of range value for remove_col provided; defaulting to middle")
        remove_col <- round(median(1:dim(basis)[2]), digits = 0)
        basis      <- basis[, -remove_col]
      }
    }
  }
  
  # Return
  if (!dframe)
  {
    return (basis)
  } else
  {
    basis <- as.data.frame(basis, row.names = xvec)
    if (demean)
    {
      colnames(basis) <- paste0("bf_demean_", setdiff(1:actCols, remove_col))
    } else
    {
      colnames(basis) <- paste0("bf_", setdiff(1:actCols, remove_col))
    }
    return(basis)
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
# This function takes, as inputs, a dataframe `data` and corresponding column
# name `varname` and the basis function `basis` returned by createBasis and
# returns the corresponding basis function values for the values in
# data$varname; can be used to map basis function values to a full dataset; note
# that this function does not perform any interpolation or extrapoliation
  require(data.table)
  
  # get basis values
  basis_values <- lapply(data[,varname], extract_basis, basis=basis)
  
  # Combine basis values as columns in dataframe
  basis_values_df <- rbindlist(lapply(basis_values, function(x) if(is.null(x[])) data.frame(matrix(ncol = dim(basis)[2], nrow = 1)) else x), use.names=FALSE)
  
  return (basis_values_df)
  
}