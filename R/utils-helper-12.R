


#' An internal function to construct b-spline basis matrix
#' 
#' @param x A numeric vector
#' @param degree An integer 
#' @param knots A vector
#' @param bknots A vector
#' @param calcderiv A real number
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_bsp <- function(x, 
                   knots, 
                   bknots, 
                   degree, 
                   calcderiv) {
  Nobs       <- length(x)
  fullknots  <- c(rep(bknots[1], degree+1), knots, rep(bknots[2], degree+1))
  Nintervals <- length(fullknots) - 1
  
  M1         <- matrix(0, Nobs, Nintervals)
  for (i in 1:Nintervals) {
    M1[, i] <- as.numeric(fullknots[i] <= x  & x < fullknots[i + 1])
  }
  lastknot_index <- which(x == bknots[2])
  if (length(lastknot_index) > 0) {
    M1[lastknot_index, (Nintervals - degree):Nintervals] <- 1
  }
  
  for (p in 1:degree) {
    M2 <- matrix(1, Nobs, Nintervals - p)
    for (i in 1:(Nintervals - p)) {
      if (fullknots[i + p] == fullknots[i]) {
        C1 <- 0
      } else {
        C1 <- (x - fullknots[i]) / (fullknots[i + p] - fullknots[i])
      }
      
      if (fullknots[i + p + 1] == fullknots[i + 1]) {
        C2 <- 0
      } else {
        C2 <- (fullknots[i + p + 1] - x) / (fullknots[i + p + 1] - 
                                              fullknots[i + 1])
      }
      M2[, i] <- C1 * M1[, i] + C2 * M1[, i + 1]
    }
    if(p != degree) M1 <- M2
  }
  
  splinevars <- M2 
  if(!calcderiv) {
    out <- splinevars
  }
 
  if(calcderiv) {
    deriv <- matrix(NA, Nobs, Nintervals - degree)
    for (i in 1:(Nintervals - degree)) {
      if (fullknots[i + degree] == fullknots[i]) {
        C1 <- 0
      } else {
        C1 <- degree / (fullknots[i + degree] - fullknots[i])
      }
      if (fullknots[i + degree + 1] == fullknots[i + 1]) {
        C2 <- 0
      } else {
        C2 <- degree / (fullknots[i + degree + 1] - fullknots[i + 1])
      }
      deriv[, i] <- C1 * M1[, i] - C2 * M1[, i + 1]
    }
    
    dsplinevars <- deriv # [, 2:ncol(deriv)]
    out <- dsplinevars
  }
  return(out)
}



#' An internal function to construct natural cubic spline basis matrix
#' 
#' @param x A numeric vector
#' @param knots A vector
#' @param bknots A vector
#' @param intercept An integer
#' @param calcderiv A real number
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_nsp_nsk <- function(x, 
                       knots, 
                       bknots, 
                       degree, 
                       intercept, 
                       calcderiv, 
                       normalize, 
                       preH) {
  
  Nintk     <- length(knots) 
  Nk        <- Nintk + 2
  allknots  <- c(rep(bknots[1], 4), knots,rep(bknots[2], 4))
  Nintervals <- length(allknots) - 1
  bs        <- matrix(0, length(x), Nk)
  bsderiv   <- matrix(0, length(x), Nk)
  # Compute B-spline values and derivatives
  if (calcderiv) {
    bs      <- GS_bsp(x = x, 
                      knots = knots, 
                      bknots = bknots, 
                      degree = degree, 
                      calcderiv = 0)
    bsderiv <- GS_bsp(x = x, 
                      knots = knots, 
                      bknots = bknots, 
                      degree = degree, 
                      calcderiv = 1)
  } else {
    bs      <- GS_bsp(x = x, 
                      knots = knots, 
                      bknots = bknots, 
                      degree = degree, 
                      calcderiv = 0)
  }
  
  # Handle boundary cases (for values below and above the boundary knots)
  x_below_boundary <- sum(x < bknots[1])
  x_above_boundary <- sum(x > bknots[2])
  if (x_below_boundary || x_above_boundary) {
    bs_bknots      <- matrix(0, 2, Nk)
    bsderiv_bknots <- matrix(0, 2, Nk)
    bs_bknots      <- GS_bsp(x = bknots, 
                             knots = knots, 
                             bknots = bknots, 
                             degree = degree, 
                             calcderiv = 0)
    bsderiv_bknots <- GS_bsp(x = bknots, 
                             knots = knots, 
                             bknots = bknots, 
                             degree = degree, 
                             calcderiv = 1)
   
    if (x_below_boundary) {
      xselect    <- which(x < bknots[1])
      Nxselect   <- length(xselect)
      below_azx1 <- matrix(rep(bs_bknots[1,], Nxselect), 
                           byrow = T, nrow = Nxselect) 
      below_azx2 <- matrix(rep(bsderiv_bknots[1,], Nxselect), 
                           byrow = T, nrow = Nxselect) 
      below_azx3 <- (x[xselect] - bknots[1]) 
      bs[xselect, ] <- below_azx1 + below_azx2 * below_azx3
      if (calcderiv) {
        for(i in xselect) bsderiv[i, ] <- bsderiv_bknots[1, ]
      }
    }
    if (x_above_boundary) {
      xselect     <- which(x >= bknots[2])
      Nxselect    <- length(xselect)
      above_azx1  <- matrix(rep(bs_bknots[2,], Nxselect), 
                           byrow = T, nrow = Nxselect)
      above_azx2  <- matrix(rep(bsderiv_bknots[2,], Nxselect), 
                           byrow = T, nrow = Nxselect) 
      above_azx3  <- (x[xselect] - bknots[2])  
      assignabove <- above_azx1 + above_azx2 * above_azx3
      bs[xselect, ] <- above_azx1 + above_azx2 * above_azx3
       if (calcderiv) {
        for(i in xselect) bsderiv[i, ] <- bsderiv_bknots[2, ]
       }
    } # if (x_above_boundary) {
  } # if (x_below_boundary || x_above_boundary) {
  
  # H <- GS_getH(allknots, normalize) # , 1 normalize T/F 1/0
  
  MatpreH <- NULL
  if(preH) {
    if(is.null(MatpreH)) {
      H <- GS_getH(allknots, normalize) 
    } else if(!is.null(MatpreH)) {
      if(MatpreH[1,1] == 1) {
        H <- GS_getH(allknots, normalize) 
      } else {
        H <- MatpreH
      }
    }
  } else if(!preH) {
    H <- GS_getH(allknots, normalize) 
  }
  
  
  out <- (bs %*% H)
  if (calcderiv) out <- (bsderiv %*% H)
  if(intercept) out <- out else out <- out[, 2:ncol(out)]
  
  return(out)
}




#' An internal function to computer H matrix
#' 
#' @param knots A vector
#' @param normalize An integer
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_getH <- function(knots, 
                    normalize) {
  # Step 1: Define constants based on the knots
  Nintk <- length(knots) - 8
  C11 <- 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]))
  C31 <- 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]))
  C21 <- -C11 - C31
  Cp22 <- 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] -
                                                          knots[Nintk + 4]))
  Cp2  <- 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - 
                                                          knots[Nintk + 4]))
  Cp12 <- -Cp22 - Cp2
  
  # Step 2: Build the matrix H depending on Nintk
  if (Nintk == 0) {
    H <- t(matrix(c(3, 0, 
                    2, 1, 
                    1, 2, 
                    0, 3), nrow = 4, byrow = TRUE))
    # print(H)
  } else if (Nintk == 1) {
    H <- matrix(c(-C21 / C11 ,        1, 0, 0, 0,
                  0, -C31 / C21, 1, -Cp22 / Cp12, 0,
                  # 0, C31 / C11, 1, C31 / Cp22, 0,
                  0, 0, 0, 1, -Cp12 / Cp2),
                nrow = 3, byrow = TRUE)
  } else {
    H1 <- rbind(matrix(1, nrow = 1, ncol = 3), 
                c(0, 1, -C21 / C31), 
                matrix(0, nrow = Nintk - 2, ncol = 3), 
                matrix(0, nrow = 2, ncol = 3))
    
    H2 <- rbind(matrix(0, nrow = 2, ncol = Nintk - 2),
                diag(Nintk - 2), 
                matrix(0, nrow = 2, ncol = Nintk - 2))
    
    H3 <- rbind(matrix(0, nrow = 2, ncol = 3),
                matrix(0, nrow = Nintk - 2, ncol = 3),
                c(-Cp12 / Cp22, 1, 0),
                matrix(1, nrow = 1, ncol = 3))
    H <- cbind(H1 , H2 , H3 )
  }
  
  if(normalize) {
    sumH <- rowSums(H)
    notzero <- which(sumH != 0)
    H[notzero, ] <- H[notzero, ]  / sumH[notzero]
  }
 
  return(t(H))
}


#' An internal function to construct natural cubic spline basis matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_nsp_call <- function(x, 
                        knots, 
                        bknots, 
                        degree,
                        intercept, 
                        derivs, 
                        centerval, 
                        normalize, 
                        preH, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL) {
  
  if(derivs > 1) {
    stop("Second and higher order derivatives are not supported yet")
  } else {
    calcderiv <- derivs
  }
  
  if(intercept) {
    if(centerval != 0) {
      stop("centerval should be '0' when intercept = TRUE")
      centerval <- 0
    }
  }
  
  if(length(knots) > 0) {
    if(bknots[1] > knots[1]) {
      stop("Lower boundary is greater than lower internal knot")
    } else if(bknots[2] < knots[length(knots)]) {
      stop("Upper Boundary knot is less than upper internal knot")
    }
  }
  
  df     <- length(knots) + 1 + intercept
  out <- GS_nsp_nsk(x, 
                    knots, 
                    bknots, 
                    degree = degree, 
                    intercept = intercept, 
                    calcderiv = calcderiv, 
                    normalize = normalize,
                    preH = preH)
  
  # if no internal knot, the out is a vector, not matrix, convert it to matrix
  # but if no internal knot but intercept TRUE, then it is already a matrix
  if(length(knots) == 0) {
   if(!intercept) out <- matrix(out, length(out), 1)
  }
  
  # Centering
  if (centerval != 0) {
    cenout <- GS_nsp_nsk(centerval, 
                         knots, bknots, 
                         degree = degree, 
                         intercept = intercept, 
                         calcderiv = calcderiv, 
                         normalize = normalize,
                         preH = preH)
    
    if(!is.matrix(cenout)) cenout <- matrix(cenout, nrow = 1) 
    # if(length(knots) == 0) {
    #   if(!intercept) cenout <- matrix(cenout, length(cenout), 1)
    # }
    if (!calcderiv) {
      if(intercept) {
        for (i in 2:ncol(cenout)) {
          out[,i] <- out[,i] - cenout[,i]
        }
      } else if(!intercept) {
        for (i in 1:ncol(cenout)) {
          out[,i] <- out[,i] - cenout[,i]
        }
      }
    } else if (calcderiv) {
      if(length(knots) == 0) {
        if(!intercept) out <- matrix(out, length(out), 1)
      }
    }
  } # if (centerval != 0) {

  return(out)
}



#' An internal function to construct a variant of natural cubic spline basis
#' matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_nsk_call <- function(x, 
                        knots, 
                        bknots, 
                        degree, 
                        intercept, 
                        derivs, 
                        centerval, 
                        normalize, 
                        preH, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL) {
  
  if(derivs > 1) {
    stop("Second and higher order derivatives are not supported yet")
  } else {
    calcderiv <- derivs
  }
  
  temp <- c(bknots, knots)
  if(length(knots) == 0) {
    kx   <- c(bknots[1], bknots[2])
  } else {
    kx   <- c(bknots[1], knots, bknots[2])
  }
  
 
  if (!calcderiv) {
    basis <- GS_nsp_call(x, 
                         knots = temp[-(1:2)], 
                         bknots = bknots, 
                         degree = degree,
                         intercept = intercept, 
                         derivs = derivs, 
                         centerval = centerval, 
                         normalize = normalize,
                         preH = preH)
    
    kbasis <- GS_nsp_call(kx, 
                          knots=knots, 
                          bknots = bknots, 
                          degree = degree,
                          intercept = intercept, 
                          derivs = derivs, 
                          centerval = centerval, 
                          normalize = normalize,
                          preH = preH)
    
    if(intercept) {
      out <- basis %*% solve(kbasis)
    } else {
      out <- (cbind(1, basis) %*% solve(cbind(1, kbasis)))[, -1]
    }
  } else if (calcderiv) {
    basis <- GS_nsp_call(x, 
                         knots = temp[-(1:2)], 
                         bknots = bknots, 
                         degree = degree,
                         intercept = intercept, 
                         derivs = derivs, 
                         centerval = centerval, 
                         normalize = normalize,
                         preH = preH)
    
    kbasis <- GS_nsp_call(kx, 
                          knots=knots, 
                          bknots = bknots, 
                          degree = degree,
                          intercept = intercept, 
                          derivs = 0, 
                          centerval = centerval, 
                          normalize = normalize,
                          preH = preH)
    
    if(intercept) {
      out <- basis %*% solve(kbasis)
    } else {
      out <- (cbind(1, basis) %*% solve(cbind(1, kbasis)))[, -1]
    }
  }
  
  if(length(knots) == 0) {
    if(!intercept) out <- matrix(out, length(out), 1)
  }
  
  return(out)
}



#' An internal function to construct a variant of natural cubic spline basis
#' matrix (truncated power basis based restricted cubic spline)
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default). This is useful when using
#'   \code{rcs_matrix} for creating design matrix for derivatives such as
#'   \code{deriv = 1} and \code{deriv = 2} where first (\code{deriv = 1}) or,
#'   the first and second (\code{deriv = 2}) columns are automatically set as
#'   \code{'0'}.
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval Ignored
#' @param normalize Ignored
#' @param preH Ignored
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. It is defined as \code{nk - 1}. An integer. The
#'   \code{df} defined only for R functions, not Stan.
#' @param inclx A logical indicating whether to include \code{x} as first term.
#' @param fullknots combined internal knots and boundary knots
#' @param fullknots.only A logical indicating whether to return knots only.
#' @param type An integer to indicate type, \code{0} for \code{'ordinary'} and
#'   \code{1} for \code{'integral'}. This string to integer for easy to code
#'   \code{'rcs_matrix_stan'}. Note that this \code{'type'} argument is not used
#'   even in \code{'rcs_matrix'}
#' @param verbose A logical (default \code{FALSE}) to indicate if infornation
#'   need to be displayed.
#'   
#' @inherit Hmisc::rcspline.eval params
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_rcs_call <- function(x, 
                        knots = NULL, 
                        bknots = NULL,
                        degree = NULL,
                        intercept = NULL, 
                        derivs = 0, 
                        centerval = NULL, 
                        normalize = NULL, 
                        preH = NULL, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL, 
                        inclx = TRUE, 
                        fullknots = NULL,
                        fullknots.only = FALSE,
                        type = 0, 
                        norm = 2, 
                        rpm = NULL, 
                        pc = FALSE,
                        fractied = 0.05,
                        verbose = FALSE) {
  
  
  if(!is.null(bknots)) {
    if(length(bknots) != 2) {
      stop2c("The 'bknots' should be length 2")
    }
  }
  
  # if everything is missing, then set df
  # next preferance should be fullknots
  if(is.null(knots) & is.null(bknots) & is.null(fullknots) & is.null(df)) {
    if(verbose) message2c("'df' is set as 3 for the rcs model")
    df <- 3
  } else if(!is.null(fullknots)) {
    if(!is.null(knots) & !is.null(bknots)) {
      stop2c("When specifying 'fullknots', both 'knots' and 'bknots' 
           should be NULL")
    }
    fullknots <- fullknots
  } else if(is.null(fullknots)) {
    if(!is.null(knots) & !is.null(bknots)) {
      fullknots <- c(bknots[1], knots, bknots[2])
    } else if(is.null(knots) | is.null(bknots)) {
      if(is.null(df)) {
        stop2c("please specify knots and bknots, fullknots, or df")
      }
    }
  } else {
    stop2c("Specify 'df' or 'knots' 
    Note that knots can be specified by using 'fullknots', 
    or else via 'knots' and 'bknots'")
  }
  
 
  if(is.null(intercept)) {
    intercept <- FALSE
  }
  
  
  if(!is.null(fullknots) | !is.null(knots)) {
    df <- NULL
    if(verbose) stop2c("Both 'df' and 'knots' were specified, 
                     keeping 'knots' and 'setting df = NULL'") 
  }
  
  # Here onward, knots = fullknots
  # knots <- fullknots
  
  if(is.null(df) & is.null(fullknots)) {
    stop("please specify either df or knots, not both")
  }
  if(!is.null(df) & !is.null(fullknots)) {
    stop("please specify either df or knots, not both")
  }
  
  if(is.null(df) & !is.null(fullknots)) {
    nk <- length(fullknots)
  } else if(!is.null(df) & is.null(fullknots)) {
    nk <- df + 1
  } else {
    stop()
  }
  

  # sourced from Hmisc::rcspline.eval
 if (!length(fullknots)) {
   if(is.data.frame(x) | data.table::is.data.table(x)) {
     stop("The argument 'x' must be a numeric vector but found data frame")
   }
    xx <- x[!is.na(x)]
    n <- length(xx)
    if (n < 6) 
      stop("fullknots not specified, and < 6 non-missing observations")
    if (nk < 3) 
      stop("nk must be >= 3")
    xu <- sort(unique(xx))
    nxu <- length(xu)
    if ((nxu - 2) <= nk) {
      warning(sprintf("%s fullknots requested with %s unique values of x.
                      fullknots set to %s interior values.", 
                      nk, nxu, nxu - 2))
      fullknots <- xu[-c(1, length(xu))]
    }
    else {
      outer <- if (nk > 3) 
        0.05
      else 0.1
      if (nk > 6) 
        outer <- 0.025
      fullknots <- numeric(nk)
      overrideFirst <- overrideLast <- FALSE
      nke <- nk
      firstknot <- lastknot <- numeric(0)
      if (fractied > 0 && fractied < 1) {
        f <- table(xx)/n
        if (max(f[-c(1, length(f))]) < fractied) {
          if (f[1] >= fractied) {
            firstknot <- min(xx[xx > min(xx)])
            xx <- xx[xx > firstknot]
            nke <- nke - 1
            overrideFirst <- TRUE
          }
          if (f[length(f)] >= fractied) {
            lastknot <- max(xx[xx < max(xx)])
            xx <- xx[xx < lastknot]
            nke <- nke - 1
            overrideLast <- TRUE
          }
        }
      }
      if (nke == 1) 
        fullknots <- median(xx)
      else {
        if (nxu <= nke) 
          fullknots <- xu
        else {
          p <- if (nke == 2) 
            seq(0.5, 1 - outer, length = nke)
          else seq(outer, 1 - outer, length = nke)
          fullknots <- quantile(xx, p)
          if (length(unique(fullknots)) < min(nke, 3)) {
            fullknots <- quantile(xx, seq(outer, 1 - outer, 
                                          length = 2 * nke))
            if (length(firstknot) && length(unique(fullknots)) < 
                3) {
              midval <- if (length(firstknot) && length(lastknot)) 
                (firstknot + lastknot)/2
              else median(xx)
              fullknots <- 
                sort(c(firstknot, 
                       midval, 
                       if (length(lastknot)) lastknot else quantile(xx, 1-outer)))
            }
            if ((nu <- length(unique(fullknots))) < 3) {
              cat("Fewer than 3 unique fullknots  Frequency table of variable:\n")
              print(table(x))
              stop()
            }
            warning(paste("could not obtain", nke, 
                          "interior fullknots with default algorithm.\n", 
                          "Used alternate algorithm to obtain", nu, 
                          "fullknots"))
          }
        }
        if (length(xx) < 100) {
          xx <- sort(xx)
          if (!overrideFirst) 
            fullknots[1] <- xx[5]
          if (!overrideLast) 
            fullknots[nke] <- xx[length(xx) - 4]
        }
      }
      fullknots <- c(firstknot, fullknots, lastknot)
    }
  }
  fullknots <- sort(unique(fullknots))
  
  
  ####
  if(!is.null(bknots)) {
    fullknots_org                <- fullknots
    fullknots[1]                 <- bknots[1]
    fullknots[length(fullknots)] <- bknots[2]
    if(verbose) {
      message2c("The 'boundary knots' replaced by the provided 'bknots'. 
                Thus, the 'knots' used now are ", 
                fullknots, ", instead of the ", fullknots_org)
    }
  }
  #######
  
  nk <- length(fullknots)
  if (nk < 3) {
    cat("fewer than 3 unique fullknots  Frequency table of variable:\n")
    print(table(x))
    stop()
  }
  if (fullknots.only) 
    return(fullknots)
  if (length(rpm)) 
    x[is.na(x)] <- rpm
  # end of borrow from Hmisc::rcspline.eval
  ##
  
  X <- x
  N <- length(X)
  nk <- length(fullknots)
  
  basis_evals <- matrix(0, N, nk-1)
  
  if(inclx) basis_evals[,1] = X
  
  if(derivs == 0) basis_evals[,1] = X;
  if(derivs == 1) basis_evals[,1] = 1;
  if(derivs == 2) basis_evals[,1] = 0;
  
  Xx <- matrix(0, N, nk)
  km1 = nk - 1;
  j = 1;
  knot1   <- fullknots[1     ]
  knotnk  <- fullknots[nk    ]
  knotnk1 <- fullknots[nk - 1]
  kd <-     (knotnk - knot1) ^ (2)
  
  for(ia in 1:N) {
    for(ja in 1:nk) {
      Xx[ia,ja] = ifelse(X[ia] - fullknots[ja] > 0, X[ia] - fullknots[ja], 0)
    }
  }
  
  if(derivs == 0) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] = 
        (Xx[,j]^3-(Xx[,km1]^3)*(fullknots[nk]-fullknots[j])/
           (fullknots[nk]-fullknots[km1]) + (Xx[,nk]^3)*(fullknots[km1]-fullknots[j])/
           (fullknots[nk]-fullknots[km1])) / (fullknots[nk]-fullknots[1])^2;
      j = j + 1;
    }
  }
  
  if(derivs == 1) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] =
        (3*Xx[,j]^2) * (1/((fullknots[nk]-fullknots[1])^2))  - 
        (3*Xx[,km1]^2)*(fullknots[nk]-fullknots[j]) /
        ((fullknots[nk]-fullknots[km1]) * (fullknots[nk]-fullknots[1])^2) + 
        (3*Xx[,nk]^2)*(fullknots[km1]-fullknots[j])/
        ((fullknots[nk]-fullknots[km1]) * (fullknots[nk]-fullknots[1])^2) ;
      j = j + 1;
    }
  }
  
  if(derivs == 2) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] =
        (6*Xx[,j]^1) * (1/((fullknots[nk]-fullknots[1])^2))  - 
        (6*Xx[,km1]^1)*(fullknots[nk]-fullknots[j]) /
        ((fullknots[nk]-fullknots[km1]) * (fullknots[nk]-fullknots[1])^2) + 
        (6*Xx[,nk]^1)*(fullknots[km1]-fullknots[j])/
        ((fullknots[nk]-fullknots[km1]) * (fullknots[nk]-fullknots[1])^2) ;
      j = j + 1;
    }
  }
  
  if(!inclx) basis_evals <- basis_evals[,-1,drop=FALSE]
  
  if(intercept) {
    if(derivs == 0) {
      mat_intercept <- matrix(1, nrow(basis_evals), 1)
      basis_evals <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept column added. Please use ~0 + formula")
    }
    if(derivs == 1) {
      mat_intercept <- matrix(0, nrow(basis_evals), 1)
      basis_evals <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept set to '0' for deriv = 1")
    }
    if(derivs == 2) {
      mat_intercept <- matrix(0, nrow(basis_evals), 2)
      basis_evals   <- basis_evals[, -1, drop = FALSE]
      basis_evals   <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept and first term (x) set to '0' for derivs = 2")
    }
  } # if(intercept) {
  
  attr(basis_evals, "knots") <- fullknots
 
  return(basis_evals)
} # GS_rcs_call



#' An internal function to perform elementwise multipication and rowsum 
#' This is used in R inverse wide in QR rcs _d1
#' @param matb A numeric matrix of betas with dim(N, N)
#' @param matx A numeric matrix of predictor with dim(N, M)
#' @param rowsum A logical whether to return rowsum
#' 
#' @return A matrix or a vector
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
prodrowsum_beta_Rinv_wide <- function(matb, 
                                      matx, 
                                      rowsum = TRUE) {
  matrix_X <- matx
  matrix_Y <- matb
  num_cols <- ncol(matrix_X)
  num_rows <- nrow(matrix_X)
  result_matrix_loop <- matrix(NA, nrow = num_rows, ncol = num_cols)
  for (j in 1:num_cols) { # Outer loop iterates through columns
    for (i in 1:num_rows) { # Inner loop iterates through rows
      # Perform element-wise multiplication for the current (row, column)
      result_matrix_loop[i, j] <- matrix_X[i, j] * matrix_Y[i, j]
    }
  }
  if(rowsum) {
    out <- Matrix::rowSums(result_matrix_loop)
  } else {
    out <- result_matrix_loop
  }
  return(out)
}




##########################################################################

#' An internal function to construct a variant of natural cubic spline basis
#' matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#' @param Boundary.knots Ignored
#' @param periodic Ignored
#' @param integral Ignored
#' @param warn.outside Ignored
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_bsp_call <- function(x, 
                        knots = NULL,
                        bknots = NULL,
                        degree, 
                        intercept, 
                        derivs, 
                        centerval, 
                        normalize, 
                        preH, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL, 
                        Boundary.knots = NULL,
                        fullknots = NULL,
                        periodic = FALSE, 
                        integral = FALSE, 
                        warn.outside = getOption("splines2.warn.outside", TRUE)
) {
  

  
  if(!is.null(fullknots)) {
    if(!is.null(knots)) stop("'knots' must be NULL if specified 'fullknots'")
    if(!is.null(bknots)) stop("'bknots' must be NULL if specified 'fullknots'")
    get_knots_bknots <- split_fullknots_knots_bknots(fullknots)
    knots  <- get_knots_bknots[['knots']]
    bknots <- get_knots_bknots[['bknots']]
  } else if(is.null(fullknots)) { 
    if(is.null(df)) {
      # for bsp msp and isp, knots = NULL  is a valid approach
      if(is.null(knots) | is.null(bknots)) {
       # stop("Both knots and bknots must be specified when fullknots/df NULL")
      }
      knots  <- knots
      bknots <- bknots
    }
  }
  
 
  out <- splines2::bsp(x = x, 
                       df = df, 
                       knots = knots, 
                       degree = degree, 
                       intercept = intercept, 
                       Boundary.knots = bknots, 
                       periodic = periodic, 
                       derivs = derivs, 
                       integral = integral, 
                       warn.outside = warn.outside)
  
  return(out)
}


#' An internal function to construct a variant of natural cubic spline basis
#' matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#' @param Boundary.knots Ignored
#' @param periodic Ignored
#' @param integral Ignored
#' @param warn.outside Ignored
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_msp_call <- function(x, 
                        knots = NULL,
                        bknots = NULL,
                        degree, 
                        intercept, 
                        derivs, 
                        centerval, 
                        normalize, 
                        preH, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL, 
                        Boundary.knots = NULL,
                        fullknots = NULL,
                        periodic = FALSE, 
                        integral = FALSE, 
                        warn.outside = getOption("splines2.warn.outside", TRUE)
                        ) {
  
  if(!is.null(fullknots)) {
    if(!is.null(knots)) stop("'knots' must be NULL if specified 'fullknots'")
    if(!is.null(bknots)) stop("'bknots' must be NULL if specified 'fullknots'")
    get_knots_bknots <- split_fullknots_knots_bknots(fullknots)
    knots  <- get_knots_bknots[['knots']]
    bknots <- get_knots_bknots[['bknots']]
  } else if(is.null(fullknots)) { 
    if(is.null(df)) {
      # for bsp msp and isp, knots = NULL  is a valid approach
      if(is.null(knots) | is.null(bknots)) {
        # stop("Both knots and bknots must be specified when fullknots/df NULL")
      }
      knots  <- knots
      bknots <- bknots
    }
  }
  
  out <- splines2::msp(x = x, 
                       df = df, 
                       knots = knots, 
                       degree = degree, 
                       intercept = intercept, 
                       Boundary.knots = bknots, 
                       periodic = periodic, 
                       derivs = derivs, 
                       integral = integral, 
                       warn.outside = warn.outside)
  
  return(out)
}



#' An internal function to construct a variant of natural cubic spline basis
#' matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#' @param Boundary.knots Ignored
#' @param periodic Ignored
#' @param integral Ignored
#' @param warn.outside Ignored
#' 
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
GS_isp_call <- function(x, 
                        knots = NULL,
                        bknots = NULL,
                        degree, 
                        intercept, 
                        derivs, 
                        centerval, 
                        normalize, 
                        preH, 
                        sfirst = FALSE, 
                        sparse = FALSE,
                        df = NULL, 
                        Boundary.knots = NULL,
                        fullknots = NULL,
                        # periodic = FALSE, 
                        integral = FALSE, 
                        warn.outside = getOption("splines2.warn.outside", TRUE)
) {
  
  
  if(!is.null(fullknots)) {
    if(!is.null(knots)) stop("'knots' must be NULL if specified 'fullknots'")
    if(!is.null(bknots)) stop("'bknots' must be NULL if specified 'fullknots'")
    get_knots_bknots <- split_fullknots_knots_bknots(fullknots)
    knots  <- get_knots_bknots[['knots']]
    bknots <- get_knots_bknots[['bknots']]
  } else if(is.null(fullknots)) { 
    if(is.null(df)) {
      # for bsp msp and isp, knots = NULL  is a valid approach
      if(is.null(knots) | is.null(bknots)) {
        # stop("Both knots and bknots must be specified when fullknots/df NULL")
      }
      knots  <- knots
      bknots <- bknots
    }
  }
  
  out <- splines2::isp(x = x, 
                       df = df, 
                       knots = knots, 
                       degree = degree, 
                       intercept = intercept, 
                       Boundary.knots = bknots, 
                       # periodic = periodic, 
                       derivs = derivs, 
                       integral = integral, 
                       warn.outside = warn.outside)
  
  return(out)
}





#' An internal function to construct a variant of natural cubic spline basis
#' matrix
#'
#' @param x A numeric vector for which basis matrix to be constructed
#' @param knots A vector specifying the internal knots
#' @param bknots A vector specifying the boundary knots
#' @param intercept An integer to indicate whether to compute complete basis
#'   along with intercept (\code{intercept = 1}) or to exclude intercept from
#'   the basis (\code{intercept = 0}, default).
#' @param derivs An integer to indicate whether to compute complete basis matrix
#'   (\code{derivs = 0}, default) or its first derivative (\code{derivs = 1})
#' @param centerval A real number to offset the intercept.
#' @param normalize An integer to indicate whether to normalize the basis matrix
#'   (\code{normalize = 1}) or not (\code{normalize = 0}, default).
#' @param preH A logical (as.integer()) indicating whether to use pre computed 
#' H matrix
#' @param sfirst Ignored
#' @param sparse Ignored
#' @param method  Internal. method 'bs'/'rs' used in knots_selection. Ignored
#'   here
#' @param smat Internal
#' @param degree Internal
#' @param nk Internal
#' @param inclx Internal
#' @param knots.only Internal
#' @param type Internal
#' @param norm Internal
#' @param rpm Internal
#' @param pc Internal
#' @param fractied Internal
#' @param bkrange Internal
#' @param fix_bknots Internal
#' @param xoffset Internal
#' @param bound Internal
#' @param userdata Internal
#' @param verbose Internal
#' @param df An integer. The \code{df} defined only for R functions, not Stan. 
#'
#' @return A matrix
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
get_knost_from_df <- function(x, 
                              df,
                              method = NULL,
                              smat,
                              knots = NULL, 
                              bknots = NULL, 
                              degree = 3, 
                              intercept = FALSE, 
                              derivs = 0, 
                              centerval = FALSE, 
                              normalize = FALSE, 
                              preH = FALSE, 
                              sfirst = FALSE, 
                              sparse = FALSE,
                              nk = NULL,
                              inclx = TRUE,
                              knots.only = TRUE,
                              type = "ordinary",
                              norm = 2,
                              rpm = NULL,
                              pc = FALSE,
                              fractied = 0.05,
                              bkrange = FALSE,
                              fix_bknots = TRUE,
                              xoffset = 0,
                              bound = NULL,
                              userdata = NULL,
                              verbose = FALSE) {
  
  
  if(is.null(userdata)) {
    if(is.data.frame(x) | data.table::is.data.table(x)) {
      stop("The argument 'x' must be a numeric vector but found data frame")
    }
  } else if(!is.null(userdata)) {
    if(is.character(x)) {
      x <- userdata[[x]]
    } else if(!is.character(x)) {
      x <- deparse(x)
      x <- userdata[[x]]
    }
  } # if(is.null(userdata)) { else if(!is.null(userdata)) {
  
  if(is.character(df))        df     <- str2lang(df)
  
  if(is.null(bound)) {
    bound <- 0
  } else {
    if(is.character(bound))  bound     <- str2lang(bound)
  }
  
  
  apply_bound <-  TRUE
  if(is.null(bknots)) {
    if(as.logical(bkrange)) {
      bknots <- range(x)
      apply_bound <-  TRUE
    }
  } else {
    bknots <- bknots
    apply_bound <-  FALSE
  }
  
  bknots_org <- bknots
  
  if(is.null(nk)) {
    nk <- df + 1
  } else {
    if(is.character(nk))  nk     <- str2lang(nk)
  }
  
  
  # bsp_msp_isp_args                <- list()
  # bsp_msp_isp_args[['x']]         <-  x
  # bsp_msp_isp_args[['knots']]     <-  knots
  # bsp_msp_isp_args[['bknots']]    <-  bknots
  # bsp_msp_isp_args[['df']]        <-  df
  # bsp_msp_isp_args[['degree']]    <-  degree
  # bsp_msp_isp_args[['intercept']] <-  intercept
  # bsp_msp_isp_args[['derivs']]    <-  derivs
  # bsp_msp_isp_args[['centerval']] <-  centerval
  # bsp_msp_isp_args[['normalize']] <-  normalize
  # bsp_msp_isp_args[['preH']]      <-  preH
  # bsp_msp_isp_args[['sfirst']]    <-  sfirst
  # bsp_msp_isp_args[['sparse']]    <-  sparse
  
  
  direct_call_args                        <- list()
  direct_call_args[['x']]                 <-  x
  direct_call_args[['knots']]             <-  knots
  direct_call_args[['Boundary.knots']]    <-  bknots
  direct_call_args[['df']]                <-  df
  direct_call_args[['intercept']]         <-  intercept
  
  direct_call_args_nsp_nsk                <- direct_call_args
  direct_call_args_nsp_nsk[['derivs']]    <-  derivs
  
  
  direct_call_args_msp_bsp                <- direct_call_args_nsp_nsk
  direct_call_args_msp_bsp[['degree']]    <-  degree
  
  
  rcspline_eval_args                 <- list()
  rcspline_eval_args[['x']]          <-  x
  rcspline_eval_args[['nk']]         <-  nk
  rcspline_eval_args[['inclx']]      <-  inclx
  rcspline_eval_args[['knots.only']] <-  knots.only
  rcspline_eval_args[['type']]       <-  type
  rcspline_eval_args[['norm']]       <-  norm
  rcspline_eval_args[['rpm']]        <-  rpm
  rcspline_eval_args[['pc']]         <-  pc
  rcspline_eval_args[['fractied']]   <-  fractied
  
  
  if(smat == 'ns') {
    temp_mat_s <- do.call(splines::ns, direct_call_args)
  } else if(smat == 'nsk') {
    temp_mat_s <- do.call(splines2::nsk, direct_call_args_nsp_nsk)
  } else if(smat == 'nsp') {
    temp_mat_s <- do.call(splines2::nsp, direct_call_args_nsp_nsk)
  } else if(smat == 'bsp') {
    temp_mat_s <- do.call(splines2::bsp, direct_call_args_msp_bsp)
  } else if(smat == 'msp') {
    temp_mat_s <- do.call(splines2::msp, direct_call_args_msp_bsp)
  } else if(smat == 'isp') {
    temp_mat_s <- do.call(splines2::isp, direct_call_args_msp_bsp)
  } else if(smat == 'rcs') {
    temp_mat_s <- do.call(Hmisc::rcspline.eval, rcspline_eval_args)
  }
  
  
  if(smat == 'rcs') {
    knots <- temp_mat_s
    temp_mat_s_bknots <- bknots
    if(apply_bound) {
      temp_mat_s_bknots <- apply_bknots_bounds(temp_mat_s_bknots, bound)
    }
    if(as.logical(fix_bknots)) {
      knots[1]               <- temp_mat_s_bknots[1]
      knots[length(knots)]   <- temp_mat_s_bknots[2]
    } else {
      knots <- knots
    }
  } else {
    temp_mat_s_knots <- attr(temp_mat_s, "knots") 
    if(is_emptyx(temp_mat_s_knots)) {
      temp_mat_s_knots <- NULL 
    }
    temp_mat_s_bknots <- attr(temp_mat_s, "Boundary.knots") 
    if(fix_bknots) {
      temp_mat_s_bknots <- bknots_org
    }
    if(apply_bound) {
      temp_mat_s_bknots <- apply_bknots_bounds(temp_mat_s_bknots, bound)
    }
    knots <- c(temp_mat_s_bknots[1], temp_mat_s_knots, temp_mat_s_bknots[2])
  }
  
  knots <- knots - ept(xoffset)
  
  return(knots)
}



