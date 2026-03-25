

###############################################################################
######################## bsp, msp & isp FUNCTIONS #############################
###############################################################################

# GS_bsp_tuple_intermediate_stan
# GS_bsp_intermediate_stan
# GS_bsp_call_stan
# 
# GS_isp_tuple_stan
# GS_isp_stan
# GS_isp_call_stan
# 
# GS_msp_tuple_stan
# GS_msp_stan
# GS_msp_call_stan


################################################################################
# GS_bsp_tuple_intermediate_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_tuple_intermediate_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bsp_intermediate_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
 if (calcderiv) {
    if (intercept) {
    return bsderiv;
  } else {
    return bsderiv[, 2:cols(bsderiv)];
  }
 } else {
    if (intercept) {
    return bs;
  } else {
    return bs[, 2:cols(bs)];
  }
 }  
}
"
  return(set_stan_str)
} # GS_bsp_tuple_intermediate_stan_R_str



################################################################################
# GS_bsp_intermediate_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_intermediate_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bsp_intermediate_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  } 
 if (calcderiv) {
    if (intercept) {
    return bsderiv;
  } else {
    return bsderiv[, 2:cols(bsderiv)];
  }
 } else {
    if (intercept) {
    return bs;
  } else {
    return bs[, 2:cols(bs)];
  }
 }  
}
"
  return(set_stan_str)
} # GS_bsp_intermediate_stan_R_str



################################################################################
# GS_bsp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_bsp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  int degree_in = degree + 0;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;  
  // int ncolselect = Nk + intercept - 1;  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;  
  int ncolselect = Nintk + degree_in + 1 + intercept - 1;  
  // generate basis of degree 0
  if(degree == 0) {
      matrix[N, Nintervals] out_degree0 = rep_matrix(0, N, Nintervals);
      if (calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            out_degree0[i, j] = 0.0;
          }
        }
      }
      if (!calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            real denom =  allknots[j+1] - allknots[j];
          //  real inv_denom = 1 / denom;
            real inv_denom = 1;
            // for last cloumn to match splines::msp, need this condition which 
            // differe by x[i] <= allknots[j + 1]
            if(j < Nintervals) {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] < allknots[j + 1]) ? inv_denom : 0;
            } else {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] <= allknots[j + 1]) ? inv_denom : 0;
            }
          }
        }
      }
     
     if(intercept) {
        return out_degree0[, 1:cols(out_degree0)];
     } else {
        return out_degree0[, 2:cols(out_degree0)];
     } 
  } // end if(degree == 0) {
  
  matrix[N, ncolselect] out = GS_bsp_intermediate_stan(x, knots, bknots, fullknots, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_bsp_intermediate_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }

  return out;
}     
"
  return(set_stan_str)
} # GS_bsp_call_stan_R_str




################################################################################
# GS_isp_tuple_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_tuple_stan_R_str <- function() {
  set_stan_str <- "
/////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
 // Main I-spline function
matrix GS_isp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
  int Psize = degree + 1 +  size(knots);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  // Calculate rowsum of consecutive bs columns  
  if (calcderiv) {
   for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bsderiv[r, (i+0):(Psize+0)]);
     }
    }
  } else {
    for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bs[r, (i+0):(Psize+0)]);
     }
    }
  }
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_isp_tuple_stan_R_str





################################################################################
# GS_isp_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
 // Main I-spline function
matrix GS_isp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  }
  int Psize = degree + 1 +  size(knots);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  // Calculate rowsum of consecutive bs columns  
  if (calcderiv) {
   for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bsderiv[r, (i+0):(Psize+0)]);
     }
    }
  } else {
    for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bs[r, (i+0):(Psize+0)]);
     }
    }
  }
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_isp_stan_R_str




################################################################################
# GS_isp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_isp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  // here in _call_stan, degree is renamed as degree_in because int degree part of fun
  int degree_in = degree + 1;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;
  
  // int ncolselect = Nk + intercept - 1;
  int Psize = degree_in + 1 + size(knots) + intercept - 1;
  int ncolselect = Psize;
  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;
  matrix[N, ncolselect] out = GS_isp_stan(x, knots, bknots, fullknots, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_isp_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }  
  return out[, 2:cols(out)];
}
"
  return(set_stan_str)
} # GS_isp_call_stan_R_str



################################################################################
# GS_msp_tuple_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_tuple_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  // Main M-spline function
matrix GS_msp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
  int Psize = cols(bs);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  if (calcderiv) {
   for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bsderiv[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  } else {
    for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bs[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  }        
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_msp_tuple_stan_R_str



################################################################################
# GS_msp_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  // Main M-spline function
matrix GS_msp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  } 
  int Psize = cols(bs);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  if (calcderiv) {
   for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bsderiv[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  } else {
    for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bs[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  }        
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_msp_stan_R_str



################################################################################
# GS_msp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_msp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  // here in _call_stan, degree is renamed as degree_in because int degree part of fun
  int degree_in = degree + 0;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;
  
  vector[2*(degree+1) + Nintk] fullknots_mspl;
  fullknots_mspl = append_row(append_row(rep_vector(bknotsx[1], (degree+1)), knotsx), rep_vector(bknotsx[2], (degree+1)));
  
  // int ncolselect = Nk + intercept - 1;
  int Psize = degree_in + 1 + size(knots) + intercept - 1;
  int ncolselect = Psize;
  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;
  
  // generate basis of degree 0
  if(degree == 0) {
      matrix[N, Nintervals] out_degree0 = rep_matrix(0, N, Nintervals);
      if (calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            out_degree0[i, j] = 0.0;
          }
        }
      }
      if (!calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            real denom =  allknots[j+1] - allknots[j];
            // for last cloumn to match splines::msp, need this condition which 
            // differe by x[i] <= allknots[j + 1]
            if(j < Nintervals) {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] < allknots[j + 1]) ? 1/denom : 0;
            } else {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] <= allknots[j + 1]) ? 1/denom : 0;
            }
          }
        }
      }
     
     if(intercept) {
        return out_degree0[, 1:cols(out_degree0)];
     } else {
        return out_degree0[, 2:cols(out_degree0)];
     } 
  } // end if(degree == 0) {
  
  
  matrix[N, ncolselect] out = GS_msp_stan(x, knots, bknots, fullknots_mspl, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_msp_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots_mspl, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }
  
  return out;
} 
"
  return(set_stan_str)
} # GS_msp_call_stan_R_str




