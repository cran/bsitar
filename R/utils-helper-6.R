

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' @description
#' Note that now all three functions nsp, nsk and rcs are set using
#' 'utils-helper-7' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it
#' yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsp_nsk_helper_stan_R_str <- function() {
  set_stan_str <- "
 ////////////////////////////////////////////////////////////////////////
  // Find elements in vector x that equal/below/above real number y
 /////////////////////////////////////////////////////////////////////////
  int num_matches(vector x, real y, real z) {
    int n = 0;
    for (i in 1:rows(x))
      if(z == 0) {
        if (x[i] == y) n += 1;
      } else if(z == 1.0) {
        if (x[i] > y) n += 1;
      } else if(z == -1.0) {
        if (x[i] < y) n += 1;
      }
    return n;
  }
  /////////////////////////////////////////////////////////////////////////    
  // Find the indexes of the elements in the vector x that equal real number y
  /////////////////////////////////////////////////////////////////////////
  array[] int which_equal(vector x, real y, real z) {
    array[num_matches(x, y, z)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if(z == 0) {
        if (x[i] == y) {
          match_positions[pos] = i;
          pos += 1;
        }
      } else if(z == 1.0) {
        if (x[i] > y) {
          match_positions[pos] = i;
          pos += 1;
        }
      } else if(z == -1.0) {
        if (x[i] < y) {
          match_positions[pos] = i;
          pos += 1;
        }
      }
    }
    return match_positions;
  }

  /////////////////////////////////////////////////////////////////////////    
  // Find the indexes of the elements in the vector x that equal real number y
  // more compact version, check it - for now, commented out
  /////////////////////////////////////////////////////////////////////////
  /*
  array[] int which_equal(vector x, real y, real z) {
    int N = num_elements(x);
    int n = num_matches(x, y, z);
    array[n] int match_positions;
    int pos = 1;
    if (z == 0) {
      for (i in 1:N) if (x[i] == y) { match_positions[pos] = i; pos += 1; }
    } else if (z == 1.0) {
      for (i in 1:N) if (x[i] > y) { match_positions[pos] = i; pos += 1; }
    } else if (z == -1.0) {
      for (i in 1:N) if (x[i] < y) { match_positions[pos] = i; pos += 1; }
    }
    return match_positions;
  }
  */
  /////////////////////////////////////////////////////////////////////////
  // Repeating input vector K times
  /////////////////////////////////////////////////////////////////////////
  vector repeat_vector(vector input, int K) {
    int N = rows(input);
    vector[N*K] repvec; // stack N-vector K times
    // assign i-th value of input to i+(k-1)*N -th value of repvec
    for (k in 1:K) {
      for (i in 1:N) {
        repvec[i+(k-1)*N] = input[i]; 
      }
    }
    return repvec;
  }
"
  return(set_stan_str)
} # GS_nsp_nsk_helper_stan_R_str



################################################################################
# GS_bsp_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bsp_stan(vector x, vector knots, vector bknots, 
                    vector fullknots, vector allknots, 
                    int N, int degree, int ord, 
                    int Nintk, int Nk, int Nintervals,
                    int intercept, int calcderiv) {    
    matrix[N, Nintervals] M1 = rep_matrix(0, N, Nintervals);    
    for(i in 1:N) {
      for (j in 1:Nintervals) {
        M1[i,j] = (allknots[j] <= x[i]  && x[i] < allknots[j + 1] ? 1 : 0);
      }
    }
    array[num_matches(x, bknots[2], 0)] int lastknot_index = which_equal(x, bknots[2], 0); 
    if (size(lastknot_index) > 0) {
      for(i in lastknot_index) {
        M1[i, (Nintervals - degree):Nintervals] = rep_row_vector(1, (Nintervals-(Nintervals - degree))+1);
      }
    }
    matrix[N, Nintervals - 0] M2x = rep_matrix(0, N, Nintervals);
    matrix[N, Nintervals - 0] M1x = M1;
    vector[N] C1;
    vector[N] C2;
    for (p in 1:degree) {
      for (i in 1:(Nintervals - p)) {
        if (allknots[i + p] == allknots[i]) {
          C1 = rep_vector(0.0, N);
        } else {
          C1 = (x - allknots[i]) / (allknots[i + p] - allknots[i]);
        }
        if (allknots[i + p + 1] == allknots[i + 1]) {
          C2 = rep_vector(0.0, N);
        } else {
          C2 = (allknots[i + p + 1] - x) / (allknots[i + p + 1] - allknots[i + 1]);
        }
        M2x[, i] = C1 .* M1x[, i] + C2 .* M1x[, i + 1];
      }
      if(p != degree) M1x = M2x;
    }
    matrix[N, Nintervals - degree]    M2 = M2x[, 1:(Nintervals - degree)];
    matrix[N, (Nintervals - degree)+1] M1xx = M1x[, 1:(Nintervals - degree)+1];
    if (calcderiv) {
      matrix[N, Nintervals - degree] deriv;
      for (i in 1:(Nintervals - degree)) {
        real dC1 = (allknots[i + degree] == allknots[i]) ? 0 : degree / (allknots[i + degree] - allknots[i]);
        real dC2 = (allknots[i + degree + 1] == allknots[i + 1]) ? 0 : degree / (allknots[i + degree + 1] - allknots[i + 1]);
        deriv[, i] = dC1 * M1xx[, i] - dC2 * M1xx[, i + 1];
      }
      return deriv;
    } else {
      return M2;
    }
  } // end matrix GS_bsp_stan  
"
  return(set_stan_str)
} # GS_bsp_stan_R_str



################################################################################
# GS_bsp_tuple_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_tuple_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Function to compute B-splines and their derivatives - Tuple version
  /////////////////////////////////////////////////////////////////////////
  tuple(matrix, matrix) GS_bsp_tuple_stan(vector x, vector knots, vector bknots, 
                    vector fullknots, vector allknots, 
                    int N, int degree, int ord, 
                    int Nintk, int Nk, int Nintervals,
                    int intercept, int calcderiv) {
    matrix[N, Nintervals] M1 = rep_matrix(0, N, Nintervals);
    for (i in 1:N)
      for (j in 1:Nintervals)
        M1[i, j] = (allknots[j] <= x[i] && x[i] < allknots[j + 1]) ? 1 : 0;
  
    array[num_matches(x, bknots[2], 0)] int lastknot_index = which_equal(x, bknots[2], 0);
    if (size(lastknot_index) > 0)
      for (ii in 1:size(lastknot_index))
        for (jj in (Nintervals - degree):(Nintervals))
          M1[lastknot_index[ii], jj] = 1;
  
    matrix[N, Nintervals] M2x;
    matrix[N, Nintervals] M1x = M1;
    for (p in 1:degree) {
      for (i in 1:(Nintervals - p)) {
        vector[N] C1 = (allknots[i + p] == allknots[i]) ? rep_vector(0.0, N) : (x - allknots[i]) / (allknots[i + p] - allknots[i]);
        vector[N] C2 = (allknots[i + p + 1] == allknots[i + 1]) ? rep_vector(0.0, N) : (allknots[i + p + 1] - x) / (allknots[i + p + 1] - allknots[i + 1]);
        M2x[, i] = C1 .* M1x[, i] + C2 .* M1x[, i + 1];
      }
      if (p != degree) M1x = M2x;
    }
    matrix[N, Nintervals - degree] M2 = M2x[, 1:(Nintervals - degree)];
    matrix[N, Nintervals - degree + 1] M1xx = M1x[, 1:(Nintervals - degree + 1)];
    matrix[N, Nintervals - degree] deriv;
    for (i in 1:(Nintervals - degree)) {
      real dC1 = (allknots[i + degree] == allknots[i]) ? 0 : degree / (allknots[i + degree] - allknots[i]);
      real dC2 = (allknots[i + degree + 1] == allknots[i + 1]) ? 0 : degree / (allknots[i + degree + 1] - allknots[i + 1]);
      deriv[, i] = dC1 * M1xx[, i] - dC2 * M1xx[, i + 1];
    }
    return (M2, deriv);
  }
"
  return(set_stan_str)
} # GS_bsp_tuple_stan_R_str




################################################################################
# GS_nsp_nsk_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsp_nsk_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Main function for B-splines and their derivatives
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsp_nsk_stan(vector x, vector knots, vector bknots, 
                    vector fullknots, vector allknots, 
                    int N, int degree, int ord, 
                    int Nintk, int Nk, int Nintervals,
                    int intercept, int calcderiv, int normalize,
                    int preH) {                    
    matrix[num_elements(x), Nintervals - degree] bs;
    matrix[num_elements(x), Nintervals - degree] bsderiv;
    if (calcderiv) {
      bs      = GS_bsp_stan(x, knots, bknots, 
                          fullknots, allknots,  
                          N, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 0);
      bsderiv = GS_bsp_stan(x, knots, bknots, 
                          fullknots, allknots,  
                          N, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 1);
    } else {
      bs      = GS_bsp_stan(x, knots, bknots, 
                          fullknots, allknots,
                          N, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 0);
    }
    array[num_matches(x, bknots[1], -1)] int x_below_boundary  = which_equal(x, bknots[1], -1); 
    array[num_matches(x, bknots[2],  1)] int x_above_boundary  = which_equal(x, bknots[2],  1); 
    int size_x_below_boundary = size(x_below_boundary);
    int size_x_above_boundary = size(x_above_boundary);
    if (size_x_below_boundary > 0 || size_x_above_boundary > 0) {
      matrix[2, Nintervals - degree] bs_bknots;
      matrix[2, Nintervals - degree] bsderiv_bknots;
      bs_bknots      = GS_bsp_stan(bknots, knots, bknots, 
                          fullknots, allknots, 
                          2, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 0);
      bsderiv_bknots = GS_bsp_stan(bknots, knots, bknots, 
                          fullknots, allknots,
                          2, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 1);      
    if (size_x_below_boundary > 0) {
      array[size_x_below_boundary] int xselect = x_below_boundary;
      int Nxselect = size_x_below_boundary;
      matrix[size_x_below_boundary, (Nintervals - degree)] below_azx1;
      matrix[size_x_below_boundary, (Nintervals - degree)] below_azx2;
      matrix[size_x_below_boundary, (Nintervals - degree)] below_azx3;      
      below_azx1 = transpose(to_matrix((repeat_vector(transpose(bs_bknots[1,]), 
                    size_x_below_boundary)), Nintervals - degree, 
                    size_x_below_boundary));
      below_azx2 = transpose(to_matrix((repeat_vector(transpose(bsderiv_bknots[1,]), 
                    size_x_below_boundary)), Nintervals - degree, 
                    size_x_below_boundary));
      below_azx3 = (to_matrix((repeat_vector(((x[xselect] - bknots[1])), 
                    Nintervals - degree)), size_x_below_boundary, 
                    Nintervals - degree));
        bs[xselect, ] = below_azx1 + (below_azx2 .* below_azx3);
        if (calcderiv) {
          for(i in xselect) {
            bsderiv[i, 1: Nintervals - degree] = (bsderiv_bknots[1,]);
          }
        }
      } // if (size_x_below_boundary > 0)      
      if (size_x_above_boundary > 0) {
        array[size_x_above_boundary] int xselect = x_above_boundary;
        int Nxselect = size_x_above_boundary; 
        matrix[size_x_above_boundary, (Nintervals - degree)] below_azx1;
        matrix[size_x_above_boundary, (Nintervals - degree)] below_azx2;
        matrix[size_x_above_boundary, (Nintervals - degree)] below_azx3;
        below_azx1 = transpose(to_matrix((repeat_vector(transpose(bs_bknots[2,]), 
                                                        size_x_above_boundary)), Nintervals - degree, 
                                         size_x_above_boundary));
        below_azx2 = transpose(to_matrix((repeat_vector(transpose(bsderiv_bknots[2,]), 
                                                        size_x_above_boundary)), Nintervals - degree, 
                                         size_x_above_boundary));
        below_azx3 = (to_matrix((repeat_vector(((x[xselect] - bknots[2])), 
                                               Nintervals - degree)), size_x_above_boundary, 
                                Nintervals - degree));        
        bs[xselect, ] = below_azx1 + (below_azx2 .* below_azx3);        
        if (calcderiv) {
            for(i in xselect) {
              bsderiv[i,1: Nintervals - degree ] = (bsderiv_bknots[2,]);
           }         
        }
      } // if (size_x_above_boundary > 0)    
    } // if (size_x_below_boundary > 0 || size_x_above_boundary > 0) {
    matrix[Nk+2, Nk] H;
    // Now MatpreH is just a place holder that will be replaced my matrix
    if(preH) {
     H = MatpreH;
    } else {
     H = GS_getH_stan(allknots, normalize);
    }
    matrix[num_elements(x), Nk] out;
    int ncolselect = Nk+intercept-1;
    matrix[num_elements(x), ncolselect] result;
    if (calcderiv) {
     out = (bsderiv * H);
      if(intercept) {
       result = out;
      } else {
       result = out[, 2:cols(out)];
      }
    } else {
     out = bs * H;
     if(intercept) {
       result = out;
      } else {
       result = out[, 2:cols(out)];
      }
    }
    return result;
  } // end matrix GS_nsp_nsk_stan  
"
  return(set_stan_str)
} # GS_bsp_tuple_stan_R_str





################################################################################
# GS_nsp_nsk_stan_fast_1_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsp_nsk_stan_fast_1_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Main function for B-splines and their derivatives
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsp_nsk_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  /*
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  }
  */
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
  // Precompute boundary indices only if needed
  int n_below = num_matches(x, bknots[1], -1);
  int n_above = num_matches(x, bknots[2], 1);
  array[n_below] int x_below_boundary;
  array[n_above] int x_above_boundary;
  if (n_below > 0) x_below_boundary = which_equal(x, bknots[1], -1);
  if (n_above > 0) x_above_boundary = which_equal(x, bknots[2], 1);

  if (n_below > 0 || n_above > 0) {
    tuple(matrix[2, n_basis], matrix[2, n_basis]) my_tuple =
      GS_bsp_tuple_stan(bknots, knots, bknots, fullknots, allknots, 2, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
    matrix[2, n_basis] bs_bknots = my_tuple.1;
    matrix[2, n_basis] bsderiv_bknots = my_tuple.2;

    if (n_below > 0) {
      vector[n_below] x_below = x[x_below_boundary];
      row_vector[n_basis] bsknot1 = bs_bknots[1, ];
      row_vector[n_basis] bsknot1_deriv = bsderiv_bknots[1, ];
      matrix[n_below, n_basis] rep_bsknot1 = rep_matrix(bsknot1, n_below);
      matrix[n_below, n_basis] rep_bsknot1_deriv = rep_matrix(bsknot1_deriv, n_below);
      matrix[n_below, n_basis] xdiff = rep_matrix(to_row_vector(x_below - bknots[1]), n_basis)';
      bs[x_below_boundary, ] = rep_bsknot1 + rep_bsknot1_deriv .* xdiff;
      if (calcderiv) {
        for (i in 1:n_below)
          bsderiv[x_below_boundary[i], ] = bsknot1_deriv;
      }
    }
    if (n_above > 0) {
      vector[n_above] x_above = x[x_above_boundary];
      row_vector[n_basis] bsknot2 = bs_bknots[2, ];
      row_vector[n_basis] bsknot2_deriv = bsderiv_bknots[2, ];
      matrix[n_above, n_basis] rep_bsknot2 = rep_matrix(bsknot2, n_above);
      matrix[n_above, n_basis] rep_bsknot2_deriv = rep_matrix(bsknot2_deriv, n_above);
      matrix[n_above, n_basis] xdiff = rep_matrix(to_row_vector(x_above - bknots[2]), n_basis)';
      bs[x_above_boundary, ] = rep_bsknot2 + rep_bsknot2_deriv .* xdiff;
      if (calcderiv) {
        for (i in 1:n_above)
          bsderiv[x_above_boundary[i], ] = bsknot2_deriv;
      }
    }
  }
    // Precompute H matrix only once
    matrix[Nk + 2, Nk] H;
    if (preH) {
      H = MatpreH;
    } else {
      H = GS_getH_stan(allknots, normalize);
    }
  
    int ncolselect = Nk + intercept - 1;
    matrix[N, Nk] out = calcderiv ? (bsderiv * H) : (bs * H);
  
    // Avoid unnecessary allocation
    if (intercept) {
      return out;
    } else {
      return out[, 2:cols(out)];
    }
  }
"
  return(set_stan_str)
} # GS_nsp_nsk_stan_fast_1_R_str


################################################################################
# GS_nsp_nsk_stan_fast_2_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsp_nsk_stan_fast_2_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Main function for B-splines and their derivatives
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsp_nsk_stan(
    vector x, vector knots, vector bknots, 
    vector fullknots, vector allknots, 
    int N, int degree, int ord, 
    int Nintk, int Nk, int Nintervals,
    int intercept, int calcderiv, int normalize,
    int preH) {
    int n_basis = Nintervals - degree;
    /*
    matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
    matrix[N, n_basis] bsderiv;
    if (calcderiv) {
      bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
    }
    */
    tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
    GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
    matrix[N, n_basis] bs = my_tuple_main.1;
    matrix[N, n_basis] bsderiv;
    if (calcderiv) {
      bsderiv = my_tuple_main.2;
    }
  
    int n_below = num_matches(x, bknots[1], -1);
    int n_above = num_matches(x, bknots[2], 1);
    array[n_below] int x_below_boundary;
    array[n_above] int x_above_boundary;
    if (n_below > 0) x_below_boundary = which_equal(x, bknots[1], -1);
    if (n_above > 0) x_above_boundary = which_equal(x, bknots[2], 1);
  
    if (n_below > 0 || n_above > 0) {
      tuple(matrix[2, n_basis], matrix[2, n_basis]) my_tuple =
        GS_bsp_tuple_stan(bknots, knots, bknots, fullknots, allknots, 2, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
      matrix[2, n_basis] bs_bknots = my_tuple.1;
      matrix[2, n_basis] bsderiv_bknots = my_tuple.2;
  
      if (n_below > 0) {
        vector[n_below] x_below = x[x_below_boundary];
        row_vector[n_basis] bsknot1 = bs_bknots[1, ];
        row_vector[n_basis] bsknot1_deriv = bsderiv_bknots[1, ];
        matrix[n_below, n_basis] rep_bsknot1 = rep_matrix(bsknot1, n_below);
        matrix[n_below, n_basis] rep_bsknot1_deriv = rep_matrix(bsknot1_deriv, n_below);
        matrix[n_below, n_basis] xdiff = rep_matrix(to_row_vector(x_below - bknots[1]), n_basis)';
        bs[x_below_boundary, ] = rep_bsknot1 + rep_bsknot1_deriv .* xdiff;
        if (calcderiv) {
          for (i in 1:n_below)
            bsderiv[x_below_boundary[i], ] = bsknot1_deriv;
        }
      }
      if (n_above > 0) {
        vector[n_above] x_above = x[x_above_boundary];
        row_vector[n_basis] bsknot2 = bs_bknots[2, ];
        row_vector[n_basis] bsknot2_deriv = bsderiv_bknots[2, ];
        matrix[n_above, n_basis] rep_bsknot2 = rep_matrix(bsknot2, n_above);
        matrix[n_above, n_basis] rep_bsknot2_deriv = rep_matrix(bsknot2_deriv, n_above);
        matrix[n_above, n_basis] xdiff = rep_matrix(to_row_vector(x_above - bknots[2]), n_basis)';
        bs[x_above_boundary, ] = rep_bsknot2 + rep_bsknot2_deriv .* xdiff;
        if (calcderiv) {
          for (i in 1:n_above)
            bsderiv[x_above_boundary[i], ] = bsknot2_deriv;
        }
      }
    }
    matrix[Nk + 2, Nk] H;
    if (preH) {
      H = MatpreH;
    } else {
      H = GS_getH_stan(allknots, normalize);
    }
  
    int ncolselect = Nk + intercept - 1;
    matrix[N, Nk] out = calcderiv ? (bsderiv * H) : (bs * H);
  
    if (intercept) {
      return out;
    } else {
      return out[, 2:cols(out)];
    }
  }
"
  return(set_stan_str)
} # GS_nsp_nsk_stan_fast_2_R_str





################################################################################
# GS_nsp_call_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsp_call_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // call GS_nsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N        = num_elements(x);
  vector[num_elements(knotsx)+2] fullknots;
  vector[num_elements(fullknots)-2] knots;
  vector[2] bknots;
  fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  knots = segment(fullknots, 2, num_elements(fullknots)-2);
  bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
   int Nintk   = num_elements(knots);
   // int degree  = 3;
   int ord     = degree + 1;
   int Nk      = Nintk + 2;
   int df      = Nintk + 1 + intercept;
   int calcderiv;
   if(derivs > 1) {
    } else {
      calcderiv = derivs;
   }
   if(intercept) {
    if(centerval != 0) {
      // centerval = 0;
    }
  }
  if(Nintk > 0) {
    if(bknots[1] > knots[1]) {
      } else if(bknots[2] < knots[Nintk]) {
    }
  }
  int ncolselect = Nk+intercept-1; 
  matrix[N, ncolselect] out;
  vector[Nintk + 2*ord] allknots;
   allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), 
                          rep_vector(bknots[2], ord));    
   int Nintervals  = num_elements(allknots) - 1;
    out = GS_nsp_nsk_stan(x, knots, bknots, fullknots, allknots, 
                        N, degree, ord, 
                        Nintk, Nk, Nintervals, 
                        intercept, calcderiv, normalize, preH);                        
    if (centerval != 0) {
    matrix[1, df] cenout;
    cenout = GS_nsp_nsk_stan(centerval+rep_vector(0.0, 1), knots, bknots, fullknots, allknots,
                        1, degree, ord, 
                        Nintk, Nk, Nintervals, 
                        intercept, calcderiv, normalize, preH);      
      if (!calcderiv) {
        if(intercept) {
          for (i in 2:cols(cenout)) {
            out[,i] = out[,i] - cenout[1,i];
          }
        } else if(!intercept) {
          for (i in 1:cols(cenout)) {
            out[,i] = out[,i] - cenout[1,i];
          }
        }
      } 
    } // if (centerval != 0) {
  return(out);
  } // end matrix GS_nsp_call_stan
"
  return(set_stan_str)
} # GS_nsp_call_stan_R_str



################################################################################
# GS_nsk_call_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_nsk_call_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // call GS_nsk_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsk_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {       
  int N        = num_elements(x);
  vector[num_elements(knotsx)+2] fullknots;
  vector[num_elements(fullknots)-2] knots;
  vector[2] bknots;
  fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  knots = segment(fullknots, 2, num_elements(fullknots)-2);
  bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
   int Nintk   = num_elements(knots);
   // int degree  = 3;
   int ord     = degree + 1;
   int Nk      = Nintk + 2;
   int df      = Nintk + 1 + intercept;
   int calcderiv;
   if(derivs > 1) {
    } else {
      calcderiv = derivs;
   }
   if(intercept) {
    if(centerval != 0) {
      // centerval = 0;
    }
  }
  if(Nintk > 0) {
    if(bknots[1] > knots[1]) {
      } else if(bknots[2] < knots[Nintk]) {
    }
  }  
  int ncolselect = Nk+intercept-1;
  matrix[N, ncolselect] basis;
  matrix[Nk, ncolselect] kbasis;
  matrix[N, ncolselect+1] kout;
  matrix[N, ncolselect] out;
   if(!calcderiv) {
     basis  = GS_nsp_call_stan(x,         knots, bknots, degree, intercept, derivs, centerval, normalize, preH);
     kbasis = GS_nsp_call_stan(fullknots, knots, bknots, degree, intercept, derivs, centerval, normalize, preH);
   } else {
     basis  = GS_nsp_call_stan(x,         knots, bknots, degree, intercept, derivs, centerval, normalize, preH);
     kbasis = GS_nsp_call_stan(fullknots, knots, bknots, degree, intercept, 0, centerval, normalize, preH);
   }
    if(intercept) {
       out = basis * inverse(kbasis);
    } else {
       kout = append_col(rep_vector(1, N), basis) * 
              inverse(append_col(rep_vector(1, Nk), kbasis));
       out = kout[, 2:cols(kout)];
    } 
    return(out);
  } // end matrix GS_nsk_call_stan
"
  return(set_stan_str)
} # GS_nsk_call_stan_R_str





################################################################################
# GS_nsk_call_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_preH_stan_R_str <- function() {
  set_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Function to calculate H matrix (as defined in GS_getH_stan)
  /////////////////////////////////////////////////////////////////////////
  matrix GS_getH_stan(vector knots, int normalize) {
  int Nintk = num_elements(knots) - 8;
  real C11 = 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]));
  real C31 = 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]));
  real C21 = -C11 - C31;
  real Cp22 = 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] - knots[Nintk + 4]));
  real Cp2 = 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - knots[Nintk + 4]));
  real Cp12 = -Cp22 - Cp2;

  if (Nintk == 0) {
    matrix[2, 4] H = to_matrix(transpose([3, 0, 2, 1, 1, 2, 0, 3]), 2, 4);
    if (normalize) {
      vector[2] sumH = H * rep_vector(1.0, 4);
      for (i in 1:2) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else if (Nintk == 1) {
    matrix[3, 5] H = transpose(to_matrix([
      -C21 / C11, 1, 0, 0, 0, 
      0, -C31 / C21, 1, -Cp22 / Cp12, 0, 
      0, 0, 0, 1, -Cp12 / Cp2
    ], 5, 3));
    if (normalize) {
      vector[3] sumH = H * rep_vector(1.0, 5);
      for (i in 1:3) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else if (Nintk == 2) {
    matrix[4, 3] H1 = append_row(
      append_row(
        append_row(rep_matrix(1, 1, 3), [0, 1, -C21 / C31]),
        rep_matrix(0, 0, 3)
      ),
      rep_matrix(0, 2, 3)
    );
    matrix[4, 3] H3 = append_row(
      append_row(
        append_row(rep_matrix(0, 2, 3), rep_matrix(0, 0, 3)),
        [-Cp12 / Cp22, 1, 0]
      ),
      rep_matrix(1, 1, 3)
    );
    matrix[4, 6] H = append_col(H1, H3);
    if (normalize) {
      vector[4] sumH = H * rep_vector(1.0, 6);
      for (i in 1:4) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else {
    int K = Nintk - 2;
    int N = Nintk + 2;
    matrix[N, 3] H1 = append_row(
      append_row(
        append_row(rep_matrix(1, 1, 3), [0, 1, -C21 / C31]),
        rep_matrix(0, K, 3)
      ),
      rep_matrix(0, 2, 3)
    );
    matrix[N, K] H2 = append_row(
      append_row(
        rep_matrix(0, 2, K), diag_matrix(rep_vector(1, K))
      ),
      rep_matrix(0, 2, K)
    );
    matrix[N, 3] H3 = append_row(
      append_row(
        append_row(rep_matrix(0, 2, 3), rep_matrix(0, K, 3)),
        [-Cp12 / Cp22, 1, 0]
      ),
      rep_matrix(1, 1, 3)
    );
    matrix[N, N + 2] H = append_col(append_col(H1, H2), H3);
    if (normalize) {
      vector[N] sumH = H * rep_vector(1.0, N + 2);
      for (i in 1:N) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
    }
  } // end matrix GS_getH_stan
"
  return(set_stan_str)
} # GS_preH_stan_R_str



################################################################################
# GS_auxillary_stan_R_str - not used anywhere
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_auxillary_stan_R_str <- function() {
  set_stan_str <- "
  ////////////////////////////////////////////////////////////////////////
  // Find elements in vector x that equal/below/above real number y
 /////////////////////////////////////////////////////////////////////////
  // Efficient vector matching/counting
  int num_matches(vector x, real y, real z) {
    int n = 0;
    int N = num_elements(x);
    if (z == 0) {
      for (i in 1:N) if (x[i] == y) n += 1;
    } else if (z == 1.0) {
      for (i in 1:N) if (x[i] > y) n += 1;
    } else if (z == -1.0) {
      for (i in 1:N) if (x[i] < y) n += 1;
    }
    return n;
  }
  /////////////////////////////////////////////////////////////////////////
   // Check if vectors VecSame
   /////////////////////////////////////////////////////////////////////////
  int VecSame(vector A, vector B) {
    int N = num_elements(A);
    for (i in 1:N) {
      if (A[i] != B[i]) return 0;
    }
    return 1;
  }
  /////////////////////////////////////////////////////////////////////////
   // Check if matrices VecSame
   /////////////////////////////////////////////////////////////////////////
  int MatSame(matrix A, matrix B) {
    int R = rows(A);
    int C = cols(A);
    for (i in 1:R) {
      for (j in 1:C) {
        if (A[i, j] != B[i, j]) return 0;
      }
    }
    return 1;
  }

  /////////////////////////////////////////////////////////////////////////
  // Repeating input matrix K times
  /////////////////////////////////////////////////////////////////////////
  matrix repeat_matrix(matrix input, int K) {
    int N = rows(input);
    int M = cols(input);
    matrix[N*K,M] repmat; // stack N*M matrix K times
    // assign i-th row of input to i+(k-1)*N -th row of repmat
    for (k in 1:K) {
      for (i in 1:N) {
        repmat[i+(k-1)*N] = input[i]; 
      }
    }
    return repmat;
  }  
  "
  return(set_stan_str)
} # GS_auxillary_stan_R_str


###############################################################################
###############################################################################
############################# RCS FUNCTION ####################################
###############################################################################
###############################################################################


################################################################################
# GS_rcs_call_stan_R_str
################################################################################

#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_rcs_call_stan_R_str <- function() {
  set_stan_str <- "
 ////////////////////////////////////////////////////////////////////////
  // Function to calculate truncated power basis based rcs matrix
  /////////////////////////////////////////////////////////////////////////
  /**
   * Calculates the Restricted Cubic Spline (RCS) basis matrix.
   *
   * This Stan function provides the core logic for generating a Restricted Cubic Spline
   * basis matrix, given an input vector `x` and pre-defined fullknots.
   *
   * IMPORTANT: Knot selection logic (e.g., quantiles, `Hmisc::rcspline.eval`)
   * must be performed OUTSIDE of Stan, in your data preparation script (R, Python, etc.).
   * The inclx is now defined within the functiin and not as an argument. 
   * The arguments are exact same as GS_nsp/nsk stan
   *
   * Args:
   * x: A vector of observations for which to calculate the basis.
   * fullknots: A vector of ordered unique fullknots. Must have at least 3 fullknots.
   * derivs: Integer representing the derivative order (0 for value, 1 for 1st derivs, 2 for 2nd derivs).
   * intercept: Integer (3). Just a placeholder.
   * intercept: Integer (0 or 1). If 1, adds an intercept column.
   * inclx: Integer (0 or 1). If 1, includes the 'x' column (for derivs = 0).
   *
   * Returns:
   * A matrix where rows correspond to observations in `x` and columns
   * correspond to the basis functions.
   */
   
  matrix GS_rcs_call_stan(vector x, 
                        vector knotsx, 
                        vector bknotsx,
                        int degree, 
                        int intercept, 
                        int derivs, 
                        real centerval,
                        int normalize,
                        int preH) {
                        
    int inclx = 1;                  
    int N = num_elements(x);
    
    vector[num_elements(knotsx)+2] fullknots;
    vector[num_elements(fullknots)-2] knots;
    vector[2] bknots;
    fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
    knots = segment(fullknots, 2, num_elements(fullknots)-2);
    bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
    
    int nk = num_elements(fullknots);
   
    matrix[N, nk - 1] basis_evals;
    real knot1 = fullknots[1];
    real knot_nk = fullknots[nk];
    real knot_nk_minus_1 = fullknots[nk - 1];
    real denom_common = (knot_nk - knot1)^2;
    real denom_spline_diff = (knot_nk - knot_nk_minus_1);

    matrix[N, nk] Xx_pos; 
    for (n_idx in 1:N) {
      for (k_idx in 1:nk) {
        Xx_pos[n_idx, k_idx] = fmax(0.0, x[n_idx] - fullknots[k_idx]);
      }
    }

    if (derivs == 0) {
      basis_evals[:, 1] = x;
    } else if (derivs == 1) {
      basis_evals[:, 1] = rep_vector(1.0, N);
    } else if (derivs == 2) {
      basis_evals[:, 1] = rep_vector(0.0, N);
    } else {
      //
    }

    for (j in 1:(nk - 2)) {
      int jp1 = j + 1; 
      if (derivs == 0) {
        basis_evals[:, jp1] =
          (pow(Xx_pos[:, j], 3) -
           pow(Xx_pos[:, nk - 1], 3) * (fullknots[nk] - fullknots[j]) / denom_spline_diff +
           pow(Xx_pos[:, nk], 3) * (fullknots[nk - 1] - fullknots[j]) / denom_spline_diff) / denom_common;
      } else if (derivs == 1) {
        basis_evals[:, jp1] =
          (3.0 * pow(Xx_pos[:, j], 2)) / denom_common -
          (3.0 * pow(Xx_pos[:, nk - 1], 2)) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (3.0 * pow(Xx_pos[:, nk], 2)) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      } else if (derivs == 2) {
        basis_evals[:, jp1] =
          (6.0 * Xx_pos[:, j]) / denom_common -
          (6.0 * Xx_pos[:, nk - 1]) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (6.0 * Xx_pos[:, nk]) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      }
    }

    matrix[N, 0] final_basis_evals_empty = rep_matrix(0.0, N, 0); 

    if (intercept == 1) {
      if (derivs == 0) {
        if (inclx == 1) {
          matrix[N, nk] temp_basis;
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk)] = basis_evals; 
          return temp_basis;
        } else {
          matrix[N, nk-1] temp_basis; 
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)];
          return temp_basis;
        }
      } else if (derivs == 1) {
        if (inclx == 1) {
            matrix[N, nk] temp_basis;
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk)] = basis_evals; 
            return temp_basis;
        } else {
            matrix[N, nk-1] temp_basis; 
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)]; 
            return temp_basis;
        }
      } else if (derivs == 2) {
        if (inclx == 1) {
             matrix[N, nk] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             print(temp_basis);
             print(basis_evals);
             temp_basis[:, 3:(nk - 0)] = basis_evals[:, 2:(nk - 1)]; 
             return temp_basis;
        } else {
             matrix[N, nk -1] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             temp_basis[:, 3:(nk-1)] = basis_evals[:, 3:(nk - 1)]; 
             return temp_basis;
        }
      }
    } else { // intercept == 0
      if (inclx == 0) {
        return basis_evals[:, 2:(nk-1)]; 
      } else {
        return basis_evals;
      }
    }
    // Should not reach here
    return final_basis_evals_empty;
  } // end matrix GS_rcs_call_stan
"
  return(set_stan_str)
} # GS_rcs_call_stan_R_str



