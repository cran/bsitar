

#############################################################
############### GS_gps_parms_stan_str_get ##################
#############################################################

GS_gps_parms_stan_str_get <- function() {
  
  GS_gps_parms_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // sequence_stan_base
  /////////////////////////////////////////////////////////////////////////
  array[] int sequence_stan_base(int start, int end) {
        array[end - start + 1] int seq;
        for (n in 1:num_elements(seq)) {
          seq[n] = (n + start - 1)  ;
        }
        return seq;
      }


  /////////////////////////////////////////////////////////////////////////
  // factorial_stan
  /////////////////////////////////////////////////////////////////////////
  real factorial_stan(int n) {
  real factorial = 1;
      for(i in 1:n) {
      factorial = factorial * i;
      }
    return factorial;
  }
  
  /////////////////////////////////////////////////////////////////////////
  // seq_fun_stan - seq.int(from, to, length.out)
  /////////////////////////////////////////////////////////////////////////
  vector seq_fun_stan(real start, real end, int N_by) {  
     real h;
     vector[N_by]  out;
     h=(end-start) / (N_by-1);
     for (i in 1:N_by) { 
       out[i]=start + (i-1)*h;
     }
    return(out);
  }
  
  /////////////////////////////////////////////////////////////////////////
  // vector_non_nan
  /////////////////////////////////////////////////////////////////////////
  vector vector_non_nan(vector x) {
    vector[size(x)] y;
    int pos = 1;
    for (i in 1:rows(x)) {
      if(!is_nan(x[i])) {
        y[pos] = x[i];
        pos = pos + 1;
      }
    }
    return y;
  }
  
  /////////////////////////////////////////////////////////////////////////
  // size_non_nan
  /////////////////////////////////////////////////////////////////////////
  int size_non_nan(vector x) {
    int n = 0;
    for (i in 1:rows(x)) {
       if(!is_nan(x[i])) {
         n += 1;
       }
    }
    return n;
  }
  
  /////////////////////////////////////////////////////////////////////////
  // index_non_nan
  /////////////////////////////////////////////////////////////////////////
  array[] int index_non_nan(vector x) {
    array[size_non_nan(x)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if(!is_nan(x[i])) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }
  
  
  /////////////////////////////////////////////////////////////////////////
  // size_nan
  /////////////////////////////////////////////////////////////////////////
  int size_nan(vector x) {
    int n = 0;
    for (i in 1:rows(x)) {
       if(is_nan(x[i])) {
         n += 1;
       }
    }
    return n;
  }
  
  
  /////////////////////////////////////////////////////////////////////////
  // index_nan
  /////////////////////////////////////////////////////////////////////////
  array[] int index_nan(vector x) {
    array[size_nan(x)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if(is_nan(x[i])) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }

  /////////////////////////////////////////////////////////////////////////
  // GS_gps_parms_stan
  /////////////////////////////////////////////////////////////////////////
   array[]  matrix  GS_gps_parms_stan(
                   vector nlp_a,
                   vector nlp_b,
                   vector nlp_c,
                   vector nlp_d,
                   matrix SParMat,
                   vector xknots,
                   array[] matrix spline_eval_array,
                   array[] vector xg_array,
                   array[] vector xg_curve_array,
                   int degree,
                   int shift_indicator,
                   int spline_subset_indicator,
                   int spline_precomputed_indicator,
                   int set_spread,
                   int return_indicator,
                   int drawni
                   ) {

    int pieces_dim = num_elements(xknots) - 1;
    int degree_dim = degree + 1;
    int deriv_root = 2;
    real solvex = 0;
    int mark_peak = 1;

    array[2] int parm_deriv_vec = {0,1};
    array[3] int curve_deriv_vec = {0,1,2};
    
    int n_parm_deriv  = num_elements(parm_deriv_vec);
    int n_curve_deriv = num_elements(curve_deriv_vec);
    
    int array_dim = num_elements(nlp_a);
    int array_mat_dim_1;
    int array_mat_dim_2;
    
    int n_piece_id = 0; 
    
    int parm_mat_dim;
    int mark_peak_dim;
    int deriv_mat_dim_1;
    int deriv_mat_dim;
    
    int degree_dim_set_spread = degree_dim * set_spread;
    
    if(drawni == 0) {
      parm_mat_dim = 6;
      mark_peak_dim = 6;
      deriv_mat_dim_1 = pieces_dim * degree_dim_set_spread;
      deriv_mat_dim   = n_curve_deriv + 1+1+1;
    } else {
      parm_mat_dim = 6+1;
      mark_peak_dim = 6+1;
      deriv_mat_dim_1 = pieces_dim * degree_dim_set_spread;
      deriv_mat_dim   = n_curve_deriv + 1+1+1+1;
    }
    
    
    array[array_dim] matrix[degree_dim, pieces_dim] coef_array;
    array[array_dim] matrix[pieces_dim, parm_mat_dim] parm_array;
    array[array_dim] matrix[deriv_mat_dim_1, deriv_mat_dim] deriv_array;
    

    matrix[1, pieces_dim] rrootsmat;
    matrix[degree_dim, pieces_dim] PiecePolyCoef;
    
    matrix[pieces_dim, parm_mat_dim] parm_mat_local;
 
    matrix[degree_dim_set_spread, deriv_mat_dim] deriv_mat_local;
    

    vector[cols(SParMat)] smat;
    
    
    // This for retun out
     if(return_indicator == 1) {
         array_mat_dim_1 = pieces_dim;
         array_mat_dim_2 = parm_mat_dim;
      }
     if(return_indicator == 2) {
         array_mat_dim_1 = pieces_dim*degree_dim_set_spread;
         array_mat_dim_2 = deriv_mat_dim;
     }
    
    vector[degree_dim] xg;
    vector[degree_dim] yg;
    matrix[degree_dim, degree_dim] Xg;
    matrix[degree_dim, degree_dim] A;
    vector[degree_dim] b;
    matrix[degree_dim, degree_dim] U;
    
    ///////////////////////////////////////////////////////////////////
    // Start elements for GS_nsp_call_stan
    ///////////////////////////////////////////////////////////////////
    
    int calcderiv = 0;
    int intercept = 0;
    int derivs = 0; 
    real centerval = 0; 
    int normalize = 1;
    int preH = 0;
    
    vector[num_elements(xknots)] fullknots;
    vector[num_elements(fullknots)-2] knots;
    vector[2] bknots;
    fullknots = xknots; // append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
    knots = segment(fullknots, 2, num_elements(fullknots)-2);
    bknots[1] = fullknots[1]; 
    bknots[2] = fullknots[num_elements(fullknots)]; 
    int Nintk   = num_elements(knots);
    int ord     = degree + 1;
    int Nk      = Nintk + 2;
    int df      = Nintk + 1 + intercept;
    
    int ncolselect = Nk+intercept-1; 
    matrix[degree_dim, ncolselect] spline_eval;
    
    
    ///////////////////////////////////////////////////////////////////
    // End elements for GS_nsp_call_stan
    ///////////////////////////////////////////////////////////////////
    
    int rowni_i = 0;
    for (rowni in 1:array_dim) {
      if (spline_subset_indicator == 1) {
        smat = SParMat[1, ]';
      } else {
        smat = SParMat[rowni, ]';
      }
      real par_a = nlp_a[rowni];
      real par_b = nlp_b[rowni];
      real par_c = nlp_c[rowni];
      real par_d = nlp_d[rowni];

      for (i in 1:pieces_dim) {
      rowni_i = rowni_i + 1;
        
        
        if(spline_precomputed_indicator == 0) {
          for(j in 1:(degree + 1)){
            xg[j] = xknots[i] + (xknots[i+1] - xknots[i]) * (j-1.0)/degree;
          }
        } else if(spline_precomputed_indicator == 1) {
           xg = xg_array[i];
        }
        
        
        if(spline_precomputed_indicator == 0) {
          spline_eval = GS_nsp_call_stan(xg, knots, bknots, 
                      intercept, derivs, centerval, normalize, preH);
        } else if(spline_precomputed_indicator == 1) {
          spline_eval = spline_eval_array[i];
        }
        
         // print(\"xg: \", xg);
         // print(\"xg_array: \", xg_array[i]);
         // print(\"spline_eval: \", spline_eval);
         // print(\"spline_eval_array: \", spline_eval_array[i]);
        
        yg = spline_eval * smat;
        for(j in 1:(degree + 1)){
        real temp_o = ((xg[j]/exp(par_c)) + par_b);
        real temp_x = (shift_indicator == 1) ? temp_o - xknots[i] : temp_o;
        for(k in 1:(degree + 1)){
          Xg[j, k] = pow(temp_x, k - 1);
          }
        }
        A = Xg' * Xg;
        b = Xg' * yg;
        U = cholesky_decompose(A);
        PiecePolyCoef[, i] = mdivide_left(U', mdivide_left_tri_low(U, b));
        vector[degree_dim - deriv_root] pc;
        if (deriv_root > 0) {
          for (j in 1:(degree_dim - deriv_root)) {
            pc[j] = PiecePolyCoef[j + deriv_root, i] * 
                    choose(deriv_root + j - 1, deriv_root) * 
                    factorial_stan(deriv_root);
          }
        } else {
          pc = PiecePolyCoef[, i];
        }
        pc[1] = pc[1] - solvex;
        real rroots;
        rroots = -pc[1] / pc[2];

        if (shift_indicator == 1) {
          rroots = rroots + xknots[i];
        }
        real xi_lower = (xknots[i] / exp(par_c)) + par_b;
        real xi_upper = (xknots[i + 1] / exp(par_c)) + par_b;
        real get_rroots;
        
        
        
        
        if (rroots >= xi_lower && rroots <= xi_upper) {
          get_rroots = rroots;
        } else {
          get_rroots = positive_infinity();
        }
        
        
        // Also replace first and last piece to NA, boundary
        // Seems no need here, later Also replace is sufficient
        if(rowni == 1) {
        //  get_rroots = positive_infinity();
        }
        
        if(rowni == array_dim) {
        //  get_rroots = positive_infinity();
        }
        
        
        // Set to NA for rounded upto 6 decimal places - to match R behaviour
        real xi_lower_round = round(xi_lower * pow(10, 6)) / pow(10, 6);
        real xi_upper_round = round(xi_upper * pow(10, 6)) / pow(10, 6);
        real rroots_round   = round(rroots   * pow(10, 6)) / pow(10, 6);
        if (rroots_round == xi_lower_round || rroots == xi_upper_round) {
         // get_rroots = positive_infinity();
        }
        
   /*
        if(rowni == 12 || rowni == 12) {
         print(\"xi_lower: \", xi_lower);
         print(\"xi_upper: \", xi_upper);
         print(\"rroots: \", rroots);
         print(\"get_rroots: \", get_rroots);
        }
    */    
    
        
        rrootsmat[1, i] = get_rroots;
      }


      if (return_indicator == 0) {
        coef_array[rowni] = PiecePolyCoef;
      }




      if (return_indicator == 1) {
         vector[pieces_dim] newx_all = rrootsmat[1,]';
         vector[pieces_dim] set_piece_id_all;
         array[size_non_nan(newx_all)] int newx_non_nan  = index_non_nan(newx_all); 
         int size_newx_non_nan = size(newx_non_nan);
         vector[size_newx_non_nan] set_piece_id;
         vector[size_newx_non_nan] newx;
          for (i in 1:pieces_dim) {
           set_piece_id_all[i] = i;
          }
          int pos = 0;
          for (i in newx_non_nan) {
           pos = pos + 1;
           set_piece_id[pos] = set_piece_id_all[i];
           newx[pos] = newx_all[i];
          }
          int n_unique_piece_id = num_elements(set_piece_id);
          for (i in 1:n_unique_piece_id) {
            int ii = newx_non_nan[i];
            real xg_newx_o = newx[ii];
            real xg_newx = xg_newx_o - shift_indicator * xknots[ii];
            vector[degree_dim] pc = PiecePolyCoef[, ii];
            int collectd_mat_i = 1;
            for (dxi in parm_deriv_vec) {
              if (dxi >= degree) {
                reject(\"'parm_deriv' should be less than 'degree' i.e., \", degree);
              }
              vector[degree - dxi + 1] pc_i;
              if (dxi == 0) {
                pc_i = pc;
              } else if (dxi > 0) {
                vector[degree - dxi + 1] pc_temp;
                for (j in 1:(degree_dim - dxi)) {
                  pc_temp[j] = PiecePolyCoef[j + dxi, i] * 
                          choose(dxi + j - 1, dxi) * 
                          factorial_stan(dxi);
                }
                pc_i = pc_temp;
              }
              vector[num_elements(pc_i)] collectd_matx;
              matrix[degree + 1 - dxi, degree + 1 - dxi] xg_pow_dxi;
              for(j in 1:(degree + 1 - dxi)){
                real temp_x = xg_newx;
                for(k in 1:(degree + 1 - dxi)){
                  xg_pow_dxi[j, k] = pow(temp_x, k - 1);
                }
              }
              collectd_matx = xg_pow_dxi * pc_i;
              
              if (dxi == 0) {
                collectd_matx = collectd_matx + par_a;
              }
              parm_mat_local[ii, 1] = xg_newx_o;
              for(ij in 1:size(collectd_matx)) {
               parm_mat_local[ii, ij+1+dxi] = collectd_matx[ij];
              }
              // placeholder for peak/trough
              // other wise previous column value carried to 2+0
              parm_mat_local[i, n_parm_deriv + 2+0] = positive_infinity(); 
              parm_mat_local[i, n_parm_deriv + 2+1] = i;
              parm_mat_local[i, n_parm_deriv + 3+1] = rowni;
              collectd_mat_i = collectd_mat_i + 1;
            } // for (dxi in parm_deriv_vec) {
          } // for (i in 1:n_unique_piece_id) {
          
           vector[pieces_dim] ya_sitar_all = parm_mat_local[, n_parm_deriv + 1];
           vector[pieces_dim] ya_sitar;
           for (mark_peaki in 1:pieces_dim) {
            if(is_nan(ya_sitar_all[mark_peaki])) {
              ya_sitar[mark_peaki] = 0;
            } else if(ya_sitar_all[mark_peaki] == positive_infinity()) {
              ya_sitar[mark_peaki] = 0;
            } else if(ya_sitar_all[mark_peaki] == negative_infinity()) {
              ya_sitar[mark_peaki] = 0;
            } else {
              ya_sitar[mark_peaki] = ya_sitar_all[mark_peaki];
            }
           } // for (mark_peaki in 1:pieces_dim) {
           int lag = 1;
           vector[pieces_dim-1] d2_test; // diff(ya_sitar);
           d2_test = ya_sitar[(1+lag):pieces_dim] - ya_sitar[1:(pieces_dim-lag)];
           
           array[num_matches(d2_test, 0, -1)] int peak_id  = which_equal(d2_test, 0, -1); 
           array[num_matches(d2_test, 0,  1)] int trough_id  = which_equal(d2_test, 0,  1); 

          for(peak_idi in peak_id) {
           parm_mat_local[peak_idi,   n_parm_deriv + 2+0] = 1;
          }
          for(trough_idi in trough_id) {
           parm_mat_local[trough_idi, n_parm_deriv + 2+0] = 0;
          }
          
          // Also replace first and last piece to NA, boundary
          int nrowsmat = rows(parm_mat_local);
          parm_mat_local[1,        n_parm_deriv + 2+0] = positive_infinity();
          parm_mat_local[nrowsmat, n_parm_deriv + 2+0] = positive_infinity();
           
          
         if (drawni > 0) {
           for (mark_peaki in 1:pieces_dim) {
             parm_mat_local[mark_peaki, n_parm_deriv + 3+1+1] = drawni;
           }
          }
        parm_array[rowni] = parm_mat_local;
      } // end if (return_indicator == 1) {
      
      
      
      
    //  print(\"ya_sitar_all: \", peak_id);
    //  print(\"ya_sitar: \", trough_id);
    //  print(\"d2_test: \", parm_mat_local);
          
    // print(\"parm_mat_local: \", parm_mat_local);
    //  print(\"rrootsmat: \", rrootsmat);
    //  print(\"xknots[i]: \", xknots[i]);
      
      
      
      if (return_indicator == 2) {
      vector[degree_dim_set_spread] xg_curve; 
      vector[degree_dim_set_spread] xg_curve_transformed_o;
      vector[degree_dim_set_spread] xg_curve_transformed;
      for (i in 1:pieces_dim) {
        if(spline_precomputed_indicator == 0) {
         xg_curve = seq_fun_stan(xknots[i], xknots[i+1], degree_dim_set_spread);
        }
        if(spline_precomputed_indicator == 1) {
         xg_curve = xg_curve_array[i];
        }
        xg_curve_transformed_o = (xg_curve / exp(par_c) + par_b);
        xg_curve_transformed = xg_curve_transformed_o - shift_indicator * xknots[i]; 
        vector[degree + 1] pc = PiecePolyCoef[, i];
        int collectd_mat_i = 1;
            for (dxi in curve_deriv_vec) {
              if (dxi >= degree) {
                reject(\"'parm_deriv' should be less than 'degree' i.e., \", degree);
              }
              vector[degree + 1 - dxi] pc_i;
              if (dxi == 0) {
                pc_i = pc;
              } else if (dxi > 0) {
                vector[degree + 1 - dxi ] pc_temp;
                for (j in 1:(degree_dim - dxi)) {
                  pc_temp[j] = PiecePolyCoef[j + dxi, i] * 
                          choose(dxi + j - 1, dxi) * 
                          factorial_stan(dxi);
                }
                pc_i = pc_temp;
              }
              vector[degree_dim_set_spread] collectd_matx;
              matrix[degree_dim_set_spread, degree + 1 - dxi] xg_pow_dxi;
              for(k in 1:(degree + 1 - dxi)) {
                xg_pow_dxi[, k] = xg_curve_transformed^(k-1);
              }
              collectd_matx = xg_pow_dxi * pc_i;
              if (dxi == 0) {
                collectd_matx = collectd_matx + par_a;
              }
              deriv_mat_local[, 1] = xg_curve_transformed_o;
              deriv_mat_local[, dxi+2] = collectd_matx;
              deriv_mat_local[, n_curve_deriv + 2]   = i      + rep_vector(0, degree_dim_set_spread);
              deriv_mat_local[, n_curve_deriv + 3]   = rowni  + rep_vector(0, degree_dim_set_spread);
              deriv_mat_local[, n_curve_deriv + 3+1] = drawni + rep_vector(0, degree_dim_set_spread);
              array[degree_dim_set_spread] int seq = 
                                                sequence_stan_base(1+degree_dim_set_spread*(i-1), 
                                                degree_dim_set_spread+degree_dim_set_spread*(i-1));
              deriv_array[rowni, seq, ] = deriv_mat_local;
              collectd_mat_i = collectd_mat_i + 1;
            } // end  for (dxi in curve_deriv_vec) {
        } // end for (i in 1:pieces_dim) {
      } // end if (return_indicator == 2) {
      
   
    
  } // end rowni
  
  
  
   // print(\"deriv_array: \", deriv_array);

   array[array_dim] matrix[array_mat_dim_1, array_mat_dim_2] out;
   if (return_indicator == 0) {
      out = coef_array; 
    } else if (return_indicator == 1) {
      out = parm_array; 
    } else if (return_indicator == 2) {
      out = deriv_array; 
    }
    
  return(out);
} // end GS_gps_parms_stan

"
  
  return(GS_gps_parms_stan_str)
} # GS_gps_parms_stan_str_get







support_GS_gps_parms_stan_str_get <- function() {
  
  support_GS_gps_parms_stan_str <- "
  /////////////////////////////////////////////////////////////////////////
  // Function to calculate H matrix (as defined in GS_ns_getH_stan)
  /////////////////////////////////////////////////////////////////////////
  matrix GS_ns_getH_stan(vector knots, int normalize) {
    int Nintk = num_elements(knots) - 8;
    real C11 = 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]));
    real C31 = 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]));
    real C21 = -C11 - C31;
    real Cp22 = 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] - knots[Nintk + 4]));
    real Cp2 = 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - knots[Nintk + 4]));
    real Cp12 = -Cp22 - Cp2;
    matrix[Nintk+1, Nintk+1] dH; // dummy
    if (Nintk == 0) {
      matrix[2, 4] H;
      H = to_matrix(transpose([3, 0, 2, 1, 1, 2, 0, 3]), 2, 4);
      if (normalize) {
        vector[2] sumH = H * rep_vector(1.0, 4);
        for (i in 1:2) {
          if (sumH[i] != 0) {
            H[i, ] = H[i, ] / sumH[i];
          }
        }
      } // if (normalize) 
      return transpose(H);
    } else if (Nintk == 1) {
      matrix[3, 5] H;
      H = transpose(to_matrix([
        -C21 / C11, 1, 0, 0, 0, 
        0, -C31 / C21, 1, -Cp22 / Cp12, 0, 
        0, 0, 0, 1, -Cp12 / Cp2
      ], 5, 3));
      if (normalize) {
        vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
        for (i in 1:(Nintk + 2)) {
          if (sumH[i] != 0) {
            H[i, ] = H[i, ] / sumH[i];
          }
        }
      } // if (normalize) 
      return transpose(H);
    } else if (Nintk == 2) {
      matrix[Nintk + 2, 3] H1;
      matrix[Nintk + 2, 3] H3;
      matrix[Nintk + 2, Nintk + 2 + 2] H;
      H1 = append_row(
            append_row(
              append_row(
                rep_matrix(1, 1, 3), [0, 1, -C21 / C31]
                ),
                rep_matrix(0,  Nintk - 2, 3)
              ),
            rep_matrix(0,  2, 3)
          );
      H3 = append_row(
            append_row(
              append_row(
                rep_matrix(0, 2, 3), rep_matrix(0, Nintk - 2, 3)
                ),
                [-Cp12 / Cp22, 1, 0]
              ),
            rep_matrix(1,  1, 3)
          );
      H = append_col(H1, H3);
      if (normalize) {
        vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
        for (i in 1:(Nintk + 2)) {
          if (sumH[i] != 0) {
            H[i, ] = H[i, ] / sumH[i];
          }
        }
      } // if (normalize) 
      return transpose(H);
    } else if (Nintk > 2) {
      matrix[Nintk + 2, 3] H1;
      matrix[Nintk + 2, 3] H3;
      matrix[Nintk + 2, Nintk - 2] H2;
      matrix[Nintk + 2, Nintk + 2 + 2] H;
      H1 = append_row(
            append_row(
              append_row(
                rep_matrix(1, 1, 3), [0, 1, -C21 / C31]
                ),
                rep_matrix(0,  Nintk - 2, 3)
              ),
            rep_matrix(0,  2, 3)
          );
      H2 = append_row(
            append_row(
              rep_matrix(0, 2, Nintk - 2), diag_matrix(rep_vector(1, Nintk - 2))
              ),
              rep_matrix(0, 2, Nintk - 2)
            );
      H3 = append_row(
            append_row(
              append_row(
                rep_matrix(0, 2, 3), rep_matrix(0, Nintk - 2, 3)
                ),
                [-Cp12 / Cp22, 1, 0]
              ),
            rep_matrix(1,  1, 3)
          );
      H = append_col(append_col(H1, H2), H3);
      if (normalize) {
        vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
        for (i in 1:(Nintk + 2)) {
          if (sumH[i] != 0) {
            H[i, ] = H[i, ] / sumH[i];
          }
        }
      } // if (normalize) 
      return transpose(H);
    } 
    return dH;
  } // end matrix GS_ns_getH_stan
/////////////////////////////////////////////////////////////////////////
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
  /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bs_stan(vector x, vector knots, vector bknots, 
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
  } // end matrix GS_bs_stan  
  /////////////////////////////////////////////////////////////////////////
  // Main function for B-splines and their derivatives
  /////////////////////////////////////////////////////////////////////////
  matrix GS_ns_stan(vector x, vector knots, vector bknots, 
                    vector fullknots, vector allknots, 
                    int N, int degree, int ord, 
                    int Nintk, int Nk, int Nintervals,
                    int intercept, int calcderiv, int normalize,
                    int preH) {                    
    matrix[num_elements(x), Nintervals - degree] bs;
    matrix[num_elements(x), Nintervals - degree] bsderiv;
    if (calcderiv) {
      bs      = GS_bs_stan(x, knots, bknots, 
                          fullknots, allknots,  
                          N, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 0);
      bsderiv = GS_bs_stan(x, knots, bknots, 
                          fullknots, allknots,  
                          N, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 1);
    } else {
      bs      = GS_bs_stan(x, knots, bknots, 
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
      bs_bknots      = GS_bs_stan(bknots, knots, bknots, 
                          fullknots, allknots, 
                          2, degree, ord,
                          Nintk, Nk, Nintervals,
                          intercept, 0);
      bsderiv_bknots = GS_bs_stan(bknots, knots, bknots, 
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
    
    
    
    
    /*
    if(preH) {
     H = GS_ns_getH_pre(Nk+2, Nk);
    } else {
     H = GS_ns_getH_stan(allknots, normalize);
    }
    */
    
    
    
    H = GS_ns_getH_stan(allknots, normalize);
    
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
  } // end matrix GS_ns_stan  
  /////////////////////////////////////////////////////////////////////////
  // call GS_nsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_nsp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N        = num_elements(x);
  vector[num_elements(knotsx)+2] fullknots;
  vector[num_elements(fullknots)-2] knots;
  vector[2] bknots;
  fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  knots = segment(fullknots, 2, num_elements(fullknots)-2);
  bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
   int Nintk   = num_elements(knots);
   int degree  = 3;
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
    out = GS_ns_stan(x, knots, bknots, fullknots, allknots, 
                        N, degree, ord, 
                        Nintk, Nk, Nintervals, 
                        intercept, calcderiv, normalize, preH);                        
    if (centerval != 0) {
    matrix[1, df] cenout;
    cenout = GS_ns_stan(centerval+rep_vector(0.0, 1), knots, bknots, fullknots, allknots,
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
  
  return(support_GS_gps_parms_stan_str)
}

