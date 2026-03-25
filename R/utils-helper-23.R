
#############################################################
############### GS_gps_parms_R ##################
#############################################################

GS_gps_parms_R <- function(nlp_a,
                           nlp_b,
                           nlp_c,
                           nlp_d,
                           SParMat,
                           xknots,
                           spline_eval_array,
                           xg_array,
                           xg_curve_array,
                           degree,
                           shift_indicator,
                           spline_subset_indicator,
                           spline_precomputed_indicator,
                           set_spread,
                           return_indicator,
                           drawni) {
  
  pieces_dim      <- (length(xknots)-1)
  degree_dim      <- (degree + 1L)
  deriv_root      <- 2
  solvex          <- 0
  PiecePolyCoef   <- matrix(NA, degree_dim, pieces_dim)
  rrootsmat       <- matrix(NA, 1, pieces_dim)
  
  piece_id        <- NULL
  mark_peak       <- TRUE
  
  parm_deriv      <- c(0,1)
  curve_deriv     <- c(0,1,2)
  
  n_rowni         <- length(nlp_a)
  
  degree_dim_set_spread <- degree_dim * set_spread
  
  if(0 %in%  return_indicator) {
    coef_array  <- array(NA, dim = c(degree_dim, pieces_dim, n_rowni))
  }
  
  if(1 %in%  return_indicator) {
    parm_mat_dim <- 6
    if(!is.null(drawni)) {
      parm_mat_dim <- parm_mat_dim + 1
    }
    parm_array  <- array(NA, dim = c(pieces_dim, parm_mat_dim, n_rowni))
  }
  
  if(2 %in%  return_indicator) {
    deriv_array_dim_1 <- pieces_dim * degree_dim * set_spread
    deriv_array_dim_2 <- (length(curve_deriv) + 1+1+1)
    if(!is.null(drawni)) {
      deriv_array_dim_2 <- deriv_array_dim_2 + 1
    }
    deriv_array <- 
      array(NA, dim = c(deriv_array_dim_1, deriv_array_dim_2, n_rowni))
  }
  
  
  ############################################################
  # Start - over 1:n_rowni
  ############################################################
  for (rowni in 1:n_rowni) {
    if(spline_subset_indicator) {
      smat <- SParMat[1, ]
    } else {
      smat <- SParMat[rowni, ]
    }
    par_a <- nlp_a[rowni]
    par_b <- nlp_b[rowni]
    par_c <- nlp_c[rowni]
    par_d <- nlp_d[rowni]
    
    ############################################################
    # Start - coefficients and roots
    ############################################################
    i <- 1L
    while (i <= pieces_dim) {
      if(spline_precomputed_indicator == 0) {
        # xg <- seq.int(xknots[i], xknots[i + 1L], length.out = degree_dim)
        xg = seq_fun_R(xknots[i], xknots[i+1], degree_dim);
      }
      if(spline_precomputed_indicator == 1) {
        xg = xg_array[,i];
      }
      
      if(spline_precomputed_indicator == 0) {
        # eval(SplineCall) %*% smat
        SplineCall[[2]] <- quote(xg)
        eval_Spline = eval(SplineCall)
      }
      if(spline_precomputed_indicator == 1) {
        eval_Spline = spline_eval_array[,,i];
      }
      
      yg <- eval_Spline %*% smat
      # Imp, adjust for b and  c after eval(SplineCall) %*% smat
      xg <- ((xg/exp(par_c)) + par_b) - shift_indicator * xknots[i]
      
      
      Xg <- outer(xg, 0:degree, "^")
      A  <- base::crossprod(Xg)
      b  <- base::crossprod(Xg, yg)
      U  <- chol.default(A)
      PiecePolyCoef[,i] <- pc <- 
        base::backsolve(U, base::forwardsolve(t.default(U), b))
      if (deriv_root > 0) {
        pc.i <- pc[-seq_len(deriv_root)] * 
          choose(deriv_root:degree, deriv_root) * 
          factorial(deriv_root)
      }
      pc.i[1] <- pc.i[1] - solvex
      
      # croots <- base::polyroot(pc.i)
      # rroots <- Re(croots)[round(Im(croots), 10) == 0]
      
      rroots <- - pc.i[1] / pc.i[2]
      
      if (shift_indicator) {
        rroots <- rroots + xknots[i]
      }
      xi_  <- (xknots[i]   / exp(par_c)) + par_b
      xi_1 <- (xknots[i+1] / exp(par_c)) + par_b
      
      
      
      if((rroots >= xi_ ) & (rroots <= xi_1 )) {
        get_rroots <- rroots # [(rroots >= xi_) & (rroots <= xi_1 )]
      } else {
        get_rroots <- NA
      }
      
      # Also replace first and last piece to NA, boundry
      if(i == 1 | i == pieces_dim) {
        get_rroots <- NA
      }
    
      rrootsmat[,i] <- get_rroots
      i <- i + 1L
    } # while (i <= pieces_dim) {
    
    if(0 %in%  return_indicator) {
      coef_array[,,rowni]   <- PiecePolyCoef
    }
    
    ############################################################
    # End - coefficients and roots
    ############################################################
    
    
    ############################################################
    # Start - parameters
    ############################################################
    newx       <- rrootsmat %>% as.vector()
    if(is.null(piece_id)) {
      set_piece_id <- seq(1, pieces_dim)
      set_piece_id <- set_piece_id[!is.na(newx)]
      newx     <- newx[!is.na(newx)]
    } else {
      set_piece_id <- piece_id
      if(length(set_piece_id) != length(newx))
        stop("lengths of 'piece_id' and 'newx' must match")
    }
    ind <- split.default(seq_len(length(newx)), set_piece_id)
    unique_piece_id <- as.integer(names(ind))
    n_pieces <- length(unique_piece_id)
    
    parm_mat_dim  <- parm_deriv + 1+1+1+1
    if(mark_peak) {
      parm_mat_dim <- parm_mat_dim + 1
      if(!1 %in% parm_deriv) stop("please set parm_deriv = 1")
    }
    if(!is.null(rowni))  parm_mat_dim <- parm_mat_dim + 1
    if(!is.null(drawni)) parm_mat_dim <- parm_mat_dim + 1
    parm_mat_dim_nrows <- pieces_dim
    parm_mat           <- matrix(NA, parm_mat_dim_nrows, parm_mat_dim)
    
    i <- 1L
    
  
    # Otherwise error object 'collectd_mat_i' not found when ann NA matrix
    if(n_pieces == 0) {
      n_piecesTF <- FALSE 
      mark_peak  <- FALSE
    } else {
      n_piecesTF <- TRUE
    }
    
    # No need, even return(NULL) sufficient when n_piecesTF = FALSE
    if(!n_piecesTF) {
      out <- list()
      parm_mat[1, ncol(parm_mat)-3] <- 1
      parm_mat[, ncol(parm_mat)-2]  <- 1:nrow(parm_mat)
      parm_mat[, ncol(parm_mat)-1]  <- n_pieces
      parm_mat[, ncol(parm_mat)]    <- drawni
      parm_array[ , ,] <- parm_mat
      if(0 %in%  return_indicator) {
        out[['coef']]    <- coef_array
      }
      if(1 %in% return_indicator) {
        out[['parm']]   <- parm_array
      }
      if(2 %in% return_indicator) {
        out[['deriv']]    <- deriv_array
      }
      if(length(out) > 1) {
        return(out)
      } else if(length(out) == 1) {
        return(out[[1]])
      }
    } # if(!n_piecesTF) {
    
    
   
    while (i <= n_pieces) {
      ii      <- unique_piece_id[i]
      xg_newx <- newx[ind[[i]]]
      xg      <- xg_newx - shift_indicator * xknots[ii]
      pc      <- PiecePolyCoef[, ii]
      collectd_mat_i <- 0
      for (dxi in parm_deriv) {
        if(!is.null(parm_deriv)) {
          if (dxi >= degree) {
            stop("'parm_deriv' should be less than 'degree' i.e., ", degree)
          }
        }
        collectd_mat_i <- collectd_mat_i + 1
        if(dxi == 0) {
          pc.i <- pc
        } else if (dxi > 0) {
          pc.i <- pc[-seq_len(dxi)] * choose(dxi:degree, dxi) * factorial(dxi)
        }
        # already above shifted -> newx[ind[[i]]] - shift_indicator * x[ii]
        collectd_matx <- c(outer(xg - 0, 0:(degree - dxi), "^") %*% pc.i)
        if(dxi == 0) {
          collectd_matx <- collectd_matx + par_a
        }
        parm_mat[ii, 1]                    <- xg_newx
        parm_mat[ii, collectd_mat_i+1]     <- collectd_matx
      }
      i <- i + 1L
    } # while (i <= n_pieces) {
    add_incre_1 <- 0
    
    if(mark_peak) {
      ya_sitar <- parm_mat[, length(parm_deriv)+1] %>% as.vector()
      ya_sitar[is.na(ya_sitar)] <- 0
      d2_test <- diff(ya_sitar)
      
      peak_id   <- which(d2_test < 0)
      trough_id <- which(d2_test > 0)
      parm_mat[peak_id,   collectd_mat_i+1+1] <- 1
      parm_mat[trough_id, collectd_mat_i+1+1] <- 0
      # Also replace first and last piece to NA, boundary
      parm_mat[1,              collectd_mat_i+1+1] <- NA
      parm_mat[nrow(parm_mat), collectd_mat_i+1+1] <- NA
      add_incre_1 <- 1
    }
    
    parm_mat[, collectd_mat_i+1+1+add_incre_1]   <- seq(1, parm_mat_dim_nrows)
    parm_mat[, collectd_mat_i+1+1+1+add_incre_1] <- rowni
    if(!is.null(drawni)) {
      parm_mat[, collectd_mat_i+1+1+1+1+add_incre_1] <- drawni
    }

    if(1 %in%  return_indicator) {
      parm_array[,,rowni] <- parm_mat
    }
    
    
    ############################################################
    # End - parameters
    ############################################################
    
    ############################################################
    # Start - derivatives
    ############################################################
    
    collectd_mat_dim                      <- length(curve_deriv) + 1+1+1
    if(!is.null(drawni)) collectd_mat_dim <- collectd_mat_dim + 1
    collectd_mat_nrows <- 0
    collectd_mat  <- matrix(NA, collectd_mat_nrows, collectd_mat_dim) 
    
    i <- 1L
    while (i <= pieces_dim) {
      if(spline_precomputed_indicator == 0) {
        xg_curve = seq_fun_R(xknots[i], xknots[i+1], degree_dim_set_spread);
      }
      if(spline_precomputed_indicator == 1) {
        xg_curve = xg_curve_array[,i];
      }
      xg_curve <- 
        seq.int(xknots[i], xknots[i + 1L], 
                length.out = set_spread * degree_dim)
      xg_curve <- ((xg_curve/exp(par_c)) + par_b)
      pc <- PiecePolyCoef[, i]
      collectd_mat_x      <- matrix(NA, length(xg_curve), collectd_mat_dim)
      collectd_mat_i <- 0
      for (dxi in curve_deriv) {
        if(!is.null(curve_deriv)) {
          if (dxi >= degree) {
            stop("'curve_deriv' should be less than 'degree' i.e., ", degree)
          }
        }
        collectd_mat_i <- collectd_mat_i + 1
        if(dxi == 0) {
          pc.i <- pc
        } else if (dxi > 0) {
          pc.i <- pc[-seq_len(dxi)] * choose(dxi:degree, dxi) * factorial(dxi)
        }
        collectd_matx <- 
          c(outer(xg_curve - 
                    shift_indicator * xknots[i], 
                  0:(degree - dxi), "^") %*% pc.i)
        
        if(dxi == 0) {
          collectd_matx <- collectd_matx + par_a
        }
        
        collectd_mat_x[, 1]                    <- xg_curve
        collectd_mat_x[, collectd_mat_i+1]     <- collectd_matx
        collectd_mat_x[, collectd_mat_i+1+1]   <- i
        collectd_mat_x[, collectd_mat_i+1+1+1] <- rowni
        if(!is.null(drawni)) collectd_mat_x[, collectd_mat_i+1+1+1+1] <- drawni
      }
      collectd_mat <- rbind(collectd_mat, collectd_mat_x)
      i <- i + 1L
    } # while (i <= pieces_dim) {
    
    if(2 %in%  return_indicator) {
      deriv_array[,,rowni] <- collectd_mat
    }
    
    ############################################################
    # End - derivatives
    ############################################################
    
  } # for (rowni in 1:n_rowni) {
  
  ############################################################
  # End - over 1:n_rowni
  ############################################################
  
  out <- list()
  if(0 %in%  return_indicator) {
    out[['coef']]    <- coef_array
  }
  if(1 %in% return_indicator) {
    out[['parm']]   <- parm_array
  }
  if(2 %in% return_indicator) {
    out[['deriv']]    <- deriv_array
  }
  if(length(out) > 1) {
    return(out)
  } else if(length(out) == 1) {
    return(out[[1]])
  }
  
} # End GS_gps_parms_R


#############################################################
############ -------- seq_fun_R -------- ##########
#############################################################


seq_fun_R <- function(start, end, N_by) { 
  h = (end - start) / (N_by - 1)
  out_c <- c()
  for (i in 1:N_by) { 
    out=start + (i - 1) * h
    out_c <- c(out_c, out)
  }
  return(out_c)
}


#############################################################
############ -------- my_counter -------- ##########
#############################################################

counter_function <- function() {
  env <- environment()
  env$counter <- 0
  increment <- function() {
    env$counter <- env$counter + 1
    return(env$counter)
  }
  get_environment <- function() {
    return(env)
  }
  reset_counter <- function(new_value = 0) {
    env$counter <- new_value
  }
  list(next_value = increment,
       get_environment = get_environment,
       reset = reset_counter)
}
my_counter <- counter_function()

#############################################################
# -------- wraper_for_drawni function for w/t future ------ #
#############################################################

wraper_for_drawni <- function(setdat_mat, 
                              drawni, 
                              callvia,
                              return_indicator, 
                              subset_data_by,
                              subset_data_by_names,
                              create_abcd_names_vector,
                              create_s_names_vector,
                              xknots,
                              degree,
                              spline_subset_indicator = 
                                spline_subset_indicator,
                              spline_precomputed_indicator = 
                                spline_precomputed_indicator,
                              shift_indicator,
                              set_spread,
                              spline_eval_array,
                              xg_array,
                              xg_curve_array,
                              call_function,
                              GS_gps_parms_assign) {
  if(callvia == 'base') {
    drawniid <- drawni
  } else if(callvia == 'future') {
    # When plan multisession, my_counter will generate seq per session
    # to get unique sequence later assign_new_sequence
    drawniid <- my_counter$next_value() 
    pid <- Sys.getpid()
    time_us <- as.numeric(Sys.time()) * 1e6
    # session_prefix <- paste0(pid, "_", floor(time_us))
    # drawniid <- paste0(session_prefix, "_", drawniid) #%>% as.numeric()
    drawniid <- drawniid+time_us+drawniid
  }
  
  if(!is.null(subset_data_by)) {
    setdat_mat <- unique(data.table::as.data.table(setdat_mat),
                                                 by= subset_data_by_names) %>% 
      as.matrix() 
  }
  set_frame_abcd <- setdat_mat[, create_abcd_names_vector, drop = FALSE]
  set_frame_smat <- setdat_mat[, create_s_names_vector,    drop = FALSE]
 

  mat_parm <-  GS_gps_parms_assign(
    nlp_a = set_frame_abcd[, 'a'],
    nlp_b = set_frame_abcd[, 'b'],
    nlp_c = set_frame_abcd[, 'c'],
    nlp_d = set_frame_abcd[, 'd'],
    SParMat = set_frame_smat,
    xknots = xknots,
    spline_eval_array = spline_eval_array,
    xg_array = xg_array,
    xg_curve_array = xg_curve_array,
    degree = degree,
    shift_indicator = shift_indicator,
    spline_subset_indicator = spline_subset_indicator,
    spline_precomputed_indicator = spline_precomputed_indicator,
    set_spread = set_spread,
    return_indicator = return_indicator,
    drawni = drawniid)
  return(mat_parm)
}



#############################################################
# ------ wraper_for_drawni_2 function for w/t future ----- #
#############################################################

wraper_for_drawni_2 <- function(drawni, 
                                nlpar_fixed, 
                                nlpar_random,
                                create_s_names_vector,
                                add_xtm,
                                xoffset,
                                d_adjusted,
                                SplineCall_d0,
                                SplineCall_d1,
                                mat.adj,
                                funx_,
                                ifunx_,
                                callvia) {
  
  if(callvia == 'base') {
    drawniid <- drawni
  } else if(callvia == 'future') {
    drawniid <- my_counter$next_value()
    pid <- Sys.getpid()
    time_us <- as.numeric(Sys.time()) * 1e6
    drawniid <- drawniid+time_us+drawniid
  }
  
  setdat_mat_fixed  <- nlpar_fixed [drawni, ,]
  setdat_mat_random <- nlpar_random[drawni, ,]
  spmat             <- setdat_mat_fixed[, create_s_names_vector]
  setx0             <- setdat_mat_fixed[, 'Xestimate']
  
  # if(all(setdat_mat_fixedx[, 'd'] == 0)) {
  #   d_parmTF <- FALSE
  # } else {
  #   d_parmTF <- TRUE
  # }
  
  # atgv might have Na
  if(all(is.na(setx0))) {
    x.adj <- NA_real_
    y.adj <- NA_real_
    v.adj <- NA_real_
  } else {
    setx              <- funx_(setx0)
    setx              <- setx - setdat_mat_fixed[, 'b']
    x.adj             <- setx/exp(setdat_mat_random[,"c"]) + 
      setdat_mat_random[,"b"] + 
      setdat_mat_fixed[, 'b']
    x.adj             <- ifunx_(x.adj)
    setx0             <- setx * exp(setdat_mat_fixed[, 'c'])
    y.adj             <- rowSums(eval(SplineCall_d0) * spmat) + 
      setdat_mat_random[,"a"] + 
      setdat_mat_random[,"d"] * x.adj +
      setdat_mat_fixed[, 'a']
    v.adj             <- rowSums(eval(SplineCall_d1) * spmat) * 
      exp(setdat_mat_fixed[, 'c'] + 
            setdat_mat_random[,"c"]) + 
      setdat_mat_random[,"d"] 
  } # end else if(all(is.na(setx0))) {
  
  if(add_xtm) {
    setxx.adj_xtm <- (setdat_mat_fixed[, 'xvar'])
    # setxx.adj_xtm <- funx_(setdat_mat_fixed[, 'xvar'])
    x.adj_xtm <- (setxx.adj_xtm - setdat_mat_random[,"b"]) * 
      exp(setdat_mat_random[,"c"])
    
    if( d_adjusted) dmultiplier <- x.adj_xtm
    if(!d_adjusted) dmultiplier <- setxx.adj_xtm
    
    y.adj_xtm <- setdat_mat_fixed[, 'yvar'] - 
      setdat_mat_random[,"a"] - 
      setdat_mat_random[,"d"] * dmultiplier # setxx.adj_xtm
    
    x.adj_xtm <- ifunx_(x.adj_xtm)
  }
  
  mat.adj[, 1] <- x.adj
  mat.adj[, 2] <- y.adj
  mat.adj[, 3] <- v.adj
  if(add_xtm) {
    mat.adj[, 4] <- x.adj_xtm
    mat.adj[, 5] <- y.adj_xtm
    mat.adj[, 6] <- drawniid
    mat.adj[, 7] <- setdat_mat_fixed[, 'fomerge']
    mat.adj[, 8] <- setdat_mat_fixed[, 'xid']
  } else {
    mat.adj[, 4] <- drawniid
    mat.adj[, 5] <- setdat_mat_fixed[, 'fomerge']
    mat.adj[, 6] <- setdat_mat_fixed[, 'xid']
  }
  return(mat.adj)
}





#############################################################
############ -------- assign_new_sequence -------- ##########
#############################################################

assign_new_sequence <- function(mat, col) {
  matrix_column <- as.matrix(mat[,col])
  new_sequence <- integer(nrow(matrix_column))
  current_value <- NULL
  sequence_number <- 0
  for (i in 1:nrow(matrix_column)) {
    if (is.null(current_value) || matrix_column[i, 1] != current_value) {
      sequence_number <- sequence_number + 1
      current_value <- matrix_column[i, 1]
    }
    new_sequence[i] <- sequence_number
  } 
  mat[,col] <- new_sequence
  return(mat)
}


#############################################################
############### brms_posterior_summary #####################
#############################################################

brms_posterior_summary <- function(.x, 
                                   probs = c(0.025, 0.975), 
                                   robust = FALSE, 
                                   ...) {
  brms::posterior_summary(.x, probs = probs, robust = robust)
}


#############################################################
############### collapse_posterior_summary ##################
#############################################################

collapse_posterior_summary <- function(.x, 
                                       probs = c(0.025, 0.975), 
                                       robust = FALSE, 
                                       setcolnames = FALSE, 
                                       ...) {
  collapse_mad <- function(x, constant = 1.4826) {
    collapse::fmedian(abs(x - collapse::fmedian(x))) * constant
  }
  if(!robust) {
    out <- c(collapse::fmean(.x), 
             collapse::fsd(.x), 
             collapse::fquantile(.x, probs = probs) 
    ) 
  } else if(robust) {
    out <-  c(collapse::fmedian(.x), 
              collapse_mad(.x), 
              collapse::fquantile(.x, probs = probs) 
    ) 
  }
  
  if(setcolnames) {
    colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 
                                                         100))
  }
  
  return(out)
} # end collapse_posterior_summary







#############################################################
# -parameter_method_loop_over_parm function for w/t future- #
#############################################################

parameter_method_loop_over_parm <- function(parm,
                                            set_loop_over_parm_last,
                                            modelbased_arguments,
                                            by,
                                            xoffset,
                                            d_adjusted,
                                            subset_data_by,
                                            set_pdrawsp,
                                            set_pdraws,
                                            newdata,
                                            draw_ids_seq,
                                            set_draws_n,
                                            method_call,
                                            set_dataf_m_collapse,
                                            set_nrows_n,
                                            xvar,
                                            yvar,
                                            future,
                                            add_xtm,
                                            nlpar_fixed,
                                            nlpar_random,
                                            create_s_names_vector,
                                            SplineCall_d0,
                                            SplineCall_d1,
                                            # mat.adj,
                                            funx_,
                                            ifunx_,
                                            callvia,
                                            get_data_cols.org,
                                            ec_agg,
                                            ei_agg,
                                            nthreads,
                                            conf,
                                            probs,
                                            na.rm = TRUE,
                                            verbose = FALSE) {
  . <- NULL;
  d0 <- NULL;
  d1 <- NULL;
  x <- NULL;
  xtm <- NULL;
  ytm <- NULL;
  xid <- NULL;
  
  if(add_xtm) {
    if(parm != set_loop_over_parm_last) add_xtm <- FALSE
    if(parm == set_loop_over_parm_last) add_xtm <- add_xtm
  }

  get_growthparameters_args <- modelbased_arguments
  # This NA because we want to plugin the population average apgv
  get_growthparameters_args[['re_formula']] <- NA
  get_growthparameters_args[['pdrawsp']]    <- set_pdrawsp
  get_growthparameters_args[['pdraws']]     <- set_pdraws
  get_growthparameters_args[['parameter']]  <- parm
  get_growthparameters_args[['newdata']]    <- newdata
  get_growthparameters_args[['draw_ids']]   <- draw_ids_seq
  get_growthparameters_args[['by']]         <- by
  
  get_growthparameters_args[['newdata_fixed']] <- 0
  get_growthparameters_args[['reformat']] <- FALSE
  
  get_growthparameters_args[['method']] <- method_call
 
  onex0 <- CustomDoCall(get_growthparameters, get_growthparameters_args)

  onex00 <- onex0
  onex00[['draw']] <- funx_(onex00[['draw']])
  onex00 <- set_dataf_m_collapse %>% collapse::join(onex00, on = by,
                                                    how = "left",
                                                    multiple = TRUE,
                                                    verbose = FALSE)
  
  onex00 <- data.table::as.data.table(onex00)
  
  array_dim     <- set_nrows_n
  pieces_dim    <- 1
  parm_mat_dim  <- 6
  
  onex00 <- onex00 %>% 
    collapse::fmutate(fomerge =  
                        collapse::finteraction(onex00 %>% 
                                                 collapse::fselect(c(by)), 
                                               factor = FALSE))
  
  xid_by_onex00 <- c("drawid", "parameter", 'fomerge')
  
  onex00$xid <- setorderv(onex00, xid_by_onex00)[, .(xid=seq_len(.N)), 
                                                 by = xid_by_onex00]$xid
  
  which_dim                  <- 3
  nlpar_fixed_names_dim3     <- attr(nlpar_fixed, "dimnames")[[which_dim]]
  nlpar_fixed_names_dim3_add <- c('Xestimate', 'fomerge', 'xid')
  
  
  extend_array <- cbind(onex00[['draw']], 
                        onex00[['fomerge']], 
                        onex00[['xid']])
  
  if(add_xtm) {
    extend_array <- cbind(extend_array, newdata[[xvar]], newdata[[yvar]])
    nlpar_fixed_names_dim3_add <- c(nlpar_fixed_names_dim3_add,
                                    "xvar", 'yvar')
    parm_mat_dim <- parm_mat_dim + 2
  }
  
  Sliced <- aperm(`dim<-`(t(extend_array), 
                          c(ncol(extend_array), 
                            dim(nlpar_fixed)[2], 
                            dim(nlpar_fixed)[1])), c(3, 2, 1))
  
  nlpar_fixed <- abind::abind(nlpar_fixed, Sliced, along = 3)
  
  attr(nlpar_fixed, "dimnames")[[which_dim]] <- c(nlpar_fixed_names_dim3, 
                                                  nlpar_fixed_names_dim3_add)
  
  mat.adj            <- matrix(NA_real_, nrow = array_dim, ncol = parm_mat_dim)
  
  mat.adj            <- data.table::as.data.table(mat.adj)
  if(!future) {
    collect_draws_parm <- list()
    for (drawni in 1:set_draws_n) {
      collect_draws_parm[[drawni]] <- 
        wraper_for_drawni_2(drawni = drawni, 
                            nlpar_fixed = nlpar_fixed,
                            nlpar_random = nlpar_random,
                            create_s_names_vector = create_s_names_vector,
                            add_xtm = add_xtm,
                            xoffset = xoffset,
                            d_adjusted = d_adjusted,
                            SplineCall_d0 = SplineCall_d0,
                            SplineCall_d1 = SplineCall_d1,
                            mat.adj = mat.adj,
                            funx_ = funx_,
                            ifunx_ = ifunx_,
                            callvia = callvia)
    } # for (drawni in 1:set_draws_n) {
  } # if(!future) {
  
  if(future) {
    # setup future
    environment(wraper_for_drawni_2) <- environment()
    future_globals_list = list( mat.adj = mat.adj,
                                `%>%` = bsitar::`%>%`,
                                my_counter = my_counter)
    # call future
    my_counter$reset()
    collect_draws_parm <- future.apply::future_lapply(
      1:set_draws_n, 
      FUN = function(drawni, ...) 
        wraper_for_drawni_2(drawni = drawni, 
                            nlpar_fixed = nlpar_fixed,
                            nlpar_random = nlpar_random,
                            create_s_names_vector = create_s_names_vector,
                            add_xtm = add_xtm,
                            xoffset = xoffset,
                            d_adjusted = d_adjusted,
                            SplineCall_d0 = SplineCall_d0,
                            SplineCall_d1 = SplineCall_d1,
                            mat.adj = mat.adj,
                            funx_ = funx_,
                            ifunx_ = ifunx_,
                            callvia = 'future'),
      future.globals = future_globals_list)
  } # end else if(future) {
  

  
  if(add_xtm) {
    names_parm      <- c("x", "d0", "d1", "xtm", "ytm", "drawid")
  } else {
    names_parm      <- c("x", "d0", "d1", "drawid")
  }
  names_parm_temp <- c(names_parm, "fomerge", "xid")
  
  bind_draws_parm <- collect_draws_parm %>% collapse::rowbind() %>%
    collapse::setrename(names_parm_temp)
  
  if(future) {
    which_cols <- 4
    if(add_xtm) {
      which_cols <- which_cols + 2
    }
    bind_draws_parm <- assign_new_sequence(mat = bind_draws_parm %>% 
                                             as.matrix(), 
                                           col = which_cols) %>% 
      data.table::as.data.table()
  }
  
  peak_data_draw <- bind_draws_parm %>% 
    collapse::join(onex00, on = c("drawid", 'fomerge', "xid"),
                   how = "right",
                   multiple = FALSE,
                   verbose = FALSE) 
  
  
  peak_data_draw_select <- c(names_parm, get_data_cols.org)
  peak_data_draw        <- collapse::fselect(peak_data_draw, 
                                             peak_data_draw_select)
  
  # One can subset peak_data_draw but not xtm_data_draw
  xtm_data_draw <- peak_data_draw
  
  # new
  if(!is.null(subset_data_by)) {
    if(!subset_data_by) subset_data_by <- NULL
  }
  
  
  if(!is.null(subset_data_by)) {
    group_by_indices <- c("drawid", subset_data_by) 
    peak_data_draw <- peak_data_draw[peak_data_draw[, .I[1:1], 
                                                    by = group_by_indices]$V1]
  }
  
  xid_by <- c("drawid", "parameter", "id") 
  

  if(parm == 'apgv') {
    parameter_names_vec <- c('apgv', 'pgv', 'spgv')
  }
  if(parm == 'atgv') {
    parameter_names_vec <- c('atgv', 'tgv', 'stgv')
  }
  if(parm == 'acgv') {
    parameter_names_vec <- c('acgv', 'cgv', 'scgv')
  }
  
  if(nrow(peak_data_draw) > 0) {
    apgv_draw    <- peak_data_draw %>% 
      collapse::fmutate(draw = x) %>% 
      collapse::fmutate(parameter = parameter_names_vec[1]) 
    pgv_draw    <- peak_data_draw %>% 
      collapse::fmutate(draw = d1) %>% 
      collapse::fmutate(parameter = parameter_names_vec[2]) 
    spgv_draw    <- peak_data_draw %>% 
      collapse::fmutate(draw = d0) %>% 
      collapse::fmutate(parameter = parameter_names_vec[3])
    
    all_peak_data_draw <- collapse::rowbind(apgv_draw, pgv_draw, spgv_draw)
    
    all_peak_data_draw$xid <- setorderv(all_peak_data_draw, 
                                        xid_by)[, .(xid=seq_len(.N)), 
                                                by = xid_by]$xid
    
    
    # apply ifunx_
    all_peak_data_draw[['draw']] <- ifunx_(all_peak_data_draw[['draw']])
    
    get_growthparameters_args <- modelbased_arguments
    get_growthparameters_args[['preparms']] <- all_peak_data_draw
    get_growthparameters_args[['by']] <- c(by, 'xid')
    # For preparms, method must be custom
    get_growthparameters_args[['method']] <- method_call
    if(add_xtm) {
      get_growthparameters_args[['pdrawsp']] <- F
    }
    get_growthparameters_args[['reformat']] <- FALSE
    
    peak_parameters <- CustomDoCall(get_growthparameters,
                                    get_growthparameters_args)
    
    
    peak_names.ors__ <- colnames(peak_parameters)
    data.table::setnames(peak_parameters, tolower(names(peak_parameters)))
    
    peak_roworderv_vars <- 'parameter'
    peak_roworderv_vars <- c(peak_roworderv_vars, by, 'xid')
    peak_parameters <- collapse::roworderv(peak_parameters, peak_roworderv_vars)
    
    if(add_xtm) {
      xtm_draw    <- xtm_data_draw %>% 
        collapse::fmutate(draw = xtm) %>% 
        collapse::fmutate(parameter = 'xtm') 
      ytm_draw    <- xtm_data_draw %>% 
        collapse::fmutate(draw = ytm) %>% 
        collapse::fmutate(parameter = 'ytm') 
      all_tm_data_draw <- collapse::rowbind(xtm_draw, ytm_draw)
      
      all_tm_data_draw$xid <- setorderv(all_tm_data_draw,
                                        xid_by)[, .(xid=seq_len(.N)),
                                                by = xid_by]$xid
      
      get_growthparameters_args <- modelbased_arguments
      get_growthparameters_args[['preparms']] <- all_tm_data_draw
      get_growthparameters_args[['by']] <- c(by, 'xid')
      
      get_growthparameters_args[['method']] <- method_call
      if(add_xtm) {
        get_growthparameters_args[['pdrawsp']] <- TRUE
      }
      get_growthparameters_args[['reformat']] <- FALSE
      
      tm_parameters <- CustomDoCall(get_growthparameters,
                                    get_growthparameters_args)
      
      tm_names.ors__ <- colnames(tm_parameters)
      data.table::setnames(tm_parameters, tolower(names(tm_parameters)))
      tm_roworderv_vars <- 'parameter'
      tm_roworderv_vars <- c(tm_roworderv_vars, by, 'xid')
      tm_parameters <- collapse::roworderv(tm_parameters, tm_roworderv_vars)
      tm_parameters <- data.table::setnames(tm_parameters, tm_names.ors__)
      # peak_parameters <- collapse::rowbind(peak_parameters, tm_parameters)
    } # if(add_xtm) {
  } # if(nrow(peak_data_draw) > 0) {
  
  if(add_xtm) {
    # tm_parameters[['xtm']] <- ifunx_(tm_parameters[['xtm']])
    setdrawidparm <- c(by, 'xid')
    namesx <- c('estimate', 'conf.low', 'conf.high')
    namesx <- paste0("x", ".", namesx)
    setdrawidparm_ <- c(setdrawidparm, namesx)
    tm_parameters_xtm <- 
      tm_parameters %>% collapse::fgroup_by(setdrawidparm) %>% 
      collapse::fsummarise(collapse::mctl(
        get_pe_ci_collapse(.data[['xtm']],
                           ec_agg = ec_agg, 
                           ei_agg = ei_agg, na.rm = TRUE, 
                           nthreads = nthreads, 
                           conf = conf, probs = probs))
      ) %>% collapse::frename(., setdrawidparm_) 
    
    namesx <- c('estimate', 'conf.low', 'conf.high')
    namesx <- paste0("y", ".", namesx)
    setdrawidparm_ <- c(setdrawidparm, namesx)
    tm_parameters_ytm <- 
      tm_parameters %>% collapse::fgroup_by(setdrawidparm) %>% 
      collapse::fsummarise(collapse::mctl(
        get_pe_ci_collapse(.data[['ytm']],
                           ec_agg = ec_agg, 
                           ei_agg = ei_agg, na.rm = TRUE, 
                           nthreads = nthreads, 
                           conf = conf, probs = probs))
      ) %>% collapse::frename(., setdrawidparm_) 
    
    tm_parameters_xtm_ytm <- 
      tm_parameters_xtm %>% collapse::join(tm_parameters_ytm, 
                                          on = setdrawidparm,
                                          how = "left",
                                          multiple = TRUE,
                                          verbose = FALSE)
    
    tm_parameters_xtm_ytm   <- tm_parameters_xtm_ytm %>% collapse::fselect(-xid)
  } # if(add_xtm) {
  
  if(nrow(peak_data_draw) == 0) {
    peak_parameters <- NULL
  }
  
  if(is.null(peak_parameters)) return(peak_parameters)
  
  peak_parameters   <- peak_parameters %>% collapse::fselect(-xid)
  
  peak_names.ors__2 <- peak_names.ors__[ !grepl('xid', peak_names.ors__)]
  
  indices_lastn   <- 2 # for est and q1 qu
  firstup_indices <- 
    (length(peak_names.ors__2)-indices_lastn):length(peak_names.ors__2)
  lower_case__2 <- 
    peak_names.ors__2[1:(length(peak_names.ors__2)-indices_lastn-1)]
  upper_case__2 <- firstup(peak_names.ors__2[firstup_indices])
  
  peak_names.ors__2 <- c(lower_case__2, upper_case__2)
  
  peak_parameters <- data.table::setnames(peak_parameters, peak_names.ors__2)

  peak_parameters <- DT_to_data_frames(peak_parameters)
  
  if(add_xtm) {
    attr(peak_parameters, 'xtm_ytm') <-     tm_parameters_xtm_ytm
  }
  
  return(peak_parameters)
} # parameter_method_loop_over_parm



