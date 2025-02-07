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
