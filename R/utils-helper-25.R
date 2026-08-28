


fast_nsk_rcsfun_str_get <- function() {
  fast_nsk_rcsfun_str <- 
    "
  vector spline_geths(vector nodes) {
    int n = size(nodes) - 1;
    vector[n] hs;
    for (i in 1:n) {
      hs[i] = nodes[i + 1] - nodes[i];
    }
    return hs;
  }

  vector spline_getcoeffs(vector nodes, vector vals) {
    int n_nodes = size(nodes);
    int n = n_nodes - 1;
    vector[n] hi;
    vector[n] bi;
    vector[n - 1] vi;
    vector[n - 1] ui;
    vector[n_nodes] ret;
    vector[n - 1] zs;
    matrix[n - 1, n - 1] M = rep_matrix(0, n - 1, n - 1);
    n = n_nodes - 1;
    for (i in 1:n) {
      hi[i] = nodes[i + 1] - nodes[i];
      bi[i] = 1 / hi[i] * (vals[i + 1] - vals[i]);
    }
    for (i in 2:n) {
      vi[i - 1] = 2 * (hi[i - 1] + hi[i]);
      ui[i - 1] = 6 * (bi[i] - bi[i - 1]);
    }
    for (i in 1:n - 1) {
      M[i, i] = vi[i];
    }
    for (i in 1:n - 2) {
      M[i + 1, i] = hi[i + 1];
      M[i, i + 1] = hi[i + 1];
    }
    zs = M \\ ui;
    ret[1] = 0;
    ret[n_nodes] = 0;
    ret[2:n_nodes - 1] = zs;
    return ret;
  }

  vector spline_eval(vector nodes,
                    vector vals, vector zs,
                    vector x, array[] int i) {
    int n_nodes = size(nodes);
    int n_dat = size(x);
    vector[n_nodes - 1] h;
    vector[n_dat] ret;
    array[n_dat] int i1;
    for (ii in 1:n_dat) {
      i1[ii] = i[ii] + 1;
    }
    h = spline_geths(nodes);
    ret = (
          zs[i1] ./ 6 ./ h[i] .* square(x - nodes[i]) .* (x - nodes[i]) +
          zs[i] ./ 6 ./ h[i] .* square(nodes[i1] - x) .* (nodes[i1] - x) +
          (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6) .* (x - nodes[i]) +
          (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6) .* (nodes[i1] - x)
    );
    return ret;
  }

  array[] int spline_findpos(vector nodes, vector x, vector b,  vector c) {
    int n_nodes = size(nodes);
    int n_dat = size(x);
    array[n_dat] int ret;
    for (i in 1:n_dat) {
      real bx =  b[i];
      real cx =  c[i];
      real setx = x[i];
      int success = 0;
      for (j in 1:n_nodes - 1) {
      real nodesj  = (nodes[j] - bx) * cx ;
      real nodesj1 = (nodes[j+1] - bx) * cx ;
        if ((setx >= nodesj) && (setx < nodesj1 )) {
          ret[i] = j;
          success = 1;
          break;
        }
      }
      if (success == 0) {
        reject(\"Point outside knot\");
      }
    }
    return ret;
  }
  "
  return(fast_nsk_rcsfun_str)
} 






