cutpoints_f = function(g, ordinal_counts, ordinal_latent_sum, Egamma_j, n_levels){
  value = 0
  # can probable vectorise this
  for(l in 1:n_levels){
    value = value + ordinal_counts[l] * log(g[l+1]-g[l]) - 
      ordinal_counts[l] * Egamma_j/6*(g[l+1]**2 + g[l]**2 + g[l+1]*g[l]) +
      0.5*Egamma_j*(g[l+1]+g[l])*ordinal_latent_sum[l]
  }
  return(-value)
}

cutpoints_grad = function(g, ordinal_counts, ordinal_latent_sum, Egamma_j, n_levels){
  grad = rep(0, n_levels+1)
  for(l in 2:n_levels){
    grad[l] = ordinal_counts[l-1] / (g[l]-g[l-1]) -
      ordinal_counts[l-1] * Egamma_j/6*(2*g[l]+g[l-1]) +
      0.5*Egamma_j*ordinal_latent_sum[l-1] + 
      ordinal_counts[l] / (g[l+1]-g[l]) -
      ordinal_counts[l] * Egamma_j/6*(2*g[l]+g[l+1]) +
      0.5*Egamma_j*ordinal_latent_sum[l]
  }
  return(-grad)
}

generate_constraints = function(n_levels, G=5){
  A = matrix(0, n_levels, n_levels+1)
  for(k in 1:n_levels){
    A[k, k:(k+1)] = c(-1, 1)
  }
  b = rep(0, n_levels)
  A2 = rbind(c(1, rep(0, n_levels)), c(rep(0, n_levels), 1))
  b2 = c(-G - 1e-8, G - 1e-8)
  A3 = -A2
  b3 = c(G - 1e-8, -G - 1e-8)
  ui = rbind(A, A2, A3)
  ci = c(b, b2, b3)
  return(list(ui = ui, ci = ci))
}
