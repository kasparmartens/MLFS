update_V_sigmainv = function(Egamma, Eww, R, M){
  sigmainv = diag(R)
  for(j in 1:M){
    sigmainv = sigmainv + Egamma[j] * Eww[[j]]
  }
  return(sigmainv)
}

update_V_mu_individual_i = function(i, Eu, Ew, Egamma, M){
  mu = 0
  for(j in 1:M){
    # same update rule for all types
    mu = mu + Egamma[j] * Eu[[j]][i, ] %*% t(Ew[[j]])
  }
  return(mu)
}

update_sigma_W = function(Ealpha, Egamma, Evv_sum, M){
  sigma = list()
  for(j in 1:M){
    sigmainv = diag(Ealpha[j, ]) + Egamma[j] * Evv_sum
    sigma[[j]] = solve(sigmainv)
  }
  return(sigma)
}

update_mu_W = function(Eu, Ev, Ew, Egamma, type, d, M, sigma_W){
  mu = list()
  for(j in 1:M){
    mu[[j]] = sigma_W[[j]] * Egamma[j]
    temp = sapply(1:ncol(Ew[[j]]), function(which_col){
      t(Ev) %*% Eu[[j]][, which_col]
    })
    mu[[j]] = mu[[j]] %*% temp
  }
  return(mu)
}

update_alpha = function(Ealpha, Eww, d, M, aAlpha, bAlpha){
  for(j in 1:M){
    aTmp = aAlpha + 0.5*d[j]
    bTmp = bAlpha + 0.5*diag(Eww[[j]])
    Ealpha[j, ] = aTmp / bTmp
  }
  return(Ealpha)
}

update_gamma_j_gaussian = function(j, Eu, Ew, Eww, Ev, Evv_sum){
  b = 0.5*trace(Eww[[j]] %*% Evv_sum) - trace(t(Ew[[j]]) %*% t(Ev) %*% Eu[[j]]) + 0.5*sum(Eu[[j]]**2)
  return(b)
}

update_gamma_j_ordinal = function(j, g, n_levels, ordinal_counts, Eu, Ev, Ew, sigma_W, Evv_sum, d){
  g = g[[j]]
  latent = Ev %*% Ew[[j]]
  b = - sum(2*Eu[[j]] * latent) # note the elementwise multiplication
  for(l in 1:n_levels[[j]]){
    b = b + 1/3*(g[l+1]**2 + g[l]**2 + g[l+1]*g[l])*ordinal_counts[[j]][l]
  }
  temp_W = Ew[[j]] %*% t(Ew[[j]]) + d[j]*sigma_W[[j]]
  b = b + trace(temp_W %*% Evv_sum)
  return(0.5*b)
}

update_gamma_j_similarity = function(j, Eu, Euu, Euu_sum, Ew, Eww, Ev, Evv_sum, X, N, scaling_const){
  b = trace(Eww[[j]] %*% Evv_sum) - 2*trace(t(Ew[[j]]) %*% t(Ev) %*% Eu[[j]]) + trace(Euu_sum[[j]])
  temp = sum(X[upper.tri(X)]**2)
  for(i in 1:(N-1)){
    sum_trace = 0
#     for(k in (i+1):N){
#       sum_trace = sum_trace + sum(Euu[[j]][[i]] * Euu[[j]][[k]]) # trace(Euu[[j]][[i]] %*% Euu[[j]][[k]])
#     }
    sum_trace = sum(do.call("rbind", Euu[[j]][rep(i, N-i)]) * do.call("rbind", Euu[[j]][(i+1):N]))
    temp = temp - 2*as.numeric(X[i, (i+1):N, drop=FALSE] %*% Eu[[j]][(i+1):N, , drop=FALSE] %*% Eu[[j]][i, ]) + sum_trace
  }
  return(0.5*(b + scaling_const*temp))
}

update_gamma = function(j, g, n_levels, ordinal_counts, Eu, Euu, Euu_sum, Ev, Ew, Eww, sigma_W, Evv_sum, X, M, N, d, type, aGamma, bGamma, scaling_const){
  a_tilde = rep(NA, M)
  b_tilde = rep(NA, M)
  for(j in 1:M){
    a_tilde[j] = aGamma + d[j] * N/2
    if(type[j] == "gaussian"){
      b = update_gamma_j_gaussian(j, Eu, Ew, Eww, Ev, Evv_sum)
    }
    else if(type[j] == "ordinal"){
      b = update_gamma_j_ordinal(j, g, n_levels, ordinal_counts, Eu, Ev, Ew, sigma_W, Evv_sum, d)
    }
    else if(type[j] == "similarity"){
      a_tilde[j] = a_tilde[j] + 0.25*N*(N-1)
      b = update_gamma_j_similarity(j, Eu, Euu, Euu_sum, Ew, Eww, Ev, Evv_sum, X, N, scaling_const)
    }
    b_tilde[j] = bGamma + b
  }
  return(list(a_tilde = a_tilde, b_tilde = b_tilde))
}

update_Z = function(y, N, C, Ez, sigma_beta, Ebetabeta, rho, n_samples = 1000){
  muZ = Ez
  lowerbound_yz = -0.5*rho*trace(Ebetabeta)+0.5*C*logdet(sigma_beta)
  for(i in 1:N){
    others = setdiff(1:C, y[i])
    z_diff = Ez[i, y[i]] - Ez[i, others]
    # update for Ez[i, -y[i]], by replacing expectation with empirical approximation
    z_rep = do.call("rbind", lapply(1:n_samples, function(x)z_diff))
    u = rnorm(n_samples)
    u_rep = do.call("cbind", lapply(1:(C-1), function(x)u))
    cdfs = pnorm(u_rep + z_rep) + 1e-16
    pdfs = dnorm(u_rep + z_rep) + 1e-16
    products = apply(cdfs, 1, function(x)prod(x)) # exp(sum(log(x)))
    denom = mean(products)
    for(k in 1:length(others)){
      num = mean(products / cdfs[, k] * pdfs[, k])
      Ez[i, others[k]] = Ez[i, others[k]] - num / denom
    }
    # update for Ez[i, y[i]]
    Ez[i, y[i]] = sum(muZ[i, ]) - sum(Ez[i, -y[i]])
    lowerbound_yz = lowerbound_yz + log(denom)
  }
  return(list(Ez = Ez, lowerbound_yz = lowerbound_yz))
}

# for ordinal data
update_ordinal_latent_sum_j = function(j, X_list, Ev, Ew, n_levels_j){
  latent = Ev %*% Ew[[j]]
  out = sapply(1:n_levels_j, function(l){
    sum(latent[X_list[[j]] == l])
  })
  return(out)
}
update_ordinal_latent_sum = function(X_list, Ev, Ew, n_levels, type){
  out = list()
  for(j in which(type == "ordinal")){
    out[[j]] = update_ordinal_latent_sum_j(j, X_list, Ev, Ew, n_levels[[j]])
  }
  return(out)
}
