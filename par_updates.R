update_V_sigmainv = function(Egamma, Eww, Ebetabeta, R, M){
  sigmainv = diag(R) + matrix_list_sum(Ebetabeta)
  for(j in 1:M){
    sigmainv = sigmainv + Egamma[j] * Eww[[j]]
  }
  return(sigmainv)
}

update_V_mu_individual_i = function(i, Eu, Ew, Egamma, Ez, Ebeta, M, type){
  mu = Ez[i, ] %*% t(Ebeta)
  for(j in 1:M){
    if(type[j] == "gaussian"){
      mu = mu + Egamma[j] * Eu[[j]][i, ] %*% t(Ew[[j]])
    }
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

update_mu_W = function(Eu, Ev, Ew, Egamma, type, d, M, sigmaW){
  mu = list()
  for(j in 1:M){
    mu[[j]] = sigmaW[[j]] * Egamma[j]
    if(type[j] == "gaussian"){
      temp = sapply(1:ncol(Ew[[j]]), function(which_col){
        t(Ev) %*% Eu[[j]][, which_col]
      })
    }
    mu[[j]] = mu[[j]] %*% temp
  }
  return(mu)
}

update_Z = function(y, N, C, Ez, sigma_beta, Ebetabeta, rho, n_samples = 1000){
  lowerbound_yz = -0.5*rho*trace(Ebetabeta)+0.5*C*logdet(sigma_beta)
  for(i in 1:N){
    others = setdiff(1:C, y[i])
    z_diff = Ez[i, y[i]] - Ez[i, others]
    # update for Ez[i, y[i]]
    Ez[i, y[i]] = sum(Ez[i, ]) - sum(Ez[i, -y[i]])
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
    lowerbound_yz = lowerbound_yz + log(denom)
  }
  return(list(Ez = Ez, lowerbound_yz = lowerbound_yz))
}
