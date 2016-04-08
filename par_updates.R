update_V_sigmainv = function(Egamma, Eww, Ebetabeta, R, M){
  sigmainv = diag(R) + matrix_list_sum(Ebetabeta)
  for(j in 1:M){
    sigmainv = sigmainv + Egamma[j] * Eww[[j]]
  }
  return(sigmainv)
}

update_V_mu_individual_i = function(i, Eu, Ew, Egamma, Ez, Ebeta, M, type){
  mu = Ez[i, ] %*% Ebeta
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
