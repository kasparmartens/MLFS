source("helpers.R")
source("par_updates.R")

MLFS = function(y, X_list, type, R, max_iter=10){
  M = length(X_list)
  d = sapply(X_list, ncol)
  
  # hyperpars
  aGamma=1e-3
  bGamma=1e-3
  aAlpha=1e-2
  bAlpha=1e-2
  
  # initialise cutpoints
  #n_levels = 3
  #g_list = lapply(1:M, function(i)seq(-5, 5, n_levels))
  
  # initialise U
  Eu = X_list
  
  # initialise V
  allU = do.call("cbind", Eu)
  svd_obj = svd(allU)
  Ev = svd_obj$u[, 1:R]
  sigma_V = 1e-3*diag(R)
  # Evv_i = lapply(1:nrow(Ev), function(i)Ev[i, ] %*% t(Ev[i, ]) + sigma_V) # for each i there is E(v_i v_i)
  Evv_sum = t(Ev) %*% Ev + sigma_V
  
  # initialise W
  temp = solve(t(Ev) %*% Ev, t(Ev))
  Ew = lapply(Eu, function(U)temp %*% U)
  Eww = lapply(Ew, function(W) W %*% t(W))
  sigma_W = lapply(Ew, function(x)1e-3*diag(R))
  
  # initialise alpha
  Ealpha = matrix(0, M, R)
  Elogalpha = matrix(0, M, R)
  for(j in 1:M){
    aTmp = aAlpha + 0.5*d[j]
    bTmp = bAlpha + 0.5*diag(Eww[[j]])
    Ealpha[j, ] = aTmp / bTmp
    Elogalpha[j, ] = digamma(aTmp) - log(bTmp)
  }
  
  # initalise gamma
  Egamma = rep(0, M)
  Eloggamma = rep(0, M)
  for(j in 1:M){
    a_tilde = (aGamma + d[j] * N/2)
    b_tilde = bGamma + 0.5*norm(Eu[[j]] - Ev %*% Ew[[j]], "F")**2
    Egamma[j] = a_tilde / b_tilde
    Eloggamma[j] = digamma(a_tilde) - log(b_tilde)
  }
  
  # initialise z and beta
  Ez = matrix(0, N, C)
  Ebeta = matrix(0, R, C)
  Ebetabeta = Ebeta %*% t(Ebeta)
  
  for(iter in 1:max_iter){
    # optimise g
    
    ### update U
    
    ### update V
    sigmainv_V = update_V_sigmainv(Egamma, Eww, Ebetabeta, R, M)
    sigma_V = solve(sigmainv_V)
    lowerbound_V = N*logdet(sigma_V)
    for(i in 1:N){
      mu_tmp = update_V_mu_individual_i(i, Eu, Ew, Egamma, Ez, Ebeta, M, type)
      Ev[i, ] = mu_tmp %*% sigma_V
    }
    # Evv_i = lapply(1:nrow(Ev), function(i)t(Ev[i, ]) %*% Ev[i, ] + sigma_V)
    # Evv_sum = matrix_list_sum(Evv_i)
    Evv_sum = t(Ev) %*% Ev + sigma_V
    
    ### update alpha
    for(j in 1:M){
      aTmp = aAlpha + 0.5*d[j]
      bTmp = bAlpha + 0.5*diag(Eww[[j]])
      Ealpha[j, ] = aTmp / bTmp
      Elogalpha[j, ] = digamma(aTmp) - log(bTmp)
    }
    
    ### update W
    sigma_W = update_sigma_W(Ealpha, Egamma, Evv_sum, M)
    # make sure that assigned Ew is of the same structure/dims as old one!
    Ew = update_mu_W(Eu, Ev, Ew, Egamma, type, d, M, sigma_W)
    Eww = lapply(1:M, function(j)Ew[[j]] %*% t(Ew[[j]]) + d[j] * sigma_W[[j]])
    lowerbound_W = 0
    for(j in 1:M){
      if(type[j] == "gaussian"){
        lowerbound_W = lowerbound_W + d[j]*logdet(sigma_W[[j]])
      }
    }
    
    ### summarize sufficient statistics for ordinal data
    
    ### update gamma
    for(j in 1:M){
      a_tilde = (aGamma + d[j] * N/2)
      b_tilde = bGamma + 0.5*trace(Eww[[j]] %*% Evv_sum) - 
        trace(t(Ew[[j]]) %*% t(Ev) %*% Eu[[j]]) + 0.5*sum(Eu[[j]]**2)
      Egamma[j] = a_tilde / b_tilde
      Eloggamma[j] = digamma(a_tilde) - log(b_tilde)
    }
    
    ### update tau
    
    ### rotate
    
    # rotation lower bound, do not rotate at the moment
    Q = diag(R)
    QQinv = solve(Q %*% t(Q))
    temp = 0
    for(j in 1:M){
      temp = temp + sum(d[j] * log(diag(t(Q) %*% Eww[[j]] %*% Q)))
    }
    neg_lowerbound_vw = 0.5*temp + 0.5*trace(QQinv %*% Evv_sum) + (N+C-sum(d)) * logdet(Q)
    
    
    ### update beta
    rho = 1e-2
    sigma_beta = solve(Evv_sum + rho * diag(R))
    Ebeta = sigma_beta %*% t(Ev) %*% Ez
    Ebetabeta = Ebeta %*% t(Ebeta) + sigma_beta
    
    ### classification
    Ez = Ev %*% Ebeta
    ypred = apply(Ez, 1, which.max)
    
    ### update z
    n_samples = 1000
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
      products = apply(cdfs, 1, function(x)exp(sum(log(x))))
      denom = mean(products)
      for(k in 1:length(others)){
        num = mean(products / cdfs[, k] * pdfs[, k])
        Ez[i, others[k]] = Ez[i, others[k]] - num / denom
      }
      lowerbound_yz = lowerbound_yz + log(denom)
    }
    
    lowerbound_vw = - neg_lowerbound_vw + 0.5*lowerbound_W + 0.5*lowerbound_V
    lowerbound = lowerbound_yz + lowerbound_vw
    
    cat(sprintf("lower bound:\t %1.3f\n", lowerbound))
  }
}

