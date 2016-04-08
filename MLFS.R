MLFS = function(y, X_list, type, R, max_iter=10, rotate=TRUE){
  N = length(y)
  M = length(X_list)
  d = sapply(X_list, ncol)
  C = max(y)
  if(sum(!(y %in% 1:max(y))) > 0) stop("y must have labels 1, 2, ..., C")
  
  # hyperpars
  aGamma=1e-3
  bGamma=1e-3
  aAlpha=1e-2
  bAlpha=1e-2
  
  rho = 1e-2
  
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
  Ealpha = update_alpha(Ealpha, Eww, d, M, aAlpha, bAlpha)
  
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
  
  pred_acc_train = rep(NA, max_iter)
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
    Ealpha = update_alpha(Ealpha, Eww, d, M, aAlpha, bAlpha)
    
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
    if(rotate){
      optimres = optim(par = as.numeric(diag(R)), fn = rotation_f, gr = rotation_grad, 
                       method="BFGS", control = list(maxit=1000), 
                       R = R, Eww = Eww, Evv_sum = Evv_sum, M = M, N = N, C = C, d = d)
      Q = matrix(optimres$par, R, R)
      Qinv = solve(Q)
      Ev = Ev %*% t(Qinv)
      Evv_sum = Qinv %*% Evv_sum %*% t(Qinv)
      for(j in 1:M){
        Ew[[j]] = t(Q) %*% Ew[[j]]
        Eww[[j]] = t(Q) %*% Eww[[j]] %*% Q
        # same update for all types
        sigma_W[[j]] = t(Q) %*% sigma_W[[j]] %*% Q        
      }
      Ealpha = update_alpha(Ealpha, Eww, d, M, aAlpha, bAlpha)
      neg_lowerbound_vw = optimres$value
    } 
    else{
      neg_lowerbound_vw = rotation_f(vecQ, R, Eww, Evv_sum, M, N, C, d)
    }
    
    
    ### update beta
    sigma_beta = solve(Evv_sum + rho * diag(R))
    Ebeta = sigma_beta %*% t(Ev) %*% Ez
    Ebetabeta = Ebeta %*% t(Ebeta) + sigma_beta
    
    ### classification
    Ez = Ev %*% Ebeta
    ypred = apply(Ez, 1, which.max)
    pred_acc_train[iter] = mean(ypred == y)
    
    ### update z
    obj = update_Z(y, N, C, Ez, sigma_beta, Ebetabeta, rho, n_samples = 1000)
    Ez = obj$Ez
    lowerbound_yz = obj$lowerbound_yz
    
    lowerbound_vw = - neg_lowerbound_vw + 0.5*lowerbound_W + 0.5*lowerbound_V
    lowerbound = lowerbound_yz + lowerbound_vw
    
    cat(sprintf("lower bound:\t %1.3f\n", lowerbound))
  }
  print(Ew)
  cat("Prediction accuracies (train):", pred_acc_train)
  
}

