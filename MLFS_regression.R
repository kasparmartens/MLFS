MLFS_regression = function(y, X_list, type, R, max_iter=10, rotate=TRUE, d_sim = 5, verbose = TRUE){
  N = length(y)
  M = length(X_list)
  d = ifelse(type != "similarity", sapply(X_list, ncol), d_sim)
  C = max(y)
  # if(sum(!(y %in% 1:max(y))) > 0) stop("y must have labels 1, 2, ..., C")
  if(sum(sapply(X_list, class) != "matrix") > 0) stop("X_list must contain matrices only!")
  
  # hyperpars
  aGamma=1e-3
  bGamma=1e-3
  aAlpha=1e-2
  bAlpha=1e-2
  
  rho = 1e-2
  scaling_const = 0.2
  G = 10
  
  Erho = 100
  Elambda = 100
  ypred = rep(mean(y), N)
  Ebeta = rep(0, R)
  Ebetabeta = Ebeta %*% t(Ebeta)
  
  # initialise U
  Eu = list()
  Euu = list()
  Euu_sum = list()
  for(j in 1:M){
    X = X_list[[j]]
    if(type[j] == "gaussian") Eu[[j]] = X
    if(type[j] == "ordinal") Eu[[j]] = (X - min(X)) / diff(range(X))*2*G-G
    if(type[j] == "similarity"){
      Euu[[j]] = list()
      svd_obj = svd(X)
      Eu[[j]] = svd_obj$u[, 1:d[j]] %*% diag(sqrt(svd_obj$d[1:d[j]]))
      for(i in 1:N){
        Euu[[j]][[i]] = Eu[[j]][i, ] %*% t(Eu[[j]][i, ]) + 1e-6*diag(d[j])
      }
      Euu_sum[[j]] = matrix_list_sum(Euu[[j]])
    }
  }
  
  # initialise V
  allU = do.call("cbind", Eu)
  svd_obj = svd(allU)
  Ev = svd_obj$u[, 1:R]
  sigma_V = 1e-3*diag(R)
  Evv_sum = t(Ev) %*% Ev + sigma_V
  
  # initialise W
  temp = solve(t(Ev) %*% Ev, t(Ev))
  Ew = lapply(Eu, function(U)temp %*% U)
  Eww = lapply(Ew, function(W) W %*% t(W))
  sigma_W = lapply(Ew, function(x)1e-3*diag(R))
  
  # initialise cutpoints for ordinal views
  g = list()
  n_levels = list()
  ordinal_counts = list()
  for(j in which(type == "ordinal")){
    n_levels_j = max(X_list[[j]])
    n_levels[[j]] = n_levels_j
    g[[j]] = seq(-G, G, length=n_levels_j+1)
    ordinal_counts[[j]] = sapply(1:n_levels_j, function(l)sum(X_list[[j]] == l))
  }
  ordinal_latent_sum = update_ordinal_latent_sum(X_list, Ev, Ew, n_levels, type)
  
  # initialise alpha
  Ealpha = matrix(0, M, R)
  Elogalpha = matrix(0, M, R)
  Ealpha = update_alpha(Ealpha, Eww, d, M, aAlpha, bAlpha)
  
  for(j in 1:M){
    if(type[j] == "ordinal"){
      for(l in 1:n_levels[[j]]){
        Eu[[j]][X_list[[j]] == l] = mean(g[[j]][l:(l+1)])
      }
    }
  }
  
  # initalise gamma
  Egamma = rep(0, M)
  Eloggamma = rep(0, M)
  for(j in 1:M){
    a_tilde = (aGamma + d[j] * N/2)
    b_tilde = bGamma + 0.5*norm(Eu[[j]] - Ev %*% Ew[[j]], "F")**2
    if(type[j] == "similarity"){
      a_tilde = a_tilde + 0.25*N*(N-1)
      temp = 0
      for(i in 1:(N-1)){
        residual = X_list[[j]][i, (i+1):N] - Eu[[j]][i, ] %*% t(Eu[[j]][(i+1):N, , drop=FALSE])
        temp = temp + sum(residual**2)
      }
      b_tilde = b_tilde + 0.5*scaling_const*temp
    }
    Egamma[j] = a_tilde / b_tilde
  }
  
  # initialise tau
  Etau = rep(0, M)
  for(j in which(type == "similarity")){
    Etau[j] = scaling_const*Egamma[j]
  }
  
  # initialise z and beta
  # Ez = matrix(0, N, C)
  # Ebeta = rep(NA, R)
  # Ebetabeta = Ebeta %*% t(Ebeta)
  
  pred_acc_train = rep(NA, max_iter)
  lowerbound = rep(NA, max_iter)
  for(iter in 1:max_iter){
    
    ### optimise g
    for(j in which(type == "ordinal")){
      constr = generate_constraints(n_levels[[j]], G)
      optimres = constrOptim(g[[j]], cutpoints_f, cutpoints_grad, constr$ui, constr$ci, 
                             ordinal_counts = ordinal_counts[[j]], ordinal_latent_sum = ordinal_latent_sum[[j]], 
                             Egamma_j = Egamma[[j]], n_levels = n_levels[[j]])
      g[[j]] = optimres$par
    }
    
    ### update U
    for(j in 1:M){
      if(type[j] == "ordinal"){
        for(l in 1:n_levels[[j]]){
          Eu[[j]][X_list[[j]] == l] = mean(g[[j]][l:(l+1)])
        }
      }
      if(type[j] == "similarity"){
        # sigma_U[[j]] = solve(Egamma[j] * diag(d[j]) + Etau[j] * Euu_sum[[j]])
        for(i in 1:N){
          Euu_sum_temp = matrix_list_sum(Euu[[j]][setdiff(1:N, i)])
          sigma_U = solve(Egamma[j] * diag(d[j]) + Etau[j] * Euu_sum_temp)
          temp = X_list[[j]][i, -i] %*% Eu[[j]][-i, ]
          Eu[[j]][i, ] = (Etau[j] * temp + Egamma[j] * Ev[i, ] %*% Ew[[j]]) %*% sigma_U
          Euu[[j]][[i]] = Eu[[j]][i, ] %*% t(Eu[[j]][i, ]) + sigma_U
        }
        Euu_sum[[j]] = matrix_list_sum(Euu[[j]])
      }
    }
    
    ### update V
    sigmainv_V = update_V_sigmainv(Egamma, Eww, R, M) + Ebetabeta
    sigma_V = solve(sigmainv_V)
    lowerbound_V = N*logdet(sigma_V)
    Ev = ypred %*% t(Ebeta)
    for(i in 1:N){
      mu_tmp = update_V_mu_individual_i(i, Eu, Ew, Egamma, M)
      Ev[i, ] = (Ev[i, ] + mu_tmp) %*% sigma_V
    }
    # Evv_i = lapply(1:nrow(Ev), function(i)t(Ev[i, ]) %*% Ev[i, ] + sigma_V)
    # Evv_sum = matrix_list_sum(Evv_i)
    Evv_sum = t(Ev) %*% Ev + N*sigma_V
    
    ### update alpha
    Ealpha = update_alpha(Ealpha, Eww, d, M, aAlpha, bAlpha)
    
    ### update W
    sigma_W = update_sigma_W(Ealpha, Egamma, Evv_sum, M)
    # make sure that assigned Ew is of the same structure/dims as old one!
    Ew = update_mu_W(Eu, Ev, Ew, Egamma, type, d, M, sigma_W)
    Eww = lapply(1:M, function(j)Ew[[j]] %*% t(Ew[[j]]) + d[j] * sigma_W[[j]])
    lowerbound_W = 0
    for(j in 1:M){
      if(type[j] %in% c("gaussian", "ordinal")){
        lowerbound_W = lowerbound_W + d[j]*logdet(sigma_W[[j]])
      }
    }
    
    ### update variable for ordinal data
    ordinal_latent_sum = update_ordinal_latent_sum(X_list, Ev, Ew, n_levels, type)
    
    ### update gamma
    obj = update_gamma(j, g, n_levels, ordinal_counts, Eu, Euu, Euu_sum, Ev, Ew, Eww, sigma_W, Evv_sum, X_list[[j]], M, N, d, type, aGamma, bGamma, scaling_const)
    Egamma = obj$a_tilde / obj$b_tilde
    lowerbound_xu = -sum(obj$a_tilde * log(obj$b_tilde))
    
    ### update tau 
    for(j in which(type == "similarity")){
      Etau[j] = scaling_const*Egamma[j]
    }
    
    ### rotate
    if(rotate){
      # it may be faster to initialise with the previous Q
      optimres = tryCatch({
        optim(par = as.numeric(diag(R)), fn = rotation_f, gr = rotation_grad, 
              method="BFGS", control = list(maxit=1000), 
              R = R, Eww = Eww, Evv_sum = Evv_sum, M = M, N = N, C = C, d = d)
      }, error = function(e){
        warning("Q %*% t(Q) is computationally singular. Switching to rotate=FALSE")
        list(par = diag(R), 
             value = rotation_f(as.numeric(diag(R)), R, Eww, Evv_sum, M, N, C, d))
      })
      Q = matrix(optimres$par, R, R)
      Qinv = solve(Q)
      Ev = Ev %*% t(Qinv)
      sigma_V = Qinv %*% sigma_V %*% t(Qinv)
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
      neg_lowerbound_vw = rotation_f(as.numeric(diag(R)), R, Eww, Evv_sum, M, N, C, d)
    }
    
    # update lowerbound XU
    for(j in 1:M){
      if(type[j] == "ordinal"){
        lowerbound_xu = lowerbound_xu + sum(ordinal_counts[[j]] * log(diff(g[[j]])))
      }
      else if(type[j] == "similarity"){
        lowerbound_xu = lowerbound_xu + sum(sapply(Euu[[j]], function(uu)0.5*logdet(uu)))
      }
    }
    
    sigma_beta = solve(Erho*diag(R) + Evv_sum)
    Ebeta = sigma_beta %*% t(Ev) %*% y
    a_rho = (1e-3 + 0.5*R)
    b_rho = as.numeric(1e-3 + 0.5*(Elambda * t(Ebeta)%*%Ebeta + trace(sigma_beta)))
    Erho = a_rho / b_rho
    
    ypred = Ev %*% Ebeta
    residuals = y - ypred
    a_lambda = (1e-3 + 0.5*N)
    b_lambda = as.numeric(1e-3 + 0.5*sum(residuals**2) + Erho * t(Ebeta) %*% Ebeta)
    Elambda = a_lambda / b_lambda

  
    pred_acc_train[iter] = cor(y, ypred)
    lowerbound_yz = 0.5*logdet(sigma_beta) - 0.5*Elambda*sum(residuals**2) - 0.5*sum(diag(Ev %*% sigma_beta %*% t(Ev))) +
      lgamma(a_lambda) - a_lambda*log(b_lambda) + a_lambda + lgamma(a_rho) - a_rho*log(b_rho)
    lowerbound_vw = - neg_lowerbound_vw + 0.5*lowerbound_W + 0.5*lowerbound_V
    lowerbound[iter] = lowerbound_yz + lowerbound_vw + lowerbound_xu
    
    # cat(sprintf("VW:\t %1.3f\tYZ:\t %1.3f\tXU:\t %1.3f\n", lowerbound_vw, lowerbound_yz, lowerbound_xu))
    # cat(sprintf("V:  %1.3f\tW: %1.3f\tnegVW: %1.3f\n", lowerbound_V, lowerbound_W, - neg_lowerbound_vw))
    if(verbose) cat(sprintf("lower bound:\t %1.3f\n", lowerbound[iter]))
    
    # if relative change in lower bound is less than 1e-6, then break
    #if(check_convergence(lowerbound, iter, 1e-6)) break
  }
  if(verbose) cat("Prediction accuracies (train):", pred_acc_train, "\n")
  return(list(Ebeta = Ebeta, Ew = Ew, Eww = Eww, sigma_W = sigma_W, sigma_V = sigma_V, 
              Egamma = Egamma, g = g, Etau = Etau, 
              lowerbound = lowerbound, 
              Eu_train = Eu, Euu_sum_train = Euu_sum, 
              n_levels = n_levels, type = type, R = R, d = d, d_sim = d_sim, 
              pred_acc_train = pred_acc_train[iter]))
}


pred_out_of_sample_regression = function(X_test, MLFSobj, X_test_train = NULL){
  
  type = MLFSobj$type
  Ntest = nrow(X_test[[1]])
  M = length(type)
  R = MLFSobj$R
  n_levels = MLFSobj$n_levels
  d = MLFSobj$d
  
  Ev = matrix(0, Ntest, R)
  
  # compute Eu
  Eu = list()
  for(j in 1:M){
    # initialise Eu (default value for gaussian view)
    if(type[j] %in% c("gaussian", "ordinal")) Eu[[j]] = X_test[[j]] else Eu[[j]] = matrix(0, Ntest, MLFSobj$d_sim)
    if(type[j] == "ordinal"){
      for(l in 1:n_levels[[j]]){
        Eu[[j]][X_test[[j]] == l] = mean(MLFSobj$g[[j]][l:(l+1)])
      }
    }
    if(type[j] == "similarity"){
      # sigma_U = solve(MLFSobj$Egamma[j] * diag(MLFSobj$d[j]) + MLFSobj$Etau[j] * MLFSobj$Euu_sum_train[[j]])
      for(i in 1){
        Eww = t(MLFSobj$Ew[[j]]) %*% MLFSobj$Ew[[j]] + diag(d[j])*sum(diag(MLFSobj$sigma_W[[j]]))
        sigma_U_inv = solve(1/MLFSobj$Egamma[j]*diag(d[j]) + Eww) + MLFSobj$Etau[j]*MLFSobj$Euu_sum_train[[j]]
        # sigma_U = solve(sigma_U_inv)
        # mu_U = MLFSobj$Etau[j] * X_test_train[[j]][i, ] %*% MLFSobj$Eu_train[[j]] %*% sigma_U
        Eu[[j]][i, ] = solve(sigma_U_inv, t(MLFSobj$Etau[j] * X_test_train[[j]][i, ] %*% MLFSobj$Eu_train[[j]]))
      }
    }
  }
  
  sigma_V = solve(update_V_sigmainv(MLFSobj$Egamma, MLFSobj$Eww, R, M))
  for(j in 1:M){
    for(i in 1:Ntest){
      mu_tmp = update_V_mu_individual_i(i, Eu, MLFSobj$Ew, MLFSobj$Egamma, M)
      Ev[i, ] = mu_tmp %*% sigma_V
    }
  }
  
  ypred = Ev %*% MLFSobj$Ebeta
  return(ypred)
}
