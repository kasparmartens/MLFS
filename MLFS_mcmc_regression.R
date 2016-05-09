MLFS_mcmc_regression = function(y, X_list, y_test, X_test, type, R, max_iter=10, d_sim = 5, verbose = TRUE, burnin = 100, label_switching = FALSE){
  N = length(y)
  M = length(X_list)
  d = ifelse(type != "similarity", sapply(X_list, ncol), d_sim)
  #if(sum(!(y %in% 1:max(y))) > 0) stop("y must have labels 1, 2, ..., C")
  if(sum(sapply(X_list, class) != "matrix") > 0) stop("X_list must contain matrices only!")
  
  if(sum(type != "gaussian") > 0) stop("MCMC version has been implemented for Gaussian views only")
  missing_values = lapply(X_list, function(x)!complete.cases(x))
  missing_values_test = lapply(X_test, function(x)!complete.cases(x))
  # imputation for variational inference
  X_imputed = lapply(X_list, impute_with_median)
  X_test_imputed = lapply(X_test, impute_with_median)
  
  Ntest = nrow(X_test[[1]])
  Utest = X_test_imputed
  Vtest = matrix(0, Ntest, R)
  
  a0 = 0.01
  b0 = 0.01
  a_gamma = 0.01
  b_gamma = 0.01
  rho = 0.01
  
  # alpha, gamma
  alpha = matrix(100, M, R)
  gamma = rep(100, M)
  
  # initialise pi, H, W
  pi = rbeta(M, 1, 1)
  H = matrix(runif(M*R) < pi, M, R)
  W = list()
  for(j in 1:M){
    fill = rep(H[j, ], d[j])*rnorm(d[j]*R)
    W[[j]] = matrix(fill, R, d[j])
  }
  
  # initialise V
  V = matrix(rnorm(N*R), N, R)
  
  # initalise U
  U = list()
  for(j in 1:M){
    if(type[j] == "gaussian") U[[j]] = ifelse(is.na(X_list[[j]]), 0, X_list[[j]])
  }
  U_initial = U
  
  # initialise z, beta
  beta = rep(0, R)
  pred_test = rep(0, Ntest)
  beta_trace = matrix(NA, max_iter, R)
  pred_trace = matrix(NA, max_iter, N)
  pred_test_trace = matrix(NA, max_iter, N)
  gamma_trace = list()
  W_trace = list()
  loglik_U_trace = rep(NA, max_iter)
  loglik_y_trace = rep(NA, max_iter)
  z_trace = matrix(NA, max_iter, N)
  
  current_indexes = 1:N
  
  if(verbose) cat("Running Variational Bayes for 20 iterations to obtain starting values ... \n")
  MLFSinit = MLFS_regression(y, X_imputed, type, R, max_iter=20, verbose=FALSE)
  if(verbose) cat("\tFinished VB. Now starting Gibbs sampling\n")
  
  V = MLFSinit$Ev
  W = MLFSinit$Ew
  alpha = MLFSinit$Ealpha
  gamma = MLFSinit$Egamma
  
  for(iter in 1:max_iter){
    temp = update_H_and_W(U, V, W, alpha, gamma, pi, d, M, R, iter)
    H = temp$H
    W = temp$W
    
    # update pi
    H_rowsums = rowSums(H)
    pi = rbeta(M, 0.01 + H_rowsums, 0.01 + R - H_rowsums)
    
    ### update alpha
    for(j in 1:M){
      temp0 = rgamma(d[j], a0, b0)
      temp1 = rgamma(d[j], a0 + 0.5*d[j], b0 + 0.5*sum(W[[j]]**2))
      alpha[j, ] = ifelse(H[j, ] == 1, temp1, temp0)
    }
    
    ### update gamma
    for(j in 1:M){
      temp = norm(U[[j]] - V %*% W[[j]], "F")
      gamma[j] = rgamma(1, a_gamma + 0.5*d[j]*N, b_gamma + 0.5*temp**2)
    }
    
    ### sample V
    sigmainv_V = diag(R) + beta %*% t(beta) + matrix_list_sum(lapply(1:M, function(j){
      gamma[j] * W[[j]] %*% t(W[[j]])
    }))
    sigma_V = solve(sigmainv_V + 1e-6)
    
    random = rmvnorm(N, rep(0, R), sigma_V)
    for(i in 1:N){
      mu_tmp = update_V_mu_individual_i(i, U, W, gamma, M)
      mu_i = (y[i]*beta + mu_tmp) %*% sigma_V
      V[i, ] = mu_i + random[i, ]
    }
    
    randomtest = rmvnorm(Ntest, rep(0, R), sigma_V)
    for(k in 1:Ntest){
      mu_tmp = update_V_mu_individual_i(k, Utest, W, gamma, M)
      mu_i = (pred_test[k]*beta + mu_tmp) %*% sigma_V
      Vtest[k, ] = mu_i + randomtest[k, ]
    }
    
    ### impute missing values in U
    for(j in 1:M){
      sigma_U = 1 / gamma[j] * diag(d[j])
      for(mis in which(missing_values[[j]])){
        mu_U = V[mis, ] %*% W[[j]]
        U[[j]][mis, ] = rnorm(d[j], mu_U, 1/sqrt(gamma[j])) # rmvnorm(1, mu_U, sigma_U)
      }
      for(mis in which(missing_values_test[[j]])){
        mu_U = Vtest[mis, ] %*% W[[j]]
        Utest[[j]][mis, ] = rnorm(d[j], mu_U, 1/sqrt(gamma[j])) # rmvnorm(1, mu_U, sigma_U)
      }
    }
    
    #     loglik_U = 0
    #     for(j in 1:M){
    #       mu_U = V %*% W[[j]]
    #       sigma_U = 1 / gamma[j] * diag(d[j])
    #       for(i in 1:N){
    #         loglik_U = loglik_U + dmvnorm(U[[j]][i, ], mu_U[i, ], sigma_U, log=TRUE)
    #       }
    #     }
    #     loglik_U_trace[iter] = loglik_U

    ### regression
    sigma_beta_inv = rho*diag(R) + t(V)%*%V
    sigma_beta = solve(sigma_beta_inv)
    beta_mu = sigma_beta %*% t(V) %*% y
    a_lambda = 1e-3 + 0.5*N
    b_lambda = as.numeric(1e-3 + 0.5*(sum(y**2) - t(beta_mu)%*%sigma_beta_inv%*%beta_mu))
    lambda = rgamma(1, a_lambda, b_lambda)
    beta = as.numeric(t(beta_mu) + rmvt(1, b_lambda/a_lambda*sigma_beta, 2*a_lambda))
    beta_trace[iter, ] = beta
    
    pred_train = V %*% beta + rt(N, df = 2*a_lambda) * sqrt(diag(b_lambda/a_lambda*(diag(N)+V%*%sigma_beta%*%t(V))))
    pred_test = Vtest %*% beta + rt(Ntest, df = 2*a_lambda) * sqrt(diag(b_lambda/a_lambda*(diag(Ntest)+Vtest%*%sigma_beta%*%t(Vtest))))
    pred_trace[iter, ] = pred_train
    pred_test_trace[iter, ] = pred_test
    
    
    if(label_switching){
      proposal_indexes = switch_two_labels(current_indexes)$x
      if(iter==1) l_current = compute_likelihood(U, V, W, gamma, M, N)
      proposal_view = sample(1:M, 1)
      U_proposal = reorder_rows_one_view(U_initial, proposal_indexes, proposal_view)
      l_proposal = compute_likelihood(U_proposal, V, W, gamma, M, N)
      accept_prob = min(1, exp(l_proposal - l_current))
      if(runif(1) < accept_prob){
        U = U_proposal
        l_current = l_proposal
        cat("Changed view", proposal_view, "indexes", which(proposal_indexes != 1:N), "\n")
      }
    }
    
    
    
    if((iter %% 100 == 0)){
      # plot_heatmap(W, main=sprintf("Iter %s", iter), cluster_rows = FALSE)
      if(verbose) cat(sprintf("Iter %d. Prediction accuracies: train %1.3f, test %1.3f \n", iter, cor(pred_train, y)**2, cor(pred_test, y_test)**2))
    }
  }
  pred_train = (apply(pred_trace[-c(1:burnin), ], 2, mean))
  pred_test = (apply(pred_test_trace[-c(1:burnin), ], 2, mean))
  pred_acc_train = cor(pred_train, y)**2
  pred_acc_test = cor(pred_test, y_test)**2
  
  return(list(pred_train = pred_train, pred_acc_train = pred_acc_train, 
              pred_test_trace = pred_test_trace[-c(1:burnin), ],  
              pred_test = pred_test, pred_acc_test = pred_acc_test, 
              beta_trace = beta_trace, z_trace = z_trace, W_trace = W_trace, 
              loglik_U_trace = loglik_U_trace, loglik_y_trace = loglik_y_trace))
}