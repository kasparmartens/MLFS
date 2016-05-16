library(mvtnorm)
library(truncnorm)
library(pheatmap)

MLFS_mcmc = function(y, X_list, y_test, X_test, type, R, max_iter=10, d_sim = 5, verbose = TRUE, burnin = 500, label_switching = FALSE){
  N = length(y)
  M = length(X_list)
  d = ifelse(type != "similarity", sapply(X_list, ncol), d_sim)
  C = max(y)
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
  rho = 100
  
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
    if(type[j] == "gaussian") U[[j]] = X_imputed[[j]]
  }
  U_initial = U
  
  # initialise z, beta
  z = rep(0, N)
  ztest = rep(0, Ntest)
  beta = rep(0, R)
  beta_trace = matrix(NA, max_iter, R)
  pred_trace = matrix(NA, max_iter, N)
  pred_test_trace = matrix(NA, max_iter, Ntest)
  gamma_trace = list()
  W_trace = list()
  loglik_trace = rep(NA, max_iter)
  z_trace = matrix(NA, max_iter, N)
  
  current_indexes = lapply(1:M, function(x)1:N)
  label_state_matrices = lapply(1:M, function(x)matrix(0, N, N))
  switch_indicator = rep(0, N)
  
  if(verbose) cat("Running Variational Bayes for 20 iterations to obtain starting values ... \n")
  MLFSinit = MLFS(y, X_imputed, type, R, max_iter=20, verbose=FALSE)
  if(verbose) cat("\tFinished VB. Now starting Gibbs sampling\n")
  
  V = MLFSinit$Ev
  W = MLFSinit$Ew
  alpha = MLFSinit$Ealpha
  gamma = MLFSinit$Egamma
  V_mean = matrix(0, N, R)
  Vtest_mean = matrix(0, Ntest, R)
  
  for(iter in 1:max_iter){
    temp = update_H_and_W(U, V, W, alpha, gamma, pi, d, M, R, iter)
    H = temp$H
    W = temp$W
    
    # update pi
    H_rowsums = rowSums(H)
    pi = rbeta(M, 1 + H_rowsums, 1 + R - H_rowsums)
    
    ### update alpha
    for(j in 1:M){
      temp0 = rgamma(d[j], a0, b0)
      temp1 = rgamma(d[j], a0 + 0.5*d[j], b0 + 0.5*rowSums(W[[j]]**2))
      alpha[j, ] = ifelse(H[j, ] == 1, temp1, temp0)
    }
    
    ### update gamma
    for(j in 1:M){
      mis = missing_values[[j]]
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
      mu_i = (z[i]*beta + mu_tmp) %*% sigma_V
      V[i, ] = mu_i + random[i, ]
    }
    
    randomtest = rmvnorm(Ntest, rep(0, R), sigma_V)
    for(k in 1:Ntest){
      mu_tmp = update_V_mu_individual_i(k, Utest, W, gamma, M)
      mu_i = (ztest[k]*beta + mu_tmp) %*% sigma_V
      Vtest[k, ] = mu_i + randomtest[k, ]
    }
    
    V_mean = V_mean + 1/max_iter*V
    Vtest_mean = Vtest_mean + 1/max_iter*Vtest
    
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
    
    #     ### regression
    #     sigma_beta_inv = rho*diag(R) + t(V)%*%V
    #     sigma_beta = solve(sigma_beta_inv)
    #     beta_mu = sigma_beta %*% t(V) %*% y
    #     a_lambda = 1e-3 + 0.5*N
    #     residuals = y - mean(y)
    #     b_lambda = as.numeric(1e-3 + 0.5*(sum(residuals**2) - t(beta_mu)%*%sigma_beta_inv%*%beta_mu))
    #     lambda = rgamma(1, a_lambda, b_lambda)
    #     beta = as.numeric(t(beta_mu) + rmvt(1, b_lambda/a_lambda*sigma_beta, 2*a_lambda))
    #     pred = V %*% beta
    #     pred_trace[iter, ] = pred
    
    ### For classification follow Albert & Chib 1993
    # sample beta
    sigma_beta = solve(1/rho * diag(R) + t(V) %*% V)
    mu_beta = sigma_beta %*% t(V) %*% z
    beta = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
    beta_trace[iter, ] = beta
    
    # sample z
    z_mu = V %*% beta
    subset = (y == 2)
    if(sum(subset)>0) z[subset] = rtruncnorm(sum(subset), a=0, mean=z_mu[subset], sd=1)
    subset = (y == 1)
    if(sum(subset)>0) z[subset] = rtruncnorm(sum(subset), b=0, mean=z_mu[subset], sd=1)
    
    z_mu_test = Vtest %*% beta
    if(iter == 1){
      ztest = z_mu_test      
    } else{
      subset = (ytest == 2)
      if(sum(subset)>0) ztest[subset] = rtruncnorm(sum(subset), a=0, mean=z_mu_test[subset], sd=1)
      subset = (ytest == 1)
      if(sum(subset)>0) ztest[subset] = rtruncnorm(sum(subset), b=0, mean=z_mu_test[subset], sd=1)
    }
    
    z_trace[iter, ] = z
    
    pred_train = rbinom(N, 1, pnorm(z)) + 1
    pred_trace[iter, ] = pred_train
    
    pred_test = rbinom(Ntest, 1, pnorm(ztest)) + 1
    ytest = pred_test
    pred_test_trace[iter, ] = pred_test
    
    loglik_trace[iter] = compute_loglikelihood(U, V, W, gamma, y, z, beta, rho, M, N)
    
    if((iter > burnin) & label_switching){
      for(kk in 1:100){
        V_proposal = V
        z_proposal = z
        proposal_view = 1 #sample(1:M, 1)
        proposal_indexes = current_indexes[[proposal_view]]
        
        two_indexes = sample((1:N)[switch_indicator %in% c(0, proposal_view)], 2)
        proposal_indexes[two_indexes] = proposal_indexes[rev(two_indexes)]
        
        U_proposal = reorder_rows_one_view(U_initial, proposal_indexes, proposal_view)
        # update V for the two proposal indexes
        for(i in two_indexes){
          mu_tmp = update_V_mu_individual_i(i, U_proposal, W, gamma, M)
          mu_i = (z[i]*beta + mu_tmp) %*% sigma_V
          V_proposal[i, ] = mu_i + random[i, ]
        }
        # update z
        z_mu = V_proposal %*% beta
        subset = (y == 2)
        if(sum(subset)>0) z_proposal[subset] = rtruncnorm(sum(subset), a=0, mean=z_mu[subset], sd=1)
        subset = (y == 1)
        if(sum(subset)>0) z_proposal[subset] = rtruncnorm(sum(subset), b=0, mean=z_mu[subset], sd=1)
        # update beta
        sigma_beta = solve(1/rho * diag(R) + t(V_proposal) %*% V_proposal)
        mu_beta = sigma_beta %*% t(V_proposal) %*% z_proposal
        beta_proposal = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
        
        l_proposal = compute_loglikelihood(U_proposal, V_proposal, W, gamma, y, z_proposal, beta_proposal, rho, M, N)
        accept_prob = min(1, exp(l_proposal - loglik_trace[iter]))
        if(runif(1) < accept_prob){
          cat("\n")
          switch_indicator[two_indexes] = proposal_view
          U = U_proposal
          current_indexes[[proposal_view]] = proposal_indexes
          V = V_proposal
          beta = beta_proposal
          z = z_proposal
          cat("Accept prob", accept_prob, "Changed view", proposal_view, "indexes", two_indexes, "\n")
        }
      }
      # Update label state matrices
      for(j in 1:M){
        label_state_matrices[[j]] = label_state_matrices[[j]] + as.matrix(table(1:N, current_indexes[[j]]))
      }
    }
    
    if((iter %% 100 == 0)){
      # plot_heatmap(W, main=sprintf("Iter %s", iter), cluster_rows = FALSE)
      if(verbose) cat(sprintf("Iter %d. Prediction accuracies: train %1.3f, test %1.3f \n", iter, mean(pred_train == y), mean(pred_test == y_test)))
      # cat("iter", iter, "pred acc train:", mean(pred_train == y), "pred acc test:", mean(pred_test == y_test), "\n")
      # print(H)
    }
  }
  pred_train = round(apply(pred_trace[-c(1:burnin), ], 2, mean))
  pred_test = round(apply(pred_test_trace[-c(1:burnin), ], 2, mean))
  pred_acc_train = mean(pred_train == y)
  pred_acc_test = mean(pred_test == y_test)
  
  return(list(pred_train = pred_train, pred_acc_train = pred_acc_train, 
              pred_test_trace = pred_test_trace[-c(1:burnin), ],  
              pred_test = pred_test, pred_acc_test = pred_acc_test, 
              beta_trace = beta_trace, z_trace = z_trace, W_trace = W_trace,
              loglik_trace = loglik_trace, label_state_mat = label_state_matrices))
}

### update H
# update_H_and_W = function(U, V, W, alpha, gamma, pi, d, M, R){
#   H = matrix(0, M, R)
#   for(j in 1:M){
#     s2 = 1/alpha[j, ] + gamma[j] * diag(t(V) %*% V) # sapply(1:ncol(V), function(r)t(V[, r]) %*% V[, r])
#     Ures = U[[j]] - V %*% W[[j]]
#     mu = gamma[j] * t(Ures) %*% V %*% diag(1 / s2)
#     z = logit(pi[j]) + 0.5*d[j]*log(s2 * alpha[j, ]) + 0.5*diag(t(mu) %*% mu) / s2
#     acceptance_probs = 1 / (1 + exp(-z))
#     print(acceptance_probs)
#     u = runif(R)
#     H[j, ] = ifelse(u < acceptance_probs, 1, 0)
#     W_fill = mu + rnorm(R*d[j])*sqrt(rep(s2, each=d[j]))
#     W[[j]] = matrix(W_fill * rep(H[j, ], d[j]), R, d[j])
#   }
#   return(list(H = H, W = W))W[[j]][r, ]
# }

update_H_and_W = function(U, V, W, alpha, gamma, pi, d, M, R, iter){
  H = matrix(0, M, R)
  for(j in 1:M){
    for(r in 1:R){
      Ures = U[[j]] - V[, -r, drop=FALSE] %*% W[[j]][-r, , drop=FALSE]
      lambda = (alpha[j, r] + gamma[j] * sum(V[, r]**2))
      mu = gamma[j] / lambda * t(Ures) %*% V[, r]
      z = logit(pi[j]) - 0.5*d[j]*log(lambda * alpha[j, r]) + 0.5*lambda*sum(mu**2)
      acceptance_prob = 1 / (1 + exp(-z))
      # if(iter %% 100 == 0) cat("lambda", lambda, "alpha", alpha[j, r], "gamma", gamma[j], "\n")
      # if(iter %% 10 == 0) cat("accept prob", acceptance_prob,"\n")
      u = runif(1)
      H[j, r] = ifelse(u < acceptance_prob, 1, 0)
      W_fill = mu + rnorm(d[j])*sqrt(1/lambda)
      W[[j]][r, ] = W_fill * H[j, r]
    }
  }
  return(list(H = H, W = W))
}






