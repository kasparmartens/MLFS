library(mvtnorm)
library(truncnorm)
library(pheatmap)

MLFS_mcmc = function(y, X_list, y_test, X_test, type, R, max_iter=10, d_sim = 5, verbose = TRUE, 
                     burnin = 500, label_switching = FALSE, marginal_sampler = FALSE, 
                     impute = TRUE, X_full = NULL, X_full_test = NULL, proposal_view = 1, new_method=TRUE){
  N = length(y)
  M = length(X_list)
  d = ifelse(type != "similarity", sapply(X_list, ncol), d_sim)
  C = max(y)
  #if(sum(!(y %in% 1:max(y))) > 0) stop("y must have labels 1, 2, ..., C")
  if(sum(sapply(X_list, class) != "matrix") > 0) stop("X_list must contain matrices only!")
  
  # if(sum(type != "gaussian") > 0) stop("MCMC version has been implemented for Gaussian views only")
  missing_values = lapply(X_list, function(x)!complete.cases(x))
  missing_values_test = lapply(X_test, function(x)!complete.cases(x))
  # imputation for starting values
  X_imputed = lapply(X_list, impute_with_median)
  X_test_imputed = impute_test_with_median(X_test, X_imputed)

  temp = 0
  for(j in 1:M){
    for(mis in which(missing_values[[j]])){
      temp = temp + mean(abs(X_full[[j]][mis, ] - X_imputed[[j]][mis, ]))
    }
  }
  imputation_baseline = temp / sum(unlist(missing_values))

  
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
  n_levels = list()
  for(j in 1:M){
    if(type[j] == "gaussian") U[[j]] = X_list[[j]]
    if(type[j] == "ordinal"){
      U[[j]] = scale(X_list[[j]], center=TRUE, scale=FALSE)
      g_list[[j]] = as.numeric(quantile(U[[j]], c(0.33, 0.66)))
      n_levels[[j]] = max(X_list[[j]])+1
    }
  }
  U_initial = U
  
  # initialise z, beta
  z = rep(0, N)
  ztest = rep(0, Ntest)
  beta = rep(0, R)
  mu_beta = rep(0, R)
  beta_trace = matrix(NA, max_iter, R)
  pred_trace = matrix(NA, max_iter, N)
  pred_test_trace = matrix(NA, max_iter, Ntest)
  gamma_trace = list()
  W_trace = list()
  H_trace = list()
  loglik_trace = rep(NA, max_iter)
  z_trace = matrix(NA, max_iter, N)
  # imputation accuracy
  imputation_acc = rep(0, max_iter)
  imputation_acc_test = rep(0, max_iter)
  
  current_indexes = lapply(1:M, function(x)1:N)
  label_state_matrices = lapply(1:M, function(x)matrix(0, N, N))
  # switch_indicator = rep(0, N)
  switched_pairs = list()
  
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
    mu = z %*% t(beta)
    for(j in 1:M){
      mu = mu + gamma[j] * U[[j]] %*% t(W[[j]])
    }
    V = mu %*% sigma_V + random
    
    
    randomtest = rmvnorm(Ntest, rep(0, R), sigma_V)
    mu = ztest %*% t(beta)
    for(j in 1:M){
      mu = mu + gamma[j] * Utest[[j]] %*% t(W[[j]])
    }
    Vtest = mu %*% sigma_V + randomtest
    
    if(iter > burnin){
      V_mean = V_mean + 1/(max_iter-burnin)*V
      Vtest_mean = Vtest_mean + 1/(max_iter-burnin)*Vtest
    }
    
    # draw new g
    for(j in which(type == "ordinal")){
      g_list[[j]][1] = runif(1, max(U[[j]][X_list[[j]] == 0]), min(min(U[[j]][X_list[[j]] == 1]), g_list[[j]][2]))
      g_list[[j]][2] = runif(1, max(max(U[[j]][X_list[[j]] == 1], g_list[[j]][1])), min(U[[j]][X_list[[j]] == 2]))
#       mean_pred = V%*%W[[j]]
#       current_likelihood = truncnorm_likelihood(X_list[[j]], U[[j]], g_list[[j]], mean_pred, gamma[j])
#       for(rep in 1:10){
#         lambda = rnorm(n_levels[[j]]-1, 0, 0.1)
#         g_proposal = cumsum(exp(lambda))
#         proposal_likelihood = truncnorm_likelihood(X_list[[j]], U[[j]], g_proposal, mean_pred, gamma[j])
#         # cat("proposal_ll", proposal_likelihood, "\n")
#         if(runif(1) < exp(proposal_likelihood - current_likelihood)){
#           g_list[[j]] = g_proposal
#           current_likelihood = proposal_likelihood
#           cat("proposed g", g_proposal, ": accepted\n")
#         }# else cat("proposed g", g_proposal, ": rejected\n")
#       }
    }
    if(iter%%100==0)print(g_list[[2]])
    # update U
    for(j in which(type == "ordinal")){
      U[[j]] = compute_U(X_list[[j]], g_list[[j]], mean_pred = V %*% W[[j]], gamma[j])
    }
    
    ### impute missing values in U
    if(impute){
      imp_acc = 0
      imp_acc_test = 0
      for(j in 1:M){
        sigma_U = 1 / gamma[j] * diag(d[j])
        for(mis in which(missing_values[[j]])){
          mu_U = V[mis, ] %*% W[[j]]
          U[[j]][mis, ] = rnorm(d[j], mu_U, 1/sqrt(gamma[j]))
          imp_acc = imp_acc + mean(abs(X_full[[j]][mis, ] - U[[j]][mis, ]))
        }
        for(mis in which(missing_values_test[[j]])){
          mu_U = Vtest[mis, ] %*% W[[j]]
          Utest[[j]][mis, ] = rnorm(d[j], mu_U, 1/sqrt(gamma[j]))
          imp_acc_test = imp_acc_test + mean(abs(X_full_test[[j]][mis, ] - Utest[[j]][mis, ]))
        }
      }
      imputation_acc[iter] = imp_acc / sum(unlist(missing_values))
      imputation_acc_test[iter] = imp_acc_test / sum(unlist(missing_values_test))
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
    sigma_beta = solve(rho * diag(R) + t(V) %*% V)
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
    # loglik_trace[iter] = compute_loglikelihood_marginal(U, V, W, gamma, y, z, rho, M, N)
    
    if((iter > burnin) & label_switching){
      if(new_method){
        for(kk in 1:100){
          proposal_indexes = 1:N
          temp = helper_switch_pairs(switched_pairs, N)
          switched_pairs_proposal = temp$switched_pairs
          proposal_indexes[unlist(switched_pairs_proposal)] = proposal_indexes[unlist(lapply(switched_pairs_proposal, rev))]
          
          U_proposal = reorder_rows_one_view(U_initial, proposal_indexes, proposal_view)
          V_proposal = V
          z_proposal = z
          
          for(i in temp$changed_indexes){
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
          sigma_beta = solve(rho * diag(R) + t(V_proposal) %*% V_proposal)
          mu_beta = sigma_beta %*% t(V_proposal) %*% z_proposal
          beta_proposal = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
          l_proposal = compute_loglikelihood(U_proposal, V_proposal, W, gamma, y, z_proposal, beta_proposal, rho, M, N)

          accept_prob = min(1, exp(l_proposal - loglik_trace[iter]))
          if(runif(1) < accept_prob){
            U = U_proposal
            current_indexes[[proposal_view]] = proposal_indexes
            V = V_proposal
            beta = beta_proposal
            z = z_proposal
            switched_pairs = switched_pairs_proposal
            loglik_trace[iter] = l_proposal
            # cat("Accept prob", accept_prob, "Changed view", proposal_view, "indexes", two_indexes, "\n")
          }
        }
      }
      else{ 
        for(kk in 1:100){
          V_proposal = V
          z_proposal = z
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
          
          if(marginal_sampler){
            sigma_beta = solve(rho * diag(R) + t(V_proposal) %*% V_proposal)
            h[two_indexes] = diag(V_proposal[two_indexes, ] %*% sigma_beta %*% t(V_proposal[two_indexes, ]))
            w = h / (1-h)
            S = sigma_beta %*% t(V_proposal)
            mu_beta = as.numeric(sigma_beta %*% t(V) %*% z)
            obj = probit_rcpp_helper(y, mu_beta, V, w, sqrt(w+1), z, S, N, R)
            z_proposal = obj$z
            mu_beta = obj$mu_beta
            beta_proposal = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
            l_proposal = compute_loglikelihood_marginal(U_proposal, V_proposal, W, gamma, y, z_proposal, rho, M, N)
          } else{
            # update z
            z_mu = V_proposal %*% beta
            subset = (y == 2)
            if(sum(subset)>0) z_proposal[subset] = rtruncnorm(sum(subset), a=0, mean=z_mu[subset], sd=1)
            subset = (y == 1)
            if(sum(subset)>0) z_proposal[subset] = rtruncnorm(sum(subset), b=0, mean=z_mu[subset], sd=1)
            # update beta
            sigma_beta = solve(rho * diag(R) + t(V_proposal) %*% V_proposal)
            mu_beta = sigma_beta %*% t(V_proposal) %*% z_proposal
            beta_proposal = as.numeric(rmvnorm(1, mu_beta, sigma_beta))
            l_proposal = compute_loglikelihood(U_proposal, V_proposal, W, gamma, y, z_proposal, beta_proposal, rho, M, N)
          }
          accept_prob = min(1, exp(l_proposal - loglik_trace[iter]))
          if(runif(1) < accept_prob){
            switch_indicator[two_indexes] = proposal_view
            U = U_proposal
            current_indexes[[proposal_view]] = proposal_indexes
            V = V_proposal
            beta = beta_proposal
            z = z_proposal
            loglik_trace[iter] = l_proposal
            # cat("Accept prob", accept_prob, "Changed view", proposal_view, "indexes", two_indexes, "\n")
          }
        }
        
      }
      
      # Update label state matrices
      for(j in 1:M){
        label_state_matrices[[j]] = label_state_matrices[[j]] + as.matrix(table(1:N, current_indexes[[j]]))
      }
    }
    
    if((iter %% 100 == 0)){
      H_trace[[iter/100]] = H
      W_trace[[iter/100]] = W
      # plot_heatmap(W, main=sprintf("Iter %s", iter), cluster_rows = FALSE)
      if(verbose) cat(sprintf("Iter %d. Prediction accuracies: train %1.3f, test %1.3f \n", iter, mean(pred_train == y), mean(pred_test == y_test)))
    }
  }
  pred_train = round(apply(pred_trace[-c(1:burnin), ], 2, mean))
  pred_test = round(apply(pred_test_trace[-c(1:burnin), ], 2, mean))
  pred_acc_train = mean(pred_train == y)
  pred_acc_test = mean(pred_test == y_test)
  
  label_state_matrices = lapply(label_state_matrices, function(x)x / (max_iter - burnin))
  
  
  
  return(list(pred_train = pred_train, pred_acc_train = pred_acc_train, 
              pred_test_trace = pred_test_trace[-c(1:burnin), ],  
              pred_test = pred_test, pred_acc_test = pred_acc_test, 
              beta_trace = beta_trace, z_trace = z_trace, W_trace = W_trace, H_trace = H_trace, 
              sigma_V = sigma_V, V_mean = V_mean, Vtest_mean = Vtest_mean, 
              imputation_acc = mean(imputation_acc[-c(1:burnin)]), imputation_acc_test = mean(imputation_acc_test[-c(1:burnin)]), imputation_baseline = imputation_baseline, 
              loglik_trace = loglik_trace, label_state_mat = label_state_matrices))
}


update_H_and_W = function(U, V, W, alpha, gamma, pi, d, M, R, iter){
  H = matrix(0, M, R)
  Vsq_colsums = colSums(V**2)
  for(j in 1:M){
    for(r in 1:R){
      Ures = U[[j]] - V[, -r, drop=FALSE] %*% W[[j]][-r, , drop=FALSE]
      lambda = alpha[j, r] + gamma[j] * Vsq_colsums[r]
      mu = gamma[j] / lambda * t(Ures) %*% V[, r]
      z = logit(pi[j]) - 0.5*d[j]*(log(lambda) + log(alpha[j, r])) + 0.5*lambda*sum(mu**2)
      acceptance_prob = 1 / (1 + exp(-z))
      # if(iter %% 100 == 0) cat("lambda", lambda, "alpha", alpha[j, r], "gamma", gamma[j], "\n")
      # if(iter %% 10 == 0) cat("accept prob", acceptance_prob,"\n")
      u = runif(1)
      if(u < acceptance_prob){
        H[j, r] = 1
        W[[j]][r, ] = mu + rnorm(d[j])*sqrt(1/lambda)
      } else{
        H[j, r] = 0
        W[[j]][r, ] = 0
      }
    }
  }
  return(list(H = H, W = W))
}



fill_truncnorm = function(X, g){
  out = X
  out[X == 0] = runif(sum(X==0), -5, g[1])
  out[X == 1] = runif(sum(X==1), g[1], g[2])
  out[X == length(g)] = runif(sum(X==length(g)), g[length(g)], 5)
  return(out)
}

truncnorm_likelihood = function(X, U, g, mean_pred, gamma){
  res = sum(log(dtruncnorm(U[X == 0], a=-Inf, b=g[1], mean=mean_pred[X==0], sd=1)))
  res = res + sum(log(dtruncnorm(U[X == 1], a=g[1], b=g[2], mean=mean_pred[X==1], sd=1)))
  res = res + sum(log(dtruncnorm(U[X == length(g)], a=g[length(g)], b=Inf, mean=mean_pred[X==2], sd=1)))
  return(res)
}

compute_U = function(X, g, mean_pred, gamma){
  out = X
  temp0 = rtruncnorm(sum(X==0), -Inf, g[1], mean=mean_pred[X==0], sd=1/sqrt(gamma))
  temp1 = rtruncnorm(sum(X==1), g[1], g[2], mean=mean_pred[X==1], sd=1/sqrt(gamma))
  temp2 = rtruncnorm(sum(X==length(g)), g[length(g)], Inf, mean=mean_pred[X==length(g)], sd=1/sqrt(gamma))
  out[X == 0] = temp0
  out[X == 1] = temp1
  out[X == length(g)] = temp2
  if(sum(is.na(temp0))>0){
    print(mean_pred[X==0])
    print(g)
  }
  if(sum(is.na(temp1))>0){
    print(mean_pred[X==1])
    print(g)
  }
  if(sum(is.na(temp2))>0){
    print(mean_pred[X==2])
    print(g)
  }
  return(out)
}
