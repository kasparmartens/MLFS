matrix_list_sum = function(lst) Reduce("+", lst)
trace = function(mat) sum(diag(mat))
logdet = function(mat) determinant(mat, logarithm=TRUE)$modulus
vec = function(mat)as.numeric(mat)
inversenormal = function(x) qnorm(ppoints(length(x))[rank(x)])
logit = function(p) log((p/(1-p+1e-16)+1e-16))


plot_heatmap = function(W, ...){
  W_concat = do.call("cbind", W)
  # colnames are needed to add the annotation color
  if(is.null(colnames(W_concat))) colnames(W_concat) = paste0("test", 1:ncol(W_concat))
  annotation = data.frame(view = factor(rep(1:length(W), sapply(W, ncol))))
  rownames(annotation) = colnames(W_concat)
  pheatmap(abs(W_concat), cluster_cols = F, annotation_col = annotation, ...)
}


check_convergence = function(likelihood, i, convergence_threshold = 1e-6){
  if(i == 1){
    return(FALSE)
  }
  else{
    relative_change = (likelihood[i-1] - likelihood[i]) / likelihood[i-1]
    return(ifelse(relative_change < convergence_threshold, TRUE, FALSE))
  }
}

impute_with_median = function(mat){
  colmedians = apply(mat, 2, median, na.rm=T)
  missing_rows = !complete.cases(mat)
  missing_entries = is.na(mat)
  out = mat
  for(row in which(missing_rows)){
    which_cols = missing_entries[row, ]
    out[row, ] = colmedians[which_cols]
  }
  out
}

impute_test_with_median = function(testX, trainX){
  for(j in 1:length(testX)){
    colmedians = apply(trainX[[j]], 2, median, na.rm=T)
    missing_rows = !complete.cases(testX[[j]])
    missing_entries = is.na(testX[[j]])
    for(row in which(missing_rows)){
      which_cols = missing_entries[row, ]
      testX[[j]][row, ] = colmedians[which_cols]
    }
  }
  return(testX)
}

traceplot = function(mat){
  if(class(mat) == "matrix"){
    par(mfrow = c(0.5*ncol(mat), 2))
    for(j in 1:ncol(mat)){
      plot(mat[, j], type="l", ylab = "")
    }
  }else{
    plot(mat, type="l", ylab="")
  }
  layout(1)
}

reorder_rows_one_view = function(U, reordering, which_view){
  U[[which_view]] = U[[which_view]][reordering, ]
  return(U)
}

compute_loglikelihood = function(U, V, W, gamma, y, z, beta, rho, M, N){
  loglik = 0
  for(j in 1:M){
    mu_U = V %*% W[[j]]
    # loglik = loglik + sum(dnorm(U[[j]], mu_U, 1/sqrt(gamma[j]), log=TRUE))
    loglik = loglik + logdnormCpp(U[[j]], V %*% W[[j]], 1/gamma[j]) 
  }
  # p(beta)
  loglik = loglik + sum(dnorm(beta, 0, 1/sqrt(rho), log=TRUE))
  # p(V)
  loglik = loglik + sum(dnorm(V, 0, 1, log=TRUE))
  # p(z | v, beta)
  z_mu = V %*% beta
  loglik = loglik + sum(dnorm(z, z_mu, 1, log=TRUE))
  # p(y|z)
  prob = pnorm(z, 0, 1)
  loglik = loglik + sum(ifelse(y==2, log(prob), log(1-prob)))
  return(loglik)
}

compute_loglikelihood_marginal = function(U, V, W, gamma, y, z, rho, M, N){
  loglik = 0
  for(j in 1:M){
    # loglik = loglik + sum(dnorm(U[[j]], mu_U, 1/sqrt(gamma[j]), log=TRUE))
    loglik = loglik + logdnormCpp(U[[j]], V %*% W[[j]], 1/gamma[j]) 
  }
  # p(z | v) = int p(beta) p(z | v, beta) dbeta
  loglik = loglik + dmvnorm(z, rep(0, N), diag(N) + 1/rho * V %*% t(V), log = TRUE)
  # p(y|z)
  prob = pnorm(z, 0, 1)
  loglik = loglik + sum(ifelse(y==2, log(prob), log(1-prob)))
  return(loglik)
}

update_state_mat = function(label_state_matrices, current_indexes){
  for(j in 1:length(label_state_matrices)){
    label_state_matrices[[j]] = current_indexes[[j]]
  }
}

helper_switch_pairs = function(switched_pairs, N){
  prev_k = length(switched_pairs)
  if(prev_k > 0) x = c(prev_k-1, prev_k, prev_k+1) else x = c(0, 1)
  proposal_k = sample(x, 1, prob=dgeom(x, prob=0.5))
  if(proposal_k == 0) return(list(switched_pairs = list(), changed_indexes = c()))
  if(proposal_k < prev_k){
    # switch two indexes back
    index = sample(seq_along(switched_pairs), 1)
    return(list(switched_pairs = switched_pairs[-index], changed_indexes = c()))
  } else if ((proposal_k == prev_k) & (prev_k > 1)){
    # pick two pairs and switch them around
    two_indexes = sample(seq_along(switched_pairs), 2)
    first = switched_pairs[[two_indexes[1]]]
    second = switched_pairs[[two_indexes[2]]]
    switched_pairs[[two_indexes[1]]] = c(first[1], second[1])
    switched_pairs[[two_indexes[2]]] = c(first[2], second[2])
    changed_indexes = c(first, second)
  } else{
    # introduce a new reordering
    new_pair = sample(setdiff(1:N, unlist(switched_pairs)), 2)
    switched_pairs[[length(switched_pairs)+1]] = new_pair
    changed_indexes = new_pair
  }
  return(list(switched_pairs = switched_pairs, changed_indexes = changed_indexes))
}
