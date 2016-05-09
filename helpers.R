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

switch_two_labels = function(x){
  ind = sample(1:length(x), 2)
  
  before = sum(x[ind] == ind)
  after = sum(x[rev(ind)] == ind)
  prob = ifelse(after == 2, 0.9, ifelse(after == 1, 0.5, 0.1))
  x[ind] = x[rev(ind)]
  return(list(x = x, prob = prob))
}

reorder_rows_one_view = function(U, reordering, which_view){
  U[[which_view]] = U[[which_view]][reordering, ]
  return(U)
}

compute_likelihood = function(U, V, W, gamma, M, N){
  loglik_U = 0
  for(j in 1:M){
    mu_U = V %*% W[[j]]
    sigma_U = 1 / gamma[j] * diag(ncol(U[[j]]))
    for(i in 1:N){
      loglik_U = loglik_U + dmvnorm(U[[j]][i, ], mu_U[i, ], sigma_U, log=TRUE)
    }
  }
  return(loglik_U)
}
