matrix_list_sum = function(lst) Reduce("+", lst)
trace = function(mat) sum(diag(mat))
logdet = function(mat) determinant(mat, logarithm=TRUE)$modulus
vec = function(mat)as.numeric(mat)
inversenormal = function(x) qnorm(ppoints(length(x))[rank(x)])
logit = function(p) log(p/(1-p))

check_convergence = function(likelihood, i, convergence_threshold = 1e-6){
  if(i == 1){
    return(FALSE)
  }
  else{
    relative_change = (likelihood[i-1] - likelihood[i]) / likelihood[i-1]
    return(ifelse(relative_change < convergence_threshold, TRUE, FALSE))
  }
}
