rotation_f = function(vecQ, R, Eww, Evv_sum, M, N, C, d){
  Q = matrix(vecQ, R, R)
  QQinv = solve(Q %*% t(Q))
  temp = 0
  for(j in 1:M){
    temp = temp + sum(d[j] * log(diag(t(Q) %*% Eww[[j]] %*% Q)))
  }
  value = 0.5*temp + 0.5*trace(QQinv %*% Evv_sum) + (N+C-sum(d)) * logdet(Q)
  return(value)
}

rotation_grad = function(vecQ, R, Eww, Evv_sum, M, N, C, d){
  Q = matrix(vecQ, R, R)
  QQinv = solve(Q %*% t(Q))
  grad = (N+C-sum(d)) * solve(t(Q)) - QQinv %*% Evv_sum %*% QQinv %*% Q
  for(j in 1:M){
    # denom=Q(:,r2)'*Eww{t}*Q(:,r2)
    denominators = diag(t(Q) %*% Eww[[j]] %*%  Q) + 1e-16
    temp = 2*(Eww[[j]] %*% Q) / rep(denominators, each=R)
    grad = grad + 0.5*d[j]*temp
  }
  vec_grad = as.numeric(grad)
  return(vec_grad)
}
