rsparsematrix = function(m, n, binary_vector){
  A = matrix(rnorm(m*n), m, n)
  return(A * binary_vector)
}


H = rbind(c(1, 1, 1, 1), 
          c(1, 1, 1, 0), 
          c(1, 0, 1, 1), 
          c(1, 1, 0, 1), 
          c(1, 0, 0, 1), 
          c(0, 1, 1, 0))

R = nrow(H)
N = 100
# dimensionalities
d = c(5, 5, 5, 5)

V = matrix(rnorm(N*R), N, R)

W_list = lapply(1:ncol(H), function(j)rsparsematrix(R, d[j], H[, j]))
type = rep("gaussian", length(d))

W = do.call("cbind", W_list)

# Generate U
U_list = lapply(W_list, function(W){
  V %*% W
})

# Generate X
X_list = list()
for(j in 1:length(type)){
  if(type[j] == "gaussian"){
    X_list[[j]] = U_list[[j]]
  }
  else if(type[j] == "ordinal"){
    u = U_list[[j]]
    temp = u
    for(k in 1:(length(g)-1)){
      condition = (u > g[k]) & (u < g[k+1])
      temp[condition] = k
    }
    X_list[[j]] = temp
  }
}

C = 3
beta = matrix(rnorm(R*C), R, C)
z = V %*% beta
y = apply(z, 1, which.max)

