matrix_list_sum = function(lst) Reduce("+", lst)
trace = function(mat) sum(diag(mat))
logdet = function(mat) determinant(mat, logarithm=TRUE)$modulus
