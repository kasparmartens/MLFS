library(randomForest)

generate_missing_data = function(X, prop){
  N = nrow(X[[1]])
  M = length(X)
  if(sum(prop) == 0) return(list(X = X, missing = missing))
  for(j in 1:length(prop)){
    if(prop[j]>0){
      selected_individuals = sample(1:N, floor(prop[j]*N))
      X[[j]][selected_individuals, ] = NA
    }
  }
  return(X)
}

imputation_experiment = function(prop_mis_train, prop_mis_test, trainy, trainX, testy, testX, R, max_iter = 1000, burnin=500){
  type = rep("gaussian", length(trainX))
  
  df = foreach(i = 1:length(prop_mis_train), .combine = "rbind") %do% {
    Xmis = generate_missing_data(trainX, prop_mis_train[[i]])
    Xtestmis = generate_missing_data(testX, prop_mis_test[[i]])
    
    MLFSobj = MLFS_mcmc(trainy, Xmis, testy, Xtestmis, type, R, max_iter=max_iter, burnin = burnin, verbose=FALSE, impute=TRUE)
    pred_acc1 = MLFSobj$pred_acc_test
    pred_acc2 = pred_available_cases_MLFS(trainy, Xmis, testy, Xtestmis, type, max_iter)
    
    X_imputed = MLFS_mcmc_impute(trainy, Xmis, type, R, n_imputations = 100, max_iter=max_iter, burnin = burnin, verbose=FALSE)$U_trace
    pred_acc3 = MLFS_with_missing_views(trainy, X_imputed, testy, Xtestmis, type, R, max_iter=max_iter, burnin = burnin, verbose=FALSE)$pred_acc_test
    
    pred_acc_rf1 = rfhelper(trainy, Xmis, testy, Xtestmis, impute=FALSE)
    pred_acc_rf2 = rfhelper(trainy, Xmis, testy, Xtestmis, impute=TRUE)
    
    data.frame(pred_acc1 = pred_acc1, pred_acc2 = pred_acc2, pred_acc3 = pred_acc3, 
               pred_acc_rf1 = pred_acc_rf1, pred_acc_rf2 = pred_acc_rf2, missing_pattern = i)
  
  }
  return(df)
}

rfhelper = function(trainy, trainX, testy, testX, impute = FALSE){
  Ntest = length(testy)
  ypred = rep(NA, Ntest)
  observed_patterns = do.call("cbind", lapply(testX, function(x)complete.cases(x)))
  unique_patterns = unique(observed_patterns)
  for(j in 1:nrow(unique_patterns)){
    selected_views = unique_patterns[j, ]
    subset_test = apply(observed_patterns, 1, function(x)sum(x != selected_views)==0)
    if(impute){
      X_rfimputed = rfImpute(do.call("cbind", trainX[selected_views]), factor(trainy), ntree=1000)[,-1]
      m = randomForest(X_rfimputed, factor(trainy), ntree=1000)
    } else{
      X_rf = do.call("cbind", trainX[selected_views])
      subset_train = complete.cases(X_rf)
      m = randomForest(X_rf[subset_train, ], factor(trainy[subset_train]), ntree=1000)
    }
    ypred[subset_test] = predict(m, newdata = do.call("cbind", testX[selected_views]))
  }
  pred_acc_test = mean(ypred == testy)
  return(pred_acc_test)
}

pred_available_cases_MLFS = function(trainy, trainX, testy, testX, type, max_iter = 1000){
  M = length(trainX)
  Ntest = length(testy)
  ypred = rep(NA, Ntest)
  observed_train = do.call("cbind", lapply(trainX, function(x)complete.cases(x)))
  observed_test = do.call("cbind", lapply(testX, function(x)complete.cases(x)))
  unique_patterns = unique(observed_test)
  for(j in 1:nrow(unique_patterns)){
    selected_views = unique_patterns[j, ]
    subset_train = apply(observed_train, 1, function(x)sum(x[x != selected_views] != TRUE)==0)
    subset_test = apply(observed_test, 1, function(x)sum(x != selected_views)==0)
    MLFSobj = MLFS_mcmc(trainy[subset_train], 
                        lapply(trainX[selected_views], function(mat)mat[subset_train, ]), 
                        testy[subset_test], 
                        lapply(testX[selected_views], function(mat)mat[subset_test, ]), 
                        type[selected_views], R, max_iter=max_iter, verbose=FALSE)
    ypred[subset_test] = MLFSobj$pred_test
  }
  acc = mean(ypred == testy)
  return(acc)
}
