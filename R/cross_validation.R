

require(data.table)
require(survival)

# cross validate models, return the best one 
fit_best_recurrent_cox_model = function(input,list_features_models = c('Nm'),n_folds = 10, seed = NA, shuffle_ids = TRUE, verbose = FALSE){
    
    input = data.table(input)
    # set seed
    if(!is.na(seed)){set.seed(seed)}
    
    # compute scores for all models in list_features_models
    all_scores = NULL
    for(features in list_features_models){
        scores = cross_validate_formula(input, features = features, n_folds, seed, shuffle_ids, verbose)
        all_scores = cbind(all_scores,scores)
    }
    rownames(all_scores) = 1:n_folds
    colnames(all_scores) = list_features_models

    # summarize CV scores
    mean_scores = colMeans(all_scores)
    best_score = min(mean_scores)
    worst_score = max(mean_scores)

    # print cv scores
    if(verbose){.print_cv_scores(list_features_models,mean_scores,best_score,worst_score)}

    # fit best model on the entire training set 
    if(verbose){cat("\n")}
    if(verbose){print('Fitting best model on the entire training set')}
    best_features = list_features_models[mean_scores == best_score]
    best_formula = as.formula(paste("Surv(start,stop,status)~ ",best_features,' + cluster(id)',sep = ""))
    best_model = coxph(best_formula ,data = input)

    return(best_model)
}


# cross validate a single model, return the scores
cross_validate_formula = function(input, features = 'Nm', n_folds = 10, seed = NA, shuffle_ids = TRUE, verbose = FALSE){
  input = data.table(input)
  # set seed
  if(!is.na(seed)){set.seed(seed)}
  # randomly shuffle ids
  if(shuffle_ids){ids = sample(unique(input$id))}
  #Create n_folds folds with equal number of individuals (not observations!)
  folds <- cut(seq(1,length(ids)),breaks=n_folds,labels=FALSE)
  
  # build formula
  str_formula = paste("Surv(start,stop,status)~ ",features,' + cluster(id)',sep = "")
  formula = as.formula(str_formula)
  if(verbose){
    print(paste('Cross validating formula:',str_formula))
    pb <- txtProgressBar(min = 0, max = n_folds, style = 3)
  }
  
  #Perform  cross validation
  scores = NULL
  for(i in 1:n_folds){
    if(verbose){setTxtProgressBar(pb, i)}
    #Segment ids by fold 
    validIndexes <- which(folds==i,arr.ind=TRUE)
    valid_ids <- ids[validIndexes]
    train_ids <- ids[-validIndexes]
    # Split train, valid
    train = input[id %in% train_ids]
    valid = input[!id %in% valid_ids]
    
    # fit
    model = coxph(formula ,data = train)
    # predict, evaluate
    fold_scores = c()
    
    # evaluate martingale residuals at validation events
    valid_model = coxph(model$formula,data = valid, init = model$coefficients, iter.max = 0)
    res = residuals(valid_model, type = "martingale")
    score = mean(abs(res)) # mean absolute residual 
    
    # store score
    scores = c(scores,score)
  }
  if(verbose){cat("\n")}
  return(scores)
}


# print CV scores
.print_cv_scores = function(list_features_models,mean_scores,best_score,worst_score){
  cat('\n\n')
  print('****************** Cross validation mean absolute Martingale residual ******************')
  cat('\n')
  print('|------------------------|--------------------------------------|')
  print('|   MEAN ABS. RESIDUAL   |               FORMULA                |')
  print('|------------------------|--------------------------------------|')
  for(i in 1:length(list_features_models)){
    score = mean_scores[i]
    if(score != best_score & score != worst_score){
      print(paste('|       ',round(mean_scores[i],3),'          | ',list_features_models[i]))
    }else if (score == best_score){
      print(paste('|       ',round(mean_scores[i],3),'          | ',list_features_models[i], '                 <------ BEST MODEL'))
    }else{
      print(paste('|       ',round(mean_scores[i],3),'          | ',list_features_models[i], '                 <------ WORST MODEL'))
    }
  }
}