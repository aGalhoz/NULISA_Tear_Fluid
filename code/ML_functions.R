library('caret')
library('ranger')
library('glmnet')
library('pROC')
library('dplyr')
library('RobustRankAggreg')
library('matrixStats') # row standard deviation
devtools::install_github('Mikata-Project/ggthemr')
library('ggthemr')
ggthemr('greyscale', layout = "scientific")
library('viridis')
library('fgsea')
library('org.Hs.eg.db')


###       boot strap function       ###
# runs machine learning on data with selected algorithm ( 'lm' (linear regression) or
#                                                         'rf' (random forest) or
#                                                         'svm lin' Support Vector Machine + linear kernel)
#                                                         'svm rad' Support Vector Machine + radial kernel))
# number of folds for cross validation (default = 10) and boot strap loops (default = 1)
#
# returns list objects with models, importance, predictions (raw and probability) and indices (samples used for training)


runML = function(data_frame,algorithm, cv = 10, BS_number = 1){
  #check input
  stopifnot('Last argument (BS_number) is not a number.' = is.numeric(BS_number),
            'Number given for cross validation is not a number ,' = is.numeric(cv),
            'Data is not a dataframe.' = is.data.frame(data_frame),
            'status column of df not correct' = all(df_ml$status == 1 | df_ml$status == 0),
            'unknown algorithm. please select \'lm\',\'rf\',\'svm rad\' or \'svm lin\' ' = algorithm %in% c('lm','rf','svm rad','svm lin'))
  
  key_output = c('linear regressio (lasso)','random forest', 'Support Vector Machine (radial kernel)', 'Support Vector Machine (linear kernel)')
  names(key_output) = c('lm','rf','svm rad','svm lin')
  # Output parameters for user check
  print(paste0('Running a ',BS_number,'x boot strap with a ',cv,' fold cross validation. Algorithm is ', key_output[[algorithm]]))
  
  data_frame$status = as.factor(make.names(data_frame$status)) # caret needs prediction variable to have name; turns 0 -> X0 and 1 -> X1
  
  # create lists to save results from bs
  models = list()
  importance = list()
  predictions = list()
  indices = list()
  predictions_raw = list()
  smp_size = floor(0.8 * nrow(data_frame)) # use 80% of data for training, 20% for testing
  
  #set the right parameters for each algo
  if(algorithm == 'lm'| algorithm == 'svm rad' | algorithm == 'svm lin'){
    ctrl = trainControl(method="cv",   
                        number = cv,        
                        summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                        classProbs=TRUE,savePredictions = TRUE)
  }
  else if(algorithm == 'rf'){
    
    n_features = ncol(data_frame[ ,!names(data_frame) == 'status'])
    ctrl = trainControl(method = "cv", 
                        number = cv, 
                        search = 'grid',classProbs = TRUE, savePredictions = TRUE, summaryFunction=twoClassSummary )
  }
  
  # run machine learning
  for(i in 1:BS_number){
  
    train_ind = createDataPartition(data_frame$status, p = 0.8, list = FALSE)
    train = data_frame[train_ind, ]
    test = data_frame[ -train_ind,!names(data_frame) == "status"]
    
    if (algorithm == 'lm'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "glmnet", family = "binomial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
                    tuneGrid=expand.grid(
                      .alpha=1, # alpha 1 == lasso
                      .lambda=seq(0, 100, by = 0.1))
      )
    }
    
    else if (algorithm == 'rf'){
      tunegrid = expand.grid(
        .mtry = c(2, 3, 4, 7, 11, 17, 27, floor(sqrt(n_features)), 41, 64, 99, 154, 237, 367, 567, 876),
        .splitrule = c("extratrees","gini"),
        .min.node.size = c(1,2,3,4,5)
      )
      
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    tuneGrid = tunegrid, 
                    method = "ranger",  tuneLength = 15, metric = "ROC",
                    num.trees = 500,
                    trControl = ctrl,
                    importance = 'impurity'
      )
    }
    
    else if (algorithm == 'svm lin'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmLinear", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
    
    else if (algorithm == 'svm rad'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmRadial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
    
    models[[i]] = model
    importance[[i]] = varImp(model)
    predictions[[i]] = predict(model, newdata = test,type = "prob")
    predictions_raw[[i]] = predict(model, newdata = test,type = "raw")
    indices[[i]] = train_ind
    print(paste0('finshed loop #',i)) # keep track of what is happening
  }
  return_list = list(models,importance,predictions,predictions_raw,indices)
  names(return_list) = c('models','importance','predictions','predictions_raw','indices')
  return(return_list)
}


runML_balanced = function(data_frame,algorithm, cv = 10, BS_number = 1,
                          condition_1, condition_2){
  #check input
  stopifnot('Last argument (BS_number) is not a number.' = is.numeric(BS_number),
            'Number given for cross validation is not a number ,' = is.numeric(cv),
            'Data is not a dataframe.' = is.data.frame(data_frame),
            'status column of df not correct' = all(df_ml$status == 1 | df_ml$status == 0),
            'unknown algorithm. please select \'lm\',\'rf\',\'svm rad\' or \'svm lin\' ' = algorithm %in% c('lm','rf','svm rad','svm lin'))
  
  key_output = c('linear regressio (lasso)','random forest', 'Support Vector Machine (radial kernel)', 'Support Vector Machine (linear kernel)')
  names(key_output) = c('lm','rf','svm rad','svm lin')
  # Output parameters for user check
  print(paste0('Running a ',BS_number,'x boot strap with a ',cv,' fold cross validation. Algorithm is ', key_output[[algorithm]]))
  
  data_frame$status = as.factor(make.names(data_frame$status)) # caret needs prediction variable to have name; turns 0 -> X0 and 1 -> X1
  
  # create lists to save results from bs
  models = list()
  importance = list()
  predictions = list()
  indices = list()
  indices_test = list()
  predictions_raw = list()
  smp_size = floor(0.8 * nrow(data_frame)) # use 80% of data for training, 20% for testing
  
  #set the right parameters for each algo
  if(algorithm == 'lm'| algorithm == 'svm rad' | algorithm == 'svm lin'){
    ctrl = trainControl(method="cv",   
                        number = cv,        
                        summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                        classProbs=TRUE,savePredictions = TRUE)
  }
  else if(algorithm == 'rf'){
    
    n_features = ncol(data_frame[ ,!names(data_frame) == 'status'])
    ctrl = trainControl(method = "cv", 
                        number = cv, 
                        search = 'grid',classProbs = TRUE, savePredictions = TRUE, summaryFunction=twoClassSummary )
  }
  
  # run machine learning
  for(i in 1:BS_number){
    
    # get indices per condition
    cond1_idx = which(data_frame$status == "X0")
    cond2_idx = which(data_frame$status == "X1")
    
    n_samples = min(length(cond1_idx),length(cond2_idx))
    
    n_train = floor(0.8 * n_samples)
    n_test = n_samples - n_train
    
    # random balanced selection of conditions
    train_cond1 <- sample(cond1_idx, n_train)
    test_cond1 <- setdiff(sample(cond1_idx, n_samples), train_cond1)
    
    train_cond2 <- sample(cond2_idx, n_train)
    test_cond2 <- setdiff(sample(cond2_idx, n_samples), train_cond2)
    
    train_ind = c(train_cond1,train_cond2)
    test_ind = c(test_cond1,test_cond2)
    
    train = data_frame[train_ind, ]
    test = data_frame[test_ind,!names(data_frame) == "status"]
    
    if (algorithm == 'lm'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "glmnet", family = "binomial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
                    tuneGrid=expand.grid(
                      .alpha=1, # alpha 1 == lasso
                      .lambda=seq(0, 100, by = 0.1))
      )
    }
    
    else if (algorithm == 'rf'){
      tunegrid = expand.grid(
        .mtry = c(2, 3, 4, 7, 11, 17, 27, floor(sqrt(n_features)), 41, 64, 99, 154, 237, 367, 567, 876),
        .splitrule = c("extratrees","gini"),
        .min.node.size = c(1,2,3,4,5)
      )
      
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    tuneGrid = tunegrid, 
                    method = "ranger",  tuneLength = 15, metric = "ROC",
                    num.trees = 500,
                    trControl = ctrl,
                    importance = 'impurity'
      )
    }
    
    else if (algorithm == 'svm lin'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmLinear", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
    
    else if (algorithm == 'svm rad'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmRadial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
    
    models[[i]] = model
    importance[[i]] = varImp(model)
    predictions[[i]] = predict(model, newdata = test,type = "prob")
    predictions_raw[[i]] = predict(model, newdata = test,type = "raw")
    indices[[i]] = train_ind
    indices_test[[i]] = test_ind
    print(paste0('finshed loop #',i)) # keep track of what is happening
  }
  return_list = list(models,importance,predictions,predictions_raw,indices,indices_test)
  names(return_list) = c('models','importance','predictions','predictions_raw','indices','indices_test')
  return(return_list)
}


###       calculate averaged ROC curve        ###
# plots and returns data frame for averaged ROC curve from data and result from runML (indices, predictions (raw and probable))
calculateROC = function(list_from_ML, data_frame, plot_path = FALSE){
  # check input
  stopifnot('Data is not a dataframe.' = is.data.frame(data_frame),
            'First argument is not a list (result of runML).' = inherits(list_from_ML, 'list') & all(names(list_from_ML) == c('models','importance','predictions','predictions_raw','indices')),
            'Please give a valid directory to save the file' = ifelse(plot_path == FALSE, TRUE, dir.exists(paste(strsplit(plot_path, split = '/')[[1]][-length(strsplit(plot_path, split = '/')[[1]])], collapse = '/'))))
  
  # extract needed lists from ML data
  indices = list_from_ML[['indices']]
  predictions = list_from_ML[['predictions']]
  predictions_raw = list_from_ML[['predictions_raw']]
  
  # create one data frame to compare predicted vs actual status
  predictions_compare = data.frame()
  auc_values = c()   # store auc values here
  
  for(i in seq_along(indices)){
    actuals = data_frame$status[-indices[[i]]]
    probs   = predictions[[i]]$X1
    
    # skip if only one class in test
    if(length(unique(actuals)) < 2){
      warning(paste("Run", i, "skipped: only one class present in test set"))
      next
    }
    
    # compute auc with pROC instead of glmnet:::auc (safer)
    roc_obj = pROC::roc(actuals, probs, direction = "<")
    auc_val = pROC::auc(roc_obj)
    auc_values = c(auc_values, auc_val)
    
    adding = data.frame(Class = actuals, 
                        predicted = probs, 
                        resample = paste0("run ", i), 
                        auc = as.numeric(auc_val))
    predictions_compare = rbind(predictions_compare, adding)
  }
  
  rocc = data.frame()
  for(run in unique(predictions_compare$resample)){
    onefold = dplyr::filter(predictions_compare, resample == run)
    if(length(unique(onefold$Class)) < 2){
      warning(paste("Run", run, "skipped: only one class in test set"))
      next
    }
    
    roc_run = pROC::roc(onefold$Class, onefold$predicted, direction = "<")
    rocc = rbind(rocc,
                 data.frame(Sp = roc_run$specificities,
                            Sn = roc_run$sensitivities,
                            n = seq_along(roc_run$sensitivities)))
  }
  
  # aggregate
  Sp = aggregate(Sp ~ n, rocc, mean)$Sp
  Sn = aggregate(Sn ~ n, rocc, mean)$Sn
  
  if (length(unique(predictions_compare$resample)) > 1) {
    errorSp = 1.96 * aggregate(Sp ~ n, rocc, sd)$Sp 
    errorSn = 1.96 * aggregate(Sn ~ n, rocc, sd)$Sn 
  } else {
    # no variability with a single bootstrap → set CI = 0
    errorSp = rep(0, length(Sp))
    errorSn = rep(0, length(Sn))
  }
  
  plotci = data.frame(Sp, Sn, errorSp, errorSn)
  plotci <- plotci[!(plotci$Sp == 1 & !(plotci$Sn %in% c(0,1))), ]
  
  # plot if path given
  if(plot_path != FALSE){
    pdf(plot_path,paper="a4r", width = 11, height = 8)
    print(
      ggplot(plotci, aes(x=(1-Sp), y=Sn)) + 
        geom_line(color = "darkblue") + 
        geom_ribbon(aes(x = 1 - Sp,
                        ymin = pmax(Sn - errorSn, 0),
                        ymax = pmin(Sn + errorSn, 1)),
                    fill = "#B2B2B2", alpha = 0.5) +
        theme_bw() +
        labs(x = "1 - Specificity",
             y = "Sensitivity",
             title = "Mean ROC curve with 95 % CI") +
        annotate("text", x = 0.4, y = 0.75, 
                 label = paste("mean AUC: ", round(mean(auc_values),2),
                               "\u00B1 ", round(sd(auc_values),2))) +
        geom_abline(slope = 1, intercept = 0, color="darkgrey", alpha = 0.3) +
        xlim(-0.1,1.1)
    )
    dev.off()
  }
  
  return(list(roc_data = plotci, auc_values = auc_values))
}

calculateROC_balanced = function(list_from_ML, data_frame, plot_path = FALSE,
                                 min_per_condition = 2){
  # check input
  stopifnot('Data is not a dataframe.' = is.data.frame(data_frame),
            'First argument is not a list (result of runML).' = inherits(list_from_ML, 'list') & all(names(list_from_ML) == c('models','importance','predictions','predictions_raw','indices','indices_test')),
            'Please give a valid directory to save the file' = ifelse(plot_path == FALSE, TRUE, dir.exists(paste(strsplit(plot_path, split = '/')[[1]][-length(strsplit(plot_path, split = '/')[[1]])], collapse = '/'))))
  
  # extract needed lists from ML data
  indices = list_from_ML[['indices']]
  indices_test = list_from_ML[['indices_test']]
  predictions = list_from_ML[['predictions']]
  predictions_raw = list_from_ML[['predictions_raw']]
  
  # create one data frame to compare predicted vs actual status
  predictions_compare = data.frame()
  auc_values = c()   # store auc values herev
  for(i in seq_along(indices)){
    
    actuals = data_frame$status[indices_test[[i]]]
    probs   = predictions[[i]]$X1
    
    counts = table(actuals)
    
    # skip if only one class in test
    if(any(counts < min_per_condition) || length(unique(actuals)) < 2){
      warning(paste("Run", i, "skipped: too few IDs per condition or only one class"))
      next
    }
    
    roc_obj = pROC::roc(actuals, probs, direction = "<")
    auc_val = pROC::auc(roc_obj)
    auc_values = c(auc_values, auc_val)
    
    adding = data.frame(Class = actuals, 
                        predicted = probs, 
                        resample = paste0("run ", i), 
                        auc = as.numeric(auc_val))
    predictions_compare = rbind(predictions_compare, adding)
  }
  
  # calculate ROC for each run
  rocc = data.frame()
  for(run in unique(predictions_compare$resample)){
    onefold = dplyr::filter(predictions_compare, resample == run)
    
    if(length(unique(onefold$Class)) < 2){
      warning(paste("Run", run, "skipped: only one class in test set"))
      next
    }
    
    roc_run = pROC::roc(onefold$Class, onefold$predicted, direction = "<")
    rocc = rbind(rocc,
                 data.frame(Sp = roc_run$specificities,
                            Sn = roc_run$sensitivities,
                            n = seq_along(roc_run$sensitivities)))
  }
  
  # aggregate
  Sp = aggregate(Sp ~ n, rocc, mean)$Sp
  Sn = aggregate(Sn ~ n, rocc, mean)$Sn
  
  if (length(unique(predictions_compare$resample)) > 1) {
    errorSp = 1.96 * aggregate(Sp ~ n, rocc, sd)$Sp 
    errorSn = 1.96 * aggregate(Sn ~ n, rocc, sd)$Sn 
  } else {
    # no variability with a single bootstrap → set CI = 0
    errorSp = rep(0, length(Sp))
    errorSn = rep(0, length(Sn))
  }
  
  plotci = data.frame(Sp, Sn, errorSp, errorSn)
  plotci <- plotci[!(plotci$Sp == 1 & !(plotci$Sn %in% c(0,1))), ]
  
  # plot if path given
  if(plot_path != FALSE){
    pdf(plot_path,paper="a4r", width = 11, height = 8)
    print(
      ggplot(plotci, aes(x=(1-Sp), y=Sn)) + 
        geom_line(color = "darkblue") + 
        geom_ribbon(aes(x = 1 - Sp,
                        ymin = pmax(Sn - errorSn, 0),
                        ymax = pmin(Sn + errorSn, 1)),
                    fill = "#B2B2B2", alpha = 0.5) +
        theme_bw() +
        labs(x = "1 - Specificity",
             y = "Sensitivity",
             title = "Mean ROC curve with 95 % CI") +
        annotate("text", x = 0.4, y = 0.75, 
                 label = paste("mean AUC: ", round(mean(auc_values),2),
                               "\u00B1 ", round(sd(auc_values),2))) +
        geom_abline(slope = 1, intercept = 0, color="darkgrey", alpha = 0.3) +
        xlim(-0.1,1.1)
    )
    dev.off()
  }
  
  return(list(roc_data = plotci, auc_values = auc_values))
}

# old calculateROC
# calculateROC = function(list_from_ML, data_frame, plot_path = FALSE){
#   # check input
#   stopifnot('Data is not a dataframe.' = is.data.frame(data_frame),
#             'First argument is not a list (result of runML).' = inherits(list_from_ML, 'list') & all(names(list_from_ML) == c('models','importance','predictions','predictions_raw','indices')),
#             'Please give a valid directory to save the file' = ifelse(plot_path == FALSE, TRUE, dir.exists(paste(strsplit(plot_path, split = '/')[[1]][-length(strsplit(plot_path, split = '/')[[1]])], collapse = '/'))))
#   
#   # extract needed lists from ML data
#   indices = list_from_ML[['indices']]
#   predictions = list_from_ML[['predictions']]
#   predictions_raw = list_from_ML[['predictions_raw']]
#   
#   # create one data frame to compare predicted vs actual status
#   predictions_compare = data.frame()
#   for(i in 1:length(indices)){
#     adding = data.frame(Class = data_frame$status[-indices[[i]]], 
#                         predicted = predictions[[i]]$X1, 
#                         resample = paste0("run ", i), 
#                         auc = glmnet:::auc(data_frame$status[-indices[[i]]], predictions[[i]]$X1))
#     predictions_compare = rbind(predictions_compare,adding)
#   }
#   
#   # calculate ROC for each run
#   rocc = data.frame()
#   auc = c()
#   for(run in unique(predictions_compare$resample)){
#     onefold = dplyr::filter(predictions_compare, resample == run)
#     auc = c(auc,onefold$auc[1])
#     roc_run = roc(onefold$Class, onefold$predicted, direction = "<")
#     rocc = rbind(rocc, data.frame(Sp = roc_run$specificities, Sn = roc_run$sensitivities, n = rep(1:length(roc_run$sensitivities))))
#   }
#   
#   # aggregate the results and create new data frame
#   Sp = aggregate(Sp ~ n, rocc, mean)$Sp
#   Sn = aggregate(Sn ~ n, rocc, mean)$Sn
#   errorSp = 1.96 * aggregate(Sp ~ n, rocc, sd)$Sp 
#   errorSn = 1.96 * aggregate(Sn ~ n, rocc, sd)$Sn 
#   plotci = data.frame(Sp,Sn,errorSp,errorSn)
#   
#   # plot if path is given
#   if(plot_path != FALSE){
#     pdf(plot_path,paper="a4r", width = 11, height = 8)
#     print(ggplot(plotci, aes(x=(1-Sp),y=Sn)) + 
#             geom_line(aes(color = "darkblue")) + 
#             geom_ribbon(aes(x = 1 - Sp,
#                             ymin = pmax(Sn - errorSn, 0),
#                             ymax = pmin(Sn + errorSn, 1),
#                             fill = "#B2B2B2"), alpha = 0.5) +
#             theme_bw() +
#             labs(x = "1 - Specificity",
#                  y = "Sensitivity",
#                  title = "Mean ROC curve with 95 % CI") +
#             #   scale_y_continuous(expand = c(0,0), limits = c(0,1.02)) + scale_x_continuous(expand = c(0,0), limits = c(-0.01,1)) +
#             scale_color_manual(name = NULL, label = "Mean AUC", values = c("darkblue")) +
#             scale_fill_manual(name = NULL, label = "95 % CI", values = c('#B2B2B2') ) +
#             annotate("text", x = 0.4, y = 0.75, 
#                      label = paste("mean AUC: ", round(mean(auc),2), "\u00B1 ", round(sd(auc),2)) ) +
#             geom_abline(slope = 1, intercept = 0, color="darkgrey", alpha = 0.3)) +
#       xlim(-0.1,1.1)
#     dev.off()
#   }
#   
#   return(plotci)
#   
# }


###       Extract weights/importance and plot their average       ###
# plots top n (default = 50) proteins by averaged weight/importance
# returns plot (if path is given) and data frame

plotWeights = function(list_from_ML, plot_path = FALSE, number = 50){
  #check input
  stopifnot('First argument is not result from ML with algorithm linear regression, random forest or SVM (linear kernel)' = list_from_ML[['models']][[1]]$method %in% c('glmnet','ranger','svmLinear'),
            'Please give a valid directory to save the file' = ifelse(plot_path == FALSE, TRUE, dir.exists(paste(strsplit(plot_path, split = '/')[[1]][-length(strsplit(plot_path, split = '/')[[1]])], collapse = '/'))))
  
  # extract models to continue based on which algorithm input is from
  models = list_from_ML[['models']]
  
  ## linear regression ##
  if(models[[1]]$method == 'glmnet'){ 
    importance = list_from_ML[['importance']]
    
    # extract weights from all models
    coefficient = list()
    for(i in 1:length(models)){
      coefficient[[i]] = coef.glmnet(models[[i]]$finalModel, models[[i]]$bestTune$lambda)
    }
    
    #sort the weights
    weights_lm = vector("list",0)
    for (i in 1:length(coefficient)) {
      for (j in 1:length(coefficient[[i]])) {
        if(j == 1 || coefficient[[i]][j,1] == 0){ # skip intercept
          next
        }
        else if (row.names(coefficient[[i]])[j] %in% names(weights_lm)) { # protein already has a list -> append vector with value
          weights_lm[[row.names(coefficient[[i]])[j]]] = c(weights_lm[[row.names(coefficient[[i]])[j]]],coefficient[[i]][j,1])
        }
        else if (! row.names(coefficient[[i]])[j] %in% names(weights_lm)) {
          weights_lm[[row.names(coefficient[[i]])[j]]] = c(coefficient[[i]][j])
        }
        else{
          print("Failed to extract weights from linear regression model")
          break
        }
      }
    }
    
    # calculate statistics and put into data frame
    avg = c()
    error = c()
    picks = c()
    for (i in 1:length(weights_lm)) {
      avg = c(avg,mean(weights_lm[[i]]))
      error = c(error, sd(weights_lm[[i]]))
      picks = c(picks,length(weights_lm[[i]]))
    }
    
    weight = data.frame(avg,error, picks,row.names=names(weights_lm))
    weight = weight[which(abs(weight$avg)>abs(weight$error) & weight$picks > 1), ]
    
    weight = weight[order(abs(weight$avg), decreasing = TRUE), ]
    name = 'weight'
  }
  
  ## random forest ##
  else if(models[[1]]$method == 'ranger'){ 
    rf_imp = list_from_ML[['importance']]
    
    # create data frame with names for all proteins
    weight = data.frame(row.names = row.names(rf_imp[[1]]$importance))
    
    # fill data frame
    for(i in 1:length(rf_imp)){
      weight = cbind(weight, rf_imp[[i]]$importance, by = "row.names")
      weight$by = NULL
      names(weight)[names(weight) == "Overall"] = paste("run",i, sep = "")
    }
    
    # calculate statistics into data frame
    weight$avg = rowMeans(weight)
    weight$error = rowSds(as.matrix(weight[,-(ncol(weight))]))
    
    weight = weight[order(weight$avg, decreasing = TRUE), ]
    name = 'importance'
  }
  
  ## svm linear kernel ##
  else if(models[[1]]$method == 'svmLinear'){ 
    # calculate weights for each run
    for(i in 1:length(models)){
      coef = models[[i]]$finalModel@coef[[1]]
      matr = models[[i]]$finalModel@xmatrix[[1]]
      
      weig = as.data.frame(coef %*% matr)
      
      if(i == 1){
        weight = weig
      }
      else if(i>1 & all(names(weig) == names(weight))){
        weight = rbind(weight,weig)
      }
      else{
        print(paste(i," ERROR: failed to extract protein weights from SVM linear; differing protein names between runs "))
      }
    }
    
    # calculate statistics and add to data frame
    weight= as.data.frame(weight)
    weight$avg = rowMeans(weight)
    weight$error = rowSds(as.matrix(weight[,-ncol(weight)]))
    
    weight = weight[order(abs(weight$avg), decreasing = TRUE), ]
    name = 'weight'
  }
  
  weight <- weight %>%
    mutate(freq = picks / length(models))
  
  # plot if path is given
  if(plot_path != FALSE){
    pdf(plot_path,paper="a4r", width = 11, height = 8)
    print(ggplot(weight[1:number, ], 
                 #aes(x = reorder(rownames(weight[1:number, ]),avg), y = avg, fill = avg > 0)) +
                 aes(x = reorder(rownames(weight[1:number, ]),avg), y = avg, fill = freq)) +
      geom_bar(stat = "identity")+
      geom_errorbar( aes(x=reorder(rownames(weight[1:number, ]),avg), ymin=avg-error, ymax=avg+error, 
                         colour="black"), width=0.2, alpha=0.9, size=1.2) +
     labs(x = "Protein",
          y = "Coefficient (average ± SD)",
          title = paste("Bootstrap coefficients;" ,models[[1]]$method),
          fill = "Times selected") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_x_discrete(labels = abbreviate) + 
      #scale_fill_manual(name = NULL, label = NULL, values = c("darkorange", "darkgrey")) +
      scale_fill_gradient(low = "lightblue", high = "darkblue",labels = scales::percent) +
      scale_color_manual(name = NULL, label = "standard deviation", values = c("black")) +
      guides(colour = 'none'))
    dev.off()
  }
  return(weight)
}



###       Run Gene Set Enrichment analysis on given background and vector of genes       ###
# plots top n (default = 20) enriched KEGG pathways
# returns plot (if path is given) and data frame

gsea = function(genes, min_size = 3,plot_path = FALSE, number = 20, entrezID = FALSE, background){
  stopifnot('minSize has to be a number' = is.numeric(min_size),
            'background has to be a vector of type character' = is.vector(background, mode = 'character'),
            'Please give a valid directory to save the file' = ifelse(plot_path == FALSE, TRUE, dir.exists(paste(strsplit(plot_path, split = '/')[[1]][-length(strsplit(plot_path, split = '/')[[1]])], collapse = '/'))),
            'genes has to be a named vector with numeric value' = is.vector(genes, mode = 'numeric') & !any(is.na(names(genes))) | missing(background) )
  
  ## load pathways
  kegg_pathways = gmtPathways("data/c2.cp.kegg.v7.5.1.entrez.gmt")
  
  ## get entrezID of background and filter pathways
  # set mapping direction
  keyt = ifelse(entrezID == FALSE, 'SYMBOL', 'ENTREZID')
  
  if(!missing(background)){
    ids_b = AnnotationDbi::select(org.Hs.eg.db, 
                   keys = background,
                   columns = c("ENTREZID", "SYMBOL"),
                   keytype = keyt)
    background_id = ids_b$ENTREZID
    # let user now if certain proteins were excluded due to not mapped EntrezID
    if(any(is.na(background_id))){
      background_id = background_id[!is.na(background_id)]
      print(paste0('Could not find EntrezID for ', ids_b$SYMBOL[which(is.na(ids_b$ENTREZID))]),'. These will be ignored in further analysis.')
    }
    
    # filter KEGG pathways
    for(i in c(1:length(kegg_pathways))){
      kegg_pathways[[i]] = kegg_pathways[[i]][which(kegg_pathways[[i]] %in% background_id)]
    }
  }
  
  ## get entrezID of genes
  ids_g = AnnotationDbi::select(org.Hs.eg.db, 
                 keys = names(genes),
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = keyt)
  
  gene_id = ids_g$ENTREZID
  names(genes) = gene_id
  
  # let user now if certain proteins were excluded due to not mapped EntrezID
  if(any(is.na(ids_g$ENTREZID))){
    genes = genes[!is.na(names(genes))]
    print(paste0('Could not find EntrezID for', ids_g$SYMBOL[which(is.na(ids_g$ENTREZID))],'. These will be ignored in further analysis.' ))
  }
  
  ## run fgsea()
  # set scoreType depending on data 
  scoret = ifelse(all(genes >= 0), 'pos','std')
  
  fgsea_results = fgsea(pathways = kegg_pathways,
                        stats    = genes,
                        minSize  = min_size,
                        maxSize  = 500,
                        scoreType = scoret)
  
  # calculate number of genes in pathways 
  fgsea_results$count = sapply(fgsea_results$leadingEdge,length)
  # with small gene sample count is sometimes 0 -> filter
  fgsea_results = fgsea_results[fgsea_results$count != 0, ]
  
  # store gene names, not just IDs in column
  if(!missing(background)){
    fgsea_results$geneNames =   sapply(fgsea_results$leadingEdge, function(x) paste0(ids_g$SYMBOL[which(ids_g$ENTREZID %in% x)], collapse = ', '))
  }
  
  # get the actual length of the pathways as column
  length_pathways = as.data.frame(lapply(kegg_pathways,length))
  length_pathways = as.data.frame(t(length_pathways))
  names(length_pathways) = 'pathwaySize'
  length_pathways$pathway = row.names(length_pathways)
  
  # merge into one df 
  fgsea_results = merge(fgsea_results,length_pathways, by = "pathway")
  
  # plot if path is given
  if(plot_path != FALSE){
    pdf(plot_path,paper="a4r", width = 11, height = 8)
    print(ggplot(fgsea_results[order(fgsea_results$padj), ][1:number, ], aes(x = -log10(padj), y = reorder(pathway, -padj), fill = padj)) +
      geom_bar(stat = "identity") +
      xlab("- log(p-value)") +
      ylab("Functional category") +
      ggtitle("KEGG terms") +
      labs(fill = "FDR value") +
      scale_fill_viridis(option = "E")+ 
      geom_text(aes(label = paste0(count, "/", pathwaySize)), hjust=-0.15, size =3.5)+ 
      xlim(0, max(-log10(fgsea_results$padj))+0.02))
    dev.off()
  }
  
  return(fgsea_results)
}



### AUC curve of linear model from previously calculated weights
# calculates AUC of a linear model on a data set (protein_data) based on a given protein combination with average weights from previous linear model boot strapping
# returns data frame for plotting combinations to auc
auccurve = function(vectornames,weight_data,protein_data, maxn, add = F){
  stopifnot(
    'gene names in vector not represented in data frame' = vectornames %in% row.names(weight_data),
    'vectornames is not a character vector' = is.vector(vectornames, mode = 'character')
  )
  
  # if no maxn is given all possible combinations will be built
  if(missing(maxn)){
      maxn = length(vectornames)
    }
    
    # prepare combinations of defined length
    if(add == F){
    res = Map(combn, list(vectornames), seq_along(c(1:maxn)), simplyfy = F)
    test = unlist(res, recursive = FALSE)
    vectors = list()
    z=1
    for(i in 1:(length(res))){
      for(j in 1:(length(res[[i]])/length(res[[i]][ ,1]) ) ){
        vectors[[z]] = res[[i]][,j]
        z = z + 1
      }
    }
    }
  # if add is set to TRUE, combinations are only adding one protein after the other in order they are given (going from legnth 1 to length of vectornames)
  else{
    vectors = list()
    z=1
    for(i in 1:(length(vectornames))){
      
      vectors[[i]] = vectornames[1:i]
      z = z + 1
    }
  }
  
  ## put combinations into empty dataframe
  savespace = data.frame(combinations = I(vectors))
  savespace$auc = NA
  
  ## add auc for the combinations
  for(i in 1:nrow(savespace)){
    interim = as.data.frame(weight_data)
    # set all weights that are not in the combination of interes to 0
    interim[! row.names(interim) %in% vectors[[i]], "V1"] = 0
    # predict 0 or 1 by using the linear regression formula (y = w1*x1+ .... w_n*x_n)
    interim_pred = as.matrix(as.matrix(protein_data[ ,!names(protein_data)=='status']) %*% as.matrix(interim))
    # use glmnet function to calculate the auc between actual status and the calculated status, i.e. prediction
    savespace$auc[i] = glmnet:::auc(interim_pred,protein_data$status)
  }

  savespace$combinations = lapply(savespace$combinations, function(X) paste0(X,collapse = ", "))
  return(savespace)
}