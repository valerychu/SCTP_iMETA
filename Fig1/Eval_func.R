library(DescTools)
library(MLmetrics)
library(ROCR)
library(pROC)
library("gridExtra")
#setwd("/Users/w435u/Documents/DrugAI")



func_Eval <- function(pred_vec, resp_vec){
  #pred_vec[is.na(pred_vec)] <- 0
 
  pred_vec <- pred_vec[which(!is.na(pred_vec))]
  if(length(pred_vec)==0){
    eval_all <- rep(NA, 6)
  } else {
    
    resp_vec <- resp_vec[which(!is.na(pred_vec))]
    
    pred_vec <-  sapply(pred_vec, function(x) min(max(x, 1E-15), 1-1E-15)) 
    bs <- BrierScore(resp_vec, pred_vec, scaled = FALSE)
    possibleError <- tryCatch(pred.mod <- prediction(pred_vec, resp_vec), error = function(e) e)
    
    if(inherits(possibleError, "error")){
      auc = aucpr = mxe = NA
    } else {
      pred.mod <- prediction(pred_vec, resp_vec)
      auc.tmp <- performance(pred.mod,"auc")
      auc <- as.numeric(auc.tmp@y.values)
      aucpr.tmp <- performance(pred.mod,"aucpr")
      aucpr <- as.numeric(aucpr.tmp@y.values)
      mxe.tmp <- performance(pred.mod, "mxe")
      mxe <- as.numeric(mxe.tmp@y.values)
      phi.tmp <- performance(pred.mod, "phi")
      phi <- as.numeric(unlist(phi.tmp@y.values)[2])
      #recall <- sensitivity(pred_vec, resp_vec, positive="1")
    }
    
    pred_class <- round(pred_vec)
    f1 <- F1_Score(resp_vec, pred_class)
    acc <- 1-mean(abs(resp_vec - pred_class))
    
    if(sum(pred_class)==0){
      sensitivity = 0
    } else{
      conf_mat <- table(Predicted = pred_class, True = resp_vec)
      TP <- conf_mat["1", "1"]
      FN <- conf_mat["0", "1"]
      sensitivity <- TP / (TP + FN)
    }
    
    #logLoss2 <- -mean(resp_vec * log(pred_vec) + (1 - resp_vec) * log(1 - pred_vec))
    eval_all <- c(auc, aucpr, 1-bs, f1, acc, mxe, sensitivity, phi)
  }
  
  
  names(eval_all) <- c("AUC", "AUCPR", "Brier", "F1","ACC", "LogLoss", "Sens", "Matthews")

  return(eval_all)
}


  