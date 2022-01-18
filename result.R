
Evaluation <- function(beta1, beta2){
  true.index <- which(beta1==1)
  false.index <- which(beta1==0)
  positive.index <- which(beta2==1)
  negative.index <- which(beta2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  res = matrix(0, nrow=1, ncol=9)
  colnames(res) <- c("Precision", "Recall", "Sensitivity", "Specificity", "MCC", "TP", "TN", "FP", "FN")
  res[1,] = c(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,MCC=MCC,TP=TP,TN=TN,FP=FP,FN=FN)
  
  return(res)
}


result <- function(fit){
#GAM = fit[-1,]; OBJ = fit[1,]

  GAM = fit$GAM; OBJ = fit$OBJ; tuning= fit$tuning
  p = nrow(GAM)
  marg.gam = rep(0,p)
  for(u in 1:ncol(GAM)){
    marg.gam = marg.gam + GAM[,u]*exp(OBJ[u]-max(OBJ))
  }
  marg.gam = marg.gam / sum(exp(OBJ-max(OBJ)))
  gam0 = GAM[,which.max(OBJ)]
  ind2 = which(gam0==1)
  post = exp(OBJ-max(OBJ))/sum(exp(OBJ-max(OBJ)))
  hppm = 1/sum(exp(OBJ-max(OBJ)))
  
  print("# of Searched Models by S5");print(length(OBJ))
  print("The MAP model is ")
  print(which(gam0==1))
  print(paste("with posterior probability",round(hppm,3) )) 
  
  return(list(hppm = which(gam0==1), hppm.prob = hppm, marg.prob = marg.gam,gam = GAM, obj = OBJ, post = post, tuning = tuning) )
}
