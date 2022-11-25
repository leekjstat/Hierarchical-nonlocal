
## hierarchical pMoM prior with the reduced SSS (RSSS) algorithm

# new argument: method, which is one of "marginal", "working", "dev" or "pear"
# It determines the criterion for the screening
# - marginal          : marginal correlations between X's and y's
# - working, dev, pear: correlations between X's and (working, dev, pear) residuals

H_RSSS = function(X, y, lambda1, lambda2, ind_fun, N=1000, C0=1, verbose=FALSE, r=1, init.mod=NULL, method="marginal"){
  n = nrow(X)
  p = ncol(X)
  
  m.n = round(n/2)
  
  # for screened \Gamma+
  K1 = 10 
  K2 = 10
  if(method=="marginal"){
    abs.corr = abs(c(t(X) %*% y))
  }
  
  requireNamespace("Matrix")
  requireNamespace("abind")
  abind = abind::abind
  
  verb = verbose # set verbose
  require(Matrix)
  library(abind)
  
  ### initialize the model
  if(is.null(init.mod)){
    gam  = rep(0,p)
    p.g  = 0
    ind2 = which(gam==1)
    curr = -1000000
  }else{
    sam      = init.mod
    gam      = rep(0,p)
    gam[sam] = 1
    ind2     = which(gam==1)
    p.g      = sum(gam)
    curr     = ind_fun(X[,ind2, drop=FALSE], y, r, lambda1, lambda2) # evaluate the objective value on the initial model
  }
  
  GAM.fin0 = NULL # to save the binary vector of the searched models 
  OBJ.fin0 = NULL # to save the objective values of the searched models 
  
  GAM = gam 
  OBJ = curr
  obj = OBJ
  
  ### search the neighborhood of the initial model
  C.p       = rep(-1000000,p)
  if(method=="marginal"){
    cut.off   = abs.corr[order(abs.corr, decreasing=TRUE)][K1]
    top.ind   = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
  }else{
    if(length(ind2) > 0){
      glm.fit   = glm(y ~ X[,ind2, drop=FALSE] - 1, family = "binomial")
      resid.vec = resid(glm.fit, type=method)
      abs.corr = abs(c(t(X) %*% resid.vec))
      cut.off = abs.corr[order(abs.corr, decreasing=TRUE)][K1]
      top.ind = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
    }else{ # if length(ind2) == 0
      top.ind = sample(x=(1:p), size=K1)
    }
  }
  inter.ind = intersect(top.ind, ind2)
  if(length(inter.ind) > 0){
    K11     = length(inter.ind)
    cut.off = abs.corr[order(abs.corr, decreasing=TRUE)][K1+K11]
    top.ind = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
    top.ind = setdiff(top.ind, ind2)
  }
  rest.ind   = sample(x=(1:p)[-c(top.ind, ind2)], size=K2)
  Gamma.plus = union(top.ind, rest.ind)
  for(j in Gamma.plus){
    gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
    int   = ind_fun(X[,ind.p, drop=FALSE], y, r, lambda1, lambda2)
    obj.p = c(int)
    if(is.na(obj.p)==TRUE){obj.p = -1000000}
    C.p[j] = obj.p
  }
  
  C.m       = rep(-1000000,p)
  if(p.g > 1){
    IND.m   = ind2
    p.ind.m = length(IND.m)
    for(i in 1:p.g){
      j     = ind2[i]
      gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)  
      int   = ind_fun(X[,ind.m, drop=FALSE], y, r, lambda1, lambda2)
      obj.m = c(int)
      if(is.na(obj.m)==TRUE){obj.m = -1000000}
      C.m[j] = obj.m  
    }
  }
  
  C.s = matrix(-1000000,round(n/2),p)
  if(p.g > 1){
    for(i in 1:p.g){
      for(u in Gamma.plus){
        j        = ind2[i]
        gam.s    = gam
        gam.s[j] = 0
        gam.s[u] = 1
        ind.s    = which(gam.s==1)  
        int      = ind_fun(X[,ind.s, drop=FALSE], y, r, lambda1, lambda2)
        obj.s    = c(int)
        if(is.na(obj.s)==TRUE){obj.s = -1000000}
        C.s[i,u] = obj.s
      }
    }
  }
  
  
  # print("#################################")
  # print("RSSS starts running")
  
  ### create the arrays to save the neighborhood of visited models 
  d0     = 200
  OBJ.s0 = array(-1000000,dim=c(round(n/2),p,d0));OBJ.s0[,,1] = C.s
  OBJ.m0 = matrix(-1000000,p,d0);OBJ.m0[,1] = C.m
  OBJ.p0 = matrix(-1000000,p,d0);OBJ.p0[,1] = C.p
  
  ID = sum(2^(3*log(ind2))) # set the id of the initial model
  
  ############################################# Start !!!
  for(uu in 1:C0){ # repeat the RSSS algorithm C0 times  
    
    GAM.total     = Matrix(0,p,1000000,sparse=TRUE)
    OBJ.total     = rep(-1000000,1000000)
    GAM.total[,1] = gam
    OBJ.total[1]  = obj
    time.total    = rep(0,1000000)
    
    pmt0 = proc.time()
    for(iter in 1:N){ 
      
      id = sum(2^(3*log(ind2))) # calculate the id of the current model
      id.ind = which(id==ID) # check whether the current model has already been visited using the id  
      leng = length(id.ind) # how many times visited?
      
      if(leng==0){ 
        ### if the current model has not been visited, search the neighborhood of the current model
        ID = c(ID,id)
        
        #########################################################
        ## (1) Search the set \Gamma+ (screened)
        ## Consider only the (top K1 + randomly selected K2) variables 
        ## based on absolute correlations between X's and residuals
        #########################################################
        C.p = rep(-1000000,p)
        p.g = sum(gam)
        if(p.g <= m.n){ # if the size of the current model(=p.g) is larger than m.n, we won't explore the \Gamma+ set
          if(method=="marginal"){
            cut.off   = abs.corr[order(abs.corr, decreasing=TRUE)][K1]
            top.ind   = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
          }else{
            if(p.g > 0){
              glm.fit   = glm(y ~ X[,ind2, drop=FALSE] - 1, family = "binomial")
              resid.vec = resid(glm.fit, type=method)
              abs.corr = abs(c(t(X) %*% resid.vec))
              cut.off = abs.corr[order(abs.corr, decreasing=TRUE)][K1]
              top.ind = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
            }else{ # if length(ind2) = p.g == 0
              top.ind = sample(x=(1:p), size=K1)
            }
          }
          inter.ind = intersect(top.ind, ind2)
          if(length(inter.ind) > 0){
            K11     = length(inter.ind)
            cut.off = abs.corr[order(abs.corr, decreasing=TRUE)][K1+K11]
            top.ind = which(abs.corr >= cut.off) # the set \Gamma+ (screened)
            top.ind = setdiff(top.ind, ind2)
          }
          rest.ind   = sample(x=(1:p)[-c(top.ind, ind2)], size=K2)
          Gamma.plus = union(top.ind, rest.ind)
          for(j in Gamma.plus){
            gam.p    = gam
            gam.p[j] = 1
            ind.p    = which(gam.p==1)
            int      = ind_fun(X[,ind.p, drop=FALSE], y, r, lambda1, lambda2)
            obj.p    = c(int)
            if(is.na(obj.p)==TRUE){obj.p = -1000000}
            C.p[j]    = obj.p
            ind.total = which(OBJ.total< -100000)[1]
            OBJ.total[ind.total]  = obj.p
            GAM.total[,ind.total] = gam.p
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
        }
        
        
        #########################################################
        ## (2) Search the set \Gamma-
        #########################################################
        C.m     = rep(-1000000,p)
        if(p.g > 1){ # if the size of the current model(=p.g) is smaller than 2, we won't explore the \Gamma- set
          IND.m   = ind2
          p.ind.m = length(IND.m)
          for(i in 1:p.g){
            j        = ind2[i]
            gam.m    = gam
            gam.m[j] = 0
            ind.m    = which(gam.m==1)  
            int      = ind_fun(X[,ind.m, drop=FALSE], y, r, lambda1, lambda2)
            obj.m    = c(int)
            if(is.na(obj.m)==TRUE){obj.m = -1000000}
            C.m[j]    = obj.m  
            ind.total = which(OBJ.total< -100000)[1]
            OBJ.total[ind.total]  = obj.m
            GAM.total[,ind.total] =  gam.m
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
        }
        
        
        #########################################################
        ## (3) Search the set \Gamma0 (screened)
        ## When adding new variable, consider only the (top K1 + randomly selected K2) variables 
        ## based on absolute correlations between X's and residuals
        #########################################################
        C.s = matrix(-1000000,round(n/2),p)
        for(i in 1:p.g){
          for(u in Gamma.plus){
            j        = ind2[i]
            gam.s    = gam
            gam.s[j] = 0
            gam.s[u] = 1
            ind.s    = which(gam.s==1)  
            int      = ind_fun(X[,ind.s, drop=FALSE], y, r, lambda1, lambda2)
            obj.s    = c(int)
            if(is.na(obj.s)==TRUE){obj.s = -1000000}
            C.s[i,u]  = obj.s 
            ind.total = which(OBJ.total< -100000)[1]
            OBJ.total[ind.total]  = obj.s
            GAM.total[,ind.total] =  gam.s
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
        }
        
        
        
        #########################################################
        d = which(apply(OBJ.p0,2,mean) < (min(OBJ.p0[,1])+1))[1]
        if(is.na(d)==TRUE){
          OBJ.s0 = abind(OBJ.s0,C.s)
          OBJ.p0 = cbind(OBJ.p0,C.p)
          OBJ.m0 = cbind(OBJ.m0,C.m)  
        }else{
          OBJ.s0[,,d] = C.s
          OBJ.p0[,d]  = C.p
          OBJ.m0[,d]  = C.m
        }
        ### if the current model has not been visited, search the neighborhood of the current model (end)
      }else{
        ### if the current model has already been visited, call the saved neighborhood 
        C.p = OBJ.p0[,(id.ind[1])];C.m = OBJ.m0[,(id.ind[1])];C.s = OBJ.s0[,,(id.ind[1])]
      }
      
      
      ### based on the neighborhood, choose a next model to move by randomly sampling proportion to the posterior probability 
      prop = exp(C.p-max(C.p))
      prop[which(is.na(prop))] = 0
      if(sum(prop)<0.1){prop[1]=1}
      sample.p = sample(1:length(prop), 1, prob=prop)
      obj.p = C.p[sample.p]
      
      prop = exp(C.m-max(C.m))
      prop[which(is.na(prop))] = 0
      if(sum(prop)<0.1){prop[1]=1}
      sample.m = sample(1:length(prop), 1, prob=prop)
      obj.m = C.m[sample.m]
      
      prop = exp(C.s-max(C.s))
      prop[which(is.na(prop))] = 0
      if(sum(prop)<0.1){prop[1]=1}
      sample.s = sample(1:length(prop), 1, prob=prop)
      obj.s = C.s[sample.s]
      
      l1 = 1/(1+exp(obj.m-obj.p)+exp(obj.s-obj.p))
      l2 = 1/(1+exp(obj.p-obj.m)+exp(obj.s-obj.m))
      l3 = 1-l1-l2
      if(l3<0){l3=0}
      z = sample(1:3, 1, prob=c(l1,l2,l3))
      
      if(z==1){gam[sample.p]=1;obj = obj.p;curr=obj.p}
      if(z==2){gam[sample.m]=0;obj = obj.m;curr=obj.m}
      if(z==3){ 
        wh = which(obj.s==C.s, arr.ind=TRUE)
        gam[ind2[wh[1]]]=0;gam[wh[2]] = 1
        obj = obj.s;curr=obj.s
      }
      ###
      
      ind2 = which(gam==1) # generate the binary variable for the new model
      p.g  = sum(gam) # the size of the new model
      
      
      gam.pr  = GAM.total[,which.max(OBJ.total)] # the MAP model among the searched models 
      obj.pr  = max(OBJ.total) # its objective value
      ind2.pr = which(gam.pr==1) # its size
      if(verb==TRUE&iter%%50==0){
        print("#################################")
        cat("# of iterations: ",iter, "\n")
        cat("The Searched MAP Model, ",ind2.pr,", with the Evaluated Object Value ",obj.pr,"\n")
        cat("The Current Model,      ",ind2,", with the Evaluated Object Value ",curr,"\n")    
        ind_fun(X[,ind2, drop=FALSE], y, r, lambda1, lambda2, display = TRUE) # print the objective value of the current model
      }
      
      # cat(iter, "th iteration is completed. \n")
    }
    time0 = proc.time()-pmt0
    # print(time0)
    
    cat("RSSS-H-pMOM: The Selected Variables in the Searched MAP Model: ", ind2.pr, "\n")
    cat("The Number of Total Searched Models",length(unique(OBJ.total)),"\n") 
    
    ### after a single SSS algorithm runs, summarize the searched models and their posterior model probability
    ind.total = which(OBJ.total> -100000) 
    OBJ.fin = unique(OBJ.total[ind.total]) # generate the unique objective values from models that have been visited multiple times
    
    w = length(OBJ.fin)
    time.fin = rep(0,w)
    GAM.fin = matrix(0,p,w)
    for(i in 1:length(OBJ.fin)){
      GAM.fin[,i] = GAM.total[,which(OBJ.total==OBJ.fin[i])[1]] # generate the unique binary variables of models from models that have been visited multiple times
    }
    
    const = sum(exp(OBJ.fin-max(OBJ.fin))) # calculate the normalizing constant
    posterior = exp(OBJ.fin-max(OBJ.fin))/const # calculate the posterior model probability
    total.size = length(OBJ.fin) # the total number of models that are searched by the SSS
    m = max(OBJ.fin) # the objective value of the MAP model
    ind.m0 = which.max(OBJ.fin) 
    gam = GAM.fin[,ind.m0] # the binary variable of the MAP model
    ind2 = which(gam==1);p.g = sum(gam) # its index and size
    GAM.fin0 = cbind(GAM.fin0,GAM.fin) # save the binary variables
    OBJ.fin0 = c(OBJ.fin0,OBJ.fin) # save the the objective values corresponding to the binary variables
  }
  
  
  ### after C0 number of single SSS algorithm runs, summarize the searched models and their posterior model probability
  if(C0 == 1){
    GAM = GAM.fin0 # the binary variables of searched models
    OBJ = OBJ.fin0 # the the objective values corresponding to the binary variables
  }
  if(C0 > 1){
    # print("#################################")
    # print("Post-process starts")
    # print("#################################")
    
    OBJ.fin1 = unique(OBJ.fin0)
    w = length(OBJ.fin1)
    time.fin = rep(0,w)
    GAM.fin1 = Matrix(0,p,w,sparse=TRUE);GAM.fin1[,1] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[1])[1]]
    for(i in 2:w){
      GAM.fin1[,i] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[i])[1]]
      #  time.fin[i] = time.total[which(OBJ.total==OBJ.fin[i])[1]]
    }
    rm(GAM.fin0)
    GAM = GAM.fin1 # the binary variables of searched models
    OBJ = OBJ.fin1 # the the objective values corresponding to the binary variables
    # print("Done!")
  }
  
  
  return(list(GAM = GAM,OBJ = OBJ, time0 = time0))
}

