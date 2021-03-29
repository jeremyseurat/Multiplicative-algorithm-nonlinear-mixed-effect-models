multiplicative_algo<-function(FIMs, version, nbParameters, designSpace, nbSubjects, lambda, delta, iteration){
  
  numberFIMs <- length(FIMs)
  w<- rep(1/numberFIMs,numberFIMs)
  
  if( version == "R")
  {
    for(k in 1:iteration){
      #Calculate the sum of matrix with weight associated
      #idea:   Mw<-Mw+w[i]*FIMs[[i]]
      MFw<-Map("*",FIMs,w)     #associate weights
      Mw<-Reduce("+",MFw)     #sum of all matrices
      dm<-dim(Mw)[1]         
      Dphi <- det(Mw)^(1/dm) * solve(Mw)/dm  # calculate derivatives of function phi_D
      d<-sapply(FIMs, FUN = function(x) sum(diag(Dphi %*% x)))         #the vector of multiplier
      w<- w*d^lambda/sum(w*d^lambda)          #develop weight
      if(max(d)<(1+delta)*sum(w*d)){break}  #stop criterion
    }
    Dcriterion <- det(Mw*nbSubjects)^(1/dm)
  }
  if( version == "C")
  {
    dyn.load("algomult.dll")
    q<- .C("multiplicativeAlgorithm", originMatrices=as.double(unlist(FIMs)), weights = as.double(w), dim=as.integer(dim(FIMs[[1]])[1]), 
           n=as.integer(numberFIMs), dimA = as.integer(nbParameters), lambda = as.double(lambda), delta = as.double(delta), 
           it = as.integer(iteration))
    dyn.unload("algomult.dll")
    w <- q$weights
    k <- q$it
    MFw<-Map("*",FIMs,w)     
    Mw<-Reduce("+",MFw)   
    Dcriterion <- det(Mw*nbSubjects)^(1/q$dim)
  }
 
  #by default, we take into count weights which give meaningful subject numlber
  v<-which(w>1/nbSubjects)   # v<-which(w>mean(w))  ##if we want to get weights bigger than the mean
  
  cat("\n\n*******************************************************\n")
  cat("*** Result of Algorithm ****************\n")
  cat("\n\n IMPORTANT SUBJECT NUMBER ******************\n", w[v]*nbSubjects)
  cat("\n\n IMPORTANT DESIGN VALUES ******************\n")
  print(designSpace[,v])
      
  cat("\n\n NUMBER OF INTERATIONS ******************\n",k,
      "\n\n D-CRITERION ****************************\n",Dcriterion)
  
 
  plot(w,ylim=c(0,1))
  return(list(w,v))
}

