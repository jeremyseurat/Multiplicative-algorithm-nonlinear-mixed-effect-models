funFIMpop <- function(equation,parameters,beta,o,sigma,t_group,Trand,d,PropSubjects,nbTot){
  #List of names of the fixed effects parameters and sampling times
  ParametersInEquation<-c(parameters,"t")
  
  PSI<-c(parameters,"sig.inter","sig.slope")
  lpsi<-length(PSI)
  lengthParameters<-length(parameters)
  
  #(Diagonal Matrix of) variance for inter-subject random effects:
  omega<-diag(o)
  
  #Parameters of standard deviation model of residual error (sig.inter+sig.slope*f)^2:
  sig.inter<-sigma[1]
  sig.slope<-sigma[2]
  
  #get derivatives of sigma
  Vmodel<- parse(text = paste0("( sig.inter + sig.slope * fixed )^2"))
  dv<-deriv(Vmodel[[1]],PSI)
  
  #gather all groups of protocol
  M_f<-list()
  M_F <- matrix(rep(0),nrow=lpsi+lengthParameters,ncol=lpsi+lengthParameters)
  for(q in 1:length(t_group)){
    #design of sampling times
    t<-t_group[[q]]
    
    #dose value
    dose<-d[q]
    
    #calculate matrix E for n individuals
    equatf <- parse(text = equation, n=-1)
    f<-function(ParametersInEquation){eval(equatf[[1]])}
    
    #Fixed effects parameters values
    for(i in 1:lengthParameters){
      assign(parameters[i],beta[i])
    }
    
    #calculate the observations with personnal parameters
    parameterValues <- c(beta,t)
    fixed<-f(parameterValues)
    
    #calculate variance Var
    var<-diag((sig.inter+sig.slope*fixed)^2)
    
    #get derivatives for fixed parameters
    df<-deriv(equatf[[1]],PSI)
    mdf<-attributes(eval(df))$gradient
    #delete the last two columns (correspond to sig.inter and sig.slope) 
    mdfi <- mdf[,-c(lpsi-1,lpsi)]
    #complete derivative for random effect model   
    #Random effect model (1) = additive  (2) = exponential 
    #------------------------------------------------------------------
    beta0 <- beta
    beta0[which(Trand==1)] <- 1 
    mdfie <- mdfi %*% diag(beta0)
    
    #calculate variance Vi
    Vi <- mdfie %*% omega %*% t(mdfie) + var
    #inverse of matrix Vi
    SVi <- solve(Vi)
    #compute derivatives of sigma
    mdv<-attributes(eval(dv))$gradient
    
    
    #calculate matrix part A
    M_A <- t(mdfi) %*% SVi %*% mdfi
    
    #complete the rest of the matrix with 0
    M_A <- cbind(M_A,matrix(rep(0),ncol=lpsi,nrow=lengthParameters))
    M_A <- rbind(M_A,matrix(rep(0),ncol=lpsi+lengthParameters,nrow=lpsi))
    
    #calculate matrix part B
    #initialize the matrix with 0
    M_B <- matrix(rep(0),nrow=lpsi+lengthParameters,ncol=lpsi+lengthParameters)
    #prepare a list of matrix of derivatives of sigma to simplify usage
    if(length(t)==1){
      m <- lapply(c((lpsi-1),lpsi),function(i,mdv) mdv[i],mdv=mdv )
    }else{
      m <- lapply(c((lpsi-1),lpsi),function(i,mdv) diag(mdv[,i]),mdv=mdv )
    }
    
    #calculate first three rows of part B
    for(i in 1:lengthParameters){
      mdfiei <- (mdfie[,i] %*% t(mdfie[,i]))
      
      M_B[lengthParameters+i,seq(lengthParameters+1,lengthParameters+lengthParameters)] <- sapply( 
        lapply(seq(1,lengthParameters), function(i,mdfie) mdfie[,i],mdfie=mdfie), 
        function(x) 1/2 * sum(diag((mdfiei) %*% SVi %*% (x %*% t(x)) %*% SVi)))
      
      M_B[lengthParameters+i,c(lengthParameters+lpsi-1,lengthParameters+lpsi)] <- sapply(
        m,function(x) 1/2 * sum(diag((mdfiei) %*% SVi %*% x %*% SVi)))
    }
    #calculate the last two rows of partB
    for(i in (lpsi-1):lpsi){
      
      M_B[lengthParameters+i,seq(lengthParameters+1,lengthParameters+lengthParameters)] <- sapply(
        lapply(seq(1,lengthParameters), function(i,mdfie) mdfie[,i],mdfie=mdfie),
        function(x) 1/2 * sum(diag((m[[i-lengthParameters]]) %*% SVi %*% (x %*% t(x)) %*% SVi)))
      
      M_B[lengthParameters+i,c(lengthParameters+lpsi-1,lengthParameters+lpsi)] <- sapply(
        m,function(x) 1/2 * sum(diag(((m[[i-lengthParameters]])) %*% SVi %*% x %*% SVi)))
    }
    
    M_f[[q]] <- (M_A+M_B)*PropSubjects[q]
    M_F <-M_F+M_f[[q]]
  }
  M_F <- M_F *nbTot
  
  #set names for vectors 
  fname <-c(sapply(1:lengthParameters, function(x) paste0("u_",parameters[x])),sapply(
    1:lengthParameters, function(x) paste0("w2_",parameters[x])),"sig.inter","sig.slope")
  rownames(M_F) <- fname
  colnames(M_F) <- fname
  
  if(0 %in% c(o,sigma)){
    M_F<-M_F[-c(lengthParameters+which(c(o,sigma)==0)),-c(lengthParameters+which(c(o,sigma)==0))]
    PSI<- PSI[-c(lengthParameters+which(c(o,sigma)==0))]
  }
  
  if(length(t)==1 ){
    return(list(M_F,det(M_F)))
  }else{
    tryCatch({
      deterFim <- det(M_F)
      SE <- sqrt(diag(solve(M_F)))
      RSE <- 100 * SE / c(beta,o,sigma)[which(c(beta,o,sigma)!=0)]
      CritereDopt <- deterFim^(1/length(c(beta,o,sigma)[which(c(beta,o,sigma)!=0)]))
      #write the output in console
      #       cat("******************* FISHER INFORMATION MATRIX *************************\n")
      #       print(M_F)
      #       
      #       cat("\n\n******************* DETERMINANT OF THE MATRIX *********************\n", deterFim,
      #           "\n\n******************* CRITERION *************************************\n",CritereDopt,
      #           "\n\n******************* STANDARD ERROR ********************************\n",SE,
      #           "\n\n******************* RELATIVE STANDARD ERROR ***********************\n",RSE,"\n\n")
      return(list(M_F,deterFim,CritereDopt,SE,RSE))
    },error=function(e){
      return(list(M_F,det(M_F)))
    })
    
    #     draw the graphical illustration of model
    #     t<-seq(min(unlist(t_group)),max(unlist(t_group)),(max(unlist(t_group))-min(unlist(t_group)))/100)
    #     x<-c(beta,t)
    #     fvalue<-f(x)
    #     lines(t,fvalue) #then plot the line of function
  }
   
}