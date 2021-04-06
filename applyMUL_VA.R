model <- "exp(VA)+(1-exp(-exp(k)*t))*(emax*dose/(exp(ED50)+dose)-exp(be)*exp(VA))"
paramName <- c("VA","k","be","emax","ED50")
mu <- c(log(55),log(0.005),log(0.2),30,log(150))
omega <- c(0.07,0.5,1,150,0)
sigma <- c(sqrt(28),0)
allTimes <- c(0,7,seq(28,672,28))
equiSpacedTimes <- c(0,224,420,672)
fixed_times <- list(c(1,0),c(4,672))
Trand <- c(1,1,1,1,1)
dose <- c(0,25,50,100,150,300,500)
lambda <- 1
delta <- 0.0001
iteration <- 10000
nbSubjects <- 300

# ############################################################################################################
# # Reference 26 times
# start <- Sys.time()
# ResultEval <- funFIMpop(model, paramName, mu, omega, sigma, list(allTimes, allTimes, allTimes, allTimes), 
#                         Trand, c(0,150,300,500),c(75,75,75,75)/nbSubjects,nbSubjects)
# end <- Sys.time()
# end-start   
# 
# # Fix 26 times and optimize doses
# # Create design space and elementary FIM corresponding
# start <- Sys.time()
# nbTimes <- 26
# Designs_ini<-combn(allTimes,nbTimes)
# FIMs<-list()
# Designs<-matrix(nrow=nbTimes+1)
# n<-1
# for(i in 1:dim(Designs_ini)[2]){
#   for(j in 1:length(dose)){
#     #the number of design groups need to change in accordance with the number of doses
#     FIMs[[n]] <-funFIMpop(model,paramName,mu,omega,sigma,list(Designs_ini[,i]),Trand,dose[j],1,1)[[1]]
#     Designs<-cbind(Designs,c(dose[j], Designs_ini[,i]))
#     n<-n+1
#   }
# }
# Designs <- Designs[,-1]
# end <- Sys.time()
# end-start   
# 
# # Optimize using multiplicative algorithm implemented in R or C
# start <- Sys.time()
# weights <- multiplicative_algo(FIMs, "C", length(paramName), Designs, nbSubjects, lambda, delta, iteration)
# end <- Sys.time()
# end-start   

############################################################################################################
# Reference 4 equi-spaced times
start <- Sys.time()
ResultEval <- funFIMpop(model, paramName, mu, omega, sigma, list(equiSpacedTimes, equiSpacedTimes, 
                                                                 equiSpacedTimes, equiSpacedTimes),
                        Trand, c(0,150,300,500),c(75,75,75,75)/nbSubjects,nbSubjects)
end <- Sys.time()
end-start 

# # Fix 4 equi-spaced times and optimize doses
# # Create design space and elementary FIM corresponding
# start <- Sys.time()
# nbTimes <- 4
# Designs_ini<-combn(equiSpacedTimes,nbTimes)
# FIMs<-list()
# Designs<-matrix(nrow=nbTimes+1)
# n<-1
# for(i in 1:dim(Designs_ini)[2]){
#   for(j in 1:length(dose)){
#     #the number of design groups need to change in accordance with the number of doses
#     FIMs[[n]] <-funFIMpop(model,paramName,mu,omega,sigma,list(Designs_ini[,i]),Trand,dose[j],1,1)[[1]]
#     Designs<-cbind(Designs,c(dose[j], Designs_ini[,i]))
#     n<-n+1
#   }
# }
# Designs <- Designs[,-1]
# end <- Sys.time()
# end-start  
# 
# # Optimize using multiplicative algorithm implemented in R or C
# start <- Sys.time()
# weights <- multiplicative_algo(FIMs, "C", length(paramName), Designs, nbSubjects, lambda, delta, iteration)
# end <- Sys.time()
# end-start   

############################################################################################################
# Optimize both doses and 4 times
# Create design space and elementary FIM corresponding
start <- Sys.time()
nbTimes <- 4
Designs_ini<-combn(allTimes,nbTimes)
if(length(fixed_times)>0){
  for(f in 1:length(fixed_times)){
    Designs_ini<-Designs_ini[,which(Designs_ini[fixed_times[[f]][1],]==fixed_times[[f]][2])]
  }
}

FIMs<-list()
Designs<-matrix(nrow=nbTimes+1)
o<-1
for(l in 1:dim(Designs_ini)[2]){
  for(m in 1:length(dose)){
    #the number of design groups need to change in accordance with the number of doses
    FIMs[[o]] <-funFIMpop(model,paramName,mu,omega,sigma,list(Designs_ini[,l]),Trand,dose[m],1,1)[[1]]
    Designs<-cbind(Designs,c(dose[m], Designs_ini[,l]))
    o<-o+1
  }
}
Designs <- Designs[,-1]
end <- Sys.time()
end-start   

# Optimize using multiplicative algorithm implemented in R or C
start <- Sys.time()
weights <- multiplicative_algo(FIMs, "C", length(paramName), Designs, nbSubjects, lambda, delta, iteration)
end <- Sys.time()
end-start 