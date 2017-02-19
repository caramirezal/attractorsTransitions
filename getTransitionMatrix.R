# getTransitionsMatrix() simulate boolean network model under noisy updating regimen
# return a matrix of transitions between attractors caused by noise

# set the noise level for the program can use it
perturbationCoin<-function(noiseDegree){
  number<-runif(1,min=0,max=1)
  if (number<noiseDegree){
    coin<-TRUE
  }
  else { coin <- FALSE }
  return(coin)
}

testPerturbationCoin<-function(noiseLevel=0.41){
  noiseVal<-noiseLevel
  repetitions<-1000
  valuesRep<-rep(0,repetitions)
  for (i in 1:repetitions){
    valuesRep[i]<-perturbationCoin(noiseVal)
  }
  numbTRUES<-0
  for (i in valuesRep){
    if (i == TRUE) {
      numbTRUES<-numbTRUES+1
    }
  }
  cat("Noise Level: ",noiseLevel," Number of heads: ",numbTRUES/repetitions,"\n")
}

# testPerturbationCoin()

###########################################################################################################################

# Perturb a state
# Take an state, give it a random perturbation and update
noisyUpdate<- function(net,state,noiseDegree,type="synchronous" ){
  # first perturbs 
  for (i in 1:length(state)) {                
    coin<-perturbationCoin(noiseDegree)               
    if (coin == TRUE) {                             
      state[i]<-(state[i]+1)%%2               
    }                             
  }
  # then update perturbed state
  if (type =="synchronous"){
    return(stateTransition(net,state,type = "synchronous"))
  }
  if (type =="asynchronous"){
    return(stateTransition(net,state,type = "asynchronous"))
  }
}

# testNoisyUpdate(net)

###########################################################################

# identifies if an allready calculated attractor has been found and return
# attractor number of FALSE if no attractor is reached
findAttractor<- function(net,state,attractors,type="asynchronous",numberOfFixedPoints){
  result<-FALSE
  if (type=="synchronous"){
    numberOfAttractors<-length(attractors[2]$attractors)
  }
  if (type=="asynchronous"){
    numberOfAttractors<-numberOfFixedPoints
  }
  for ( i in 1:numberOfAttractors ) {
    for ( j in 1:length( attractors[2]$attractors[[i]][[1]] ) ){
      attractor <- attractors[2]$attractors[[i]][[1]][[j]]
      decimalState<-bitsToInt(rev(state))
      attractorFound<- attractor == decimalState
      if ( attractorFound == TRUE ) break                                         
    }
    # cat(attractor,decimalState,i,j,"\n")
    if ( attractorFound == TRUE ) { result <- i}
    if ( attractorFound == TRUE ) break
  }
  return(result)
}

validateFindAttractor<-function(attractors){
  for ( i in 1:length( attractors[2]$attractors ) ) {
    for ( j in 1:length( attractors[2]$attractors[[i]][[1]] ) ){
      state<-attractors[2]$attractors[[i]][[1]][j]
      cat( "Attractor", i,": ",state, "fits attractor -> ",findAttractor(net,
                                decimalToBinary(state,length(net$genes)),
                                attractors),"\n" ) 
    }
  }
}

#validateFindAttractor(attractors)

########################################################################

# Iterate n times untill some already calcuted attractor is found
# return new attractor and time necesary to reach it
randomWalk<-function(net,state,attractors,noiseDegree,numberOfIterations,
                     type="asynchronous",numberOfFixedPoints){
  # first iteration
  i<-1
  if (type == "synchronous"){
    newState<-noisyUpdate(net,state,noiseDegree,type="synchronous")
    attractorFound <- findAttractor(net,newState,attractors,type="synchronous")
  }
  if (type=="asynchronous"){
    newState<-noisyUpdate(net,state,noiseDegree,type="asynchronous")
    attractorFound <- findAttractor(net,newState,attractors,type="asynchronous",
                                    numberOfFixedPoints = numberOfFixedPoints)
  }
  while ( ( i < numberOfIterations ) & attractorFound==FALSE ) {
    if (type=="synchronous"){
      newState<-noisyUpdate(net,newState,noiseDegree,type="synchronous")
      attractorFound <- findAttractor(net,newState,attractors,type="synchronous")
    }
    if (type=="asynchronous"){
      newState<-noisyUpdate(net,newState,noiseDegree,type="asynchronous")
      attractorFound <- findAttractor(net,newState,attractors,type="asynchronous",
                                      numberOfFixedPoints = numberOfFixedPoints)
    }
    #print(c(attractorFound,i))
    i<-i+1
  }
  return(attractorFound)
}



#################################################################################################
# transtionsFrecuencies() Repeat previous funtions n times and return frecuency of
# transition events.
# Take an an state (attractors), simulate a random walk and return frecuency
# vector of transitions to other attractors
transitionsFrecuencies<-function(net,state,attractors,noiseDegree,numberOfIterations,
                                 numberOfWalks,type="synchronous",numberOfFixedPoints){
  if (type=="synchronous"){
    transitions<-vector(mode="numeric",length(attractors$attractors))
    for (n in 1:numberOfWalks){
      walkResult<-randomWalk(net,state,attractors,noiseDegree,
                             numberOfIterations,type = "synchronous")
      if ( walkResult != FALSE ) { 
        transitions[walkResult]<-transitions[walkResult] + 1  
      }
    }
  }
  if (type=="asynchronous"){
    transitions<-vector(mode="numeric",numberOfFixedPoints)
    for (n in 1:numberOfWalks){
      walkResult<-randomWalk(net,state,attractors,noiseDegree,
                             numberOfIterations,type = "asynchronous",
                             numberOfFixedPoints = numberOfFixedPoints)
      if ( walkResult != FALSE ) { 
        transitions[walkResult]<-transitions[walkResult] + 1  
      }
    }
  }
  return((transitions)/numberOfWalks)
}

validateTransitionsFrecuencies<-function(net,attractors,numberOfWalks,
                                         noiseDegree,numberOfIterations,
                                         type,numberOfFixedPoints) {
  for ( i in 1:length( attractors[2]$attractors ) ) {
    for ( j in 1:length( attractors[2]$attractors[[i]][[1]] ) ){
      state<-decimalToBinary(attractors[2]$attractors[[i]][[1]][j],length(net$genes))
      transitions<-transitionsFrecuencies(net,state,attractors,noiseDegree,
                                          numberOfIterations,
                                          numberOfWalks,type,numberOfFixedPoints)
      cat("attractor", i,"\n",length(transitions),"\n")
      cat( transitions,"\n" ) 
    }
  }
}

#validateTransitionsFrecuencies(net,attractors,numberOfWalks=5,noiseDegree = 0.001,
#                              numberOfIterations = 10000,type="synchronous",
#                              numberOfFixedPoints = getNumberOfFixedPoints(attractors))

#result<-transitionsFrecuencies(net,state,attractors,noiseDegree,numberOfIterations,numberOfWalks)
#cat(result)

######################################################################################################################
# Construct a matrix of transitions between attractors caused by noise (synchronous regimen). 
# Where the entry A_ij represent the frecuency of transitions from the
# i-th to j-th attractors.
getTransitionsMatrixSynchronous<-function(net,attractors,numberOfWalks,noiseDegree,numberOfUpdates){
  transitionsMatrix<-matrix(0,length(attractors$attractors),length(attractors$attractors))
  for (i in 1:length(attractors$attractors)){
    periodSize<-length(attractors$attractors[[i]][[1]])
    #cat(i,periodSize,"\n")
    if (1==periodSize) {
      state<-decimalToBinary(attractors$attractors[[i]][[1]][1],length(net$genes))
      frecuencies<-transitionsFrecuencies(net,state,attractors,noiseDegree,
                                          numberOfIterations = numberOfUpdates,
                                          numberOfWalks,type = "synchronous")
      transitionsMatrix[i,]<-frecuencies
      cat("transition frecuencies of attractor",i,"was calculated","\n")
    }
    else {
      representative<-sample(1:periodSize,1)
      state<-decimalToBinary(attractors$attractors[[i]][[1]][representative],length(net$genes))
      frecuencies<-transitionsFrecuencies(net,state,attractors,noiseDegree,
                                          numberOfIterations = numberOfUpdates,
                                          numberOfWalks,type = "synchronous")
      transitionsMatrix[i,]<-frecuencies
      cat("transition frecuencies of attractor",i,"was calculated","\n")
    }
  }
  #heatmap(transitionsMatrix,Colv = NA,Rowv = NA)
  return(transitionsMatrix)
}


######################################################################################################################
# Construct a matrix of transitions between attractors caused by noise (asynchronous regumen). 
# Where the entry A_ij represent the frecuency of transitions from the
# i-th to j-th attractors.
getTransitionsMatrixAsynchronous<-function(net,attractors,numberOfWalks,noiseDegree,
                                           numberOfUpdates){
  numberOfFixedPoints<-getNumberOfFixedPoints(attractors)
  transitionsMatrix<-matrix(0,numberOfFixedPoints,numberOfFixedPoints)
  for (i in 1:numberOfFixedPoints) {
    state<-decimalToBinary(attractors$attractors[[i]][[1]][1],length(net$genes))
    frecuencies<-transitionsFrecuencies(net,state,attractors,
                                        noiseDegree=noiseDegree,
                                        numberOfIterations=numberOfUpdates,
                                        numberOfWalks=numberOfWalks,
                                        type = "asynchronous",
                                        numberOfFixedPoints=numberOfFixedPoints)
    cat(frecuencies,"\n")
    transitionsMatrix[i,]<-frecuencies
    cat("transition frecuencies of attractor",i,"was calculated","\n")
  }
  #heatmap(transitionsMatrix,Colv = NA,Rowv = NA)
  return(transitionsMatrix)
}

###############################################################################################

getTransitionsMatrix<-function(net,attractors,numberOfWalks,numberOfUpdates,noiseDegree,
                               type="asynchronous",tag="no",patternsList,simplified="no"){
  
  if (type=="synchronous"){
    cfm<-getTransitionsMatrixSynchronous(net,attractors,numberOfWalks=numberOfWalks,
                                         numberOfUpdates=numberOfUpdates,
                                         noiseDegree=noiseDegree)
  }
  if (type=="asynchronous"){
    cfm<-getTransitionsMatrixAsynchronous(net,attractors,numberOfWalks=numberOfWalks,
                                          numberOfUpdates=numberOfUpdates,
                                           noiseDegree=noiseDegree)
  }
  if (tag == "yes") {
    attractorsMatrix<-tagAttractorsMatrix(net,attractors,patternsList)
    cfm<-cfm[ 1:length(attractorsMatrix[1,]) , 1:length(attractorsMatrix[1,]) ]
    colnames(cfm)<-colnames(attractorsMatrix)
    rownames(cfm)<-colnames(attractorsMatrix)
    if (simplified=="yes") {
      cfm<-simplifyCellFateMap(cfm)
    }
  }
  return(cfm)
}

#getTransitionsMatrix(net,attractors,numberOfWalks = 30,numberOfUpdates = 30,
#                     noiseDegree = 0.03,type = "asynchronous",tag="yes",patternsList = phenotypes)


