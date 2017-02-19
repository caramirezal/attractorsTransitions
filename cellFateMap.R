# cellFateMap() give all possible bitflip nodes value perturbations to attractors and evaluate if this 
# cause a transition to another attractor




# flip the value of a node 
flipNode<-function(state,index){
  state[index] <- ! state[index]  
  state
}


################################################################################################################

# give every possible single transient perturbations to attractors
# return a heatmap of transitions in a synchronous updating regimen
cellFateMapSynchronous<-function(net,attractors){
  numberOfAttractors<-length(attractors$attractors)
  cellFateMap<-matrix(0,numberOfAttractors,numberOfAttractors)           # create a matrix of length the number attractor including cycles
  for(i in 1:length(attractors$attractors)){                               # traverse every state
    periodSize <- length(attractors$attractors[[i]][[1]])
    for(j in 1:periodSize){                                    # traverse cycles
      binaryAttractor<-decimalToBinary(attractors$attractors[[i]][[1]][[j]],length(net$genes))     
      for(k in 1:length(binaryAttractor)){                     # traverse nodes
        pulsedState<-flipNode(binaryAttractor,k)               # gives a single transient perturbations
        newAttractor<-getAttractors(net,method="chosen",startStates=list(pulsedState))    # find new attractor
        newAttractor<-decimalToBinary(newAttractor$attractors[[1]][[1]][[1]],length(net$genes) )
        newAttractor<-findAttractor(net,newAttractor,attractors,
                                    type = "synchronous")   # return attractor indexes
        #cat(i,"->",newAttractor,"\n")                                  # modify this. send output to cellFateMap[findAttractor]
        if (newAttractor != FALSE) {
          cellFateMap[i,newAttractor]<-cellFateMap[i,newAttractor] + ( 1 / ( length(net$genes) * periodSize ) )
        }
      } 
    }                
  }
  #heatmap(cellFateMap,Rowv = NA,Colv = NA,symm=TRUE,col=gray.colors(40))
  return(cellFateMap)
}


###################################################################################################

# give every possible single transient perturbations to attractors
# return a heatmap of transitions in a synchronous updating regimen
cellFateMapAsynchronous<-function(net,attractors,numberOfIterations=1){
  numberOfFixedAttractors<-getNumberOfFixedPoints(attractors)
  simulations<-matrix(0,numberOfFixedAttractors,numberOfFixedAttractors)
  j<-1
  while (j <= numberOfIterations) {
    cellFateMap<-matrix(0,numberOfFixedAttractors,numberOfFixedAttractors)           # create a matrix of length the number attractor including cycles
    for(i in 1:numberOfFixedAttractors){                               # traverse every state                                   # traverse cycles
      binaryAttractor<-decimalToBinary(attractors$attractors[[i]][[1]][[1]],length(net$genes))     
      for(k in 1:length(binaryAttractor)){                     # traverse nodes
        pulsedState<-flipNode(binaryAttractor,k)               # gives a single transient perturbations
        newAttractor<-getAttractors(net,method="chosen",type = "asynchronous",
                                    startStates=list(pulsedState))    # find new attractor
        newAttractor<-decimalToBinary(newAttractor$attractors[[1]][[1]][[1]],length(net$genes) )
        newAttractor<-findAttractor(net,newAttractor,attractors,
                                    type = "asynchronous",
                                    numberOfFixedPoints =numberOfFixedAttractors )   
        if (newAttractor != FALSE) {
          cellFateMap[i,newAttractor]<-cellFateMap[i,newAttractor] + 1  
        } 
      }                
    }
    simulations<-simulations+cellFateMap
    cat("Iteration ",j," of",numberOfIterations,"\n")
    j<-j+1
  }
  cfm<- simulations / ( length(net$genes) * numberOfIterations )
  return(cfm)
}

#################################################################################################

cellFateMap<-function(net,attractors,type=c("synchronous","asynchronous"),
                      numberOfIterations=1,simpliflied="no",tag="no",patternsList){
  if (type == "synchronous"){
    cfm<-cellFateMapSynchronous(net,attractors = attractors)
  }
  if (type == "asynchronous"){
    cfm<-cellFateMapAsynchronous(net,attractors = attractors,
                           numberOfIterations=numberOfIterations)
  }
  if (tag == "yes") {
    attractorsMatrix<-tagAttractorsMatrix(net,attractors,patternsList)
    cfm<-cfm[ 1:length(attractorsMatrix[1,]) , 1:length(attractorsMatrix[1,]) ]
    colnames(cfm)<-colnames(attractorsMatrix)
    rownames(cfm)<-colnames(attractorsMatrix)
    if (simpliflied=="yes") {
      cfm<-simplifyCellFateMap(cfm)
    }
  }
  return(cfm)
}


