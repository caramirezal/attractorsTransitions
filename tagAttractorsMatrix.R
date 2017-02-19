# August 26. Written by Carlos Ramirez.

# tagAttractors() take the attractors of the model and a list of pattern definitions. Macth both attractors to
# patterns and returns a matrix with tagged attractors
# find the pattern defined at the begining of the script in the attractors of the model. Return a matrix of the tagged
# attractors.
tagAttractorsMatrix<-function(net,attractors,patternsList){
  # traverse attractors BoolNet object and obtain binary vectors of the vectors
  numberOfFixedPoints<-0
  for ( i in 1:length(attractors$attractors) ) {
    periodSize<-length(attractors[2]$attractors[[i]][[1]])
    if ( 1 == periodSize ) {
      numberOfFixedPoints<-numberOfFixedPoints + 1  
    }
  }
  attractorsMatrix<-matrix(0,length(net$genes),numberOfFixedPoints)
  for (i in 1:numberOfFixedPoints) {
    state<-attractors[2]$attractors[[i]][[1]][[1]]
    state<-decimalToBinary( state,length(net$genes) )
    attractorsMatrix[,i]<-state
  }
  rownames(attractorsMatrix)<-net$genes
  namesList<-list()
  for (j in 1:length(attractorsMatrix[1,])) {
    namesList[[j]]<-"NS"
    for (i in 1:length(patternsList)){
      markers<-names(patternsList[[i]])
      if ( all(patternsList[[i]]==attractorsMatrix[,j][markers]) ) {
        namesList[[j]]<-c(namesList[[j]],names(patternsList)[i])
      }
    }
  }
  namesVect<-rep("NS",length(namesList))
  for (k in 1:length(namesList)){
    if (length(namesList[[k]]) > 1 ) {
      namesList[[k]]<-namesList[[k]][2:length(namesList[[k]])]
      namesList[[k]]<-paste(substr(namesList[[k]],1,3),collapse = "/")
      namesVect[k]<-namesList[[k]]
    } 
  }
  rownames(attractorsMatrix)<-net$genes
  colnames(attractorsMatrix)<-namesVect
  return(attractorsMatrix)
}