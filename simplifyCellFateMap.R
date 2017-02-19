# condenseMatrix() simplify the matrix of transitions between attractors
# merges rows and colums with the same names adding the values of the
# merged entries.

##################################################################################################

# Merge the columns of the tagged cellFateMap matrix that have the same names, that is, if
# applied before tag attractors, it merge attractors that posses the same patterns (class).
mergeCol<-function(taggedCellFateMap,attPatterns){
  simplifiedCFM<-matrix(0,length(taggedCellFateMap[,1]),length(attPatterns))
  colnames(simplifiedCFM)<-names(attPatterns)
  rownames(simplifiedCFM)<-rownames(taggedCellFateMap)
  for (m in 1:length(attPatterns)){
    for (n in 1:length(simplifiedCFM[,1])){
      simplifiedCFM[n,m]<- sum( taggedCellFateMap[n,attPatterns[[m]]] ) 
    }
  }
  return(simplifiedCFM)
}

##################################################################################################

# simplify a matrix with repeated names.
# Figure 3B. Simplified cell fate map. Merge attractors that posses the same patterns (classes) and define a transition
# between classes as any transtion from attractor from one class to the other class.
simplifyCellFateMap<-function(taggedCellFateMap){
  nombres<-colnames(taggedCellFateMap)
  attPatterns<-c() 
  # record non repeated names
  for (i in 1:length(nombres)){
    if ( ! (nombres[i]%in%attPatterns) ){
      attPatterns<-c(attPatterns,nombres[i])
    }
  }
  # record the column indexes in which every pattern appear in the matrix of attractors
  # return a list with patterns and the indexes in which they appear in the matri
  counters<-list()
  for (k in 1:length(attPatterns)){
    indexes<-c()
    for (l in 1:length(nombres)){
      if ( attPatterns[k] == nombres[l]){
        indexes<-c(indexes,l)
      }
    }
    counters[[k]]<-indexes
    #cat(k,": ",counters[[k]],"\n")
  }
  names(counters)<-attPatterns
  # simplify columns
  cfm<-mergeCol(taggedCellFateMap,counters)
  # simplify rows using the above defined funtion
  cfm<-t(cfm)
  cfm<-mergeCol(cfm,counters)
  cfm<-t(cfm)
  rownames(cfm)<-attPatterns
  colnames(cfm)<-attPatterns
  for (i in 1:length(cfm[,1])){
    if ( sum(cfm[i,]) > 0 ) { 
    cfm[i,]<-cfm[i,]/sum(cfm[i,])
    }
  }
  return(cfm)
}



