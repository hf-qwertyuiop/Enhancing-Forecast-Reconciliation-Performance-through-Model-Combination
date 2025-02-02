setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Read in and process the S-matrix
readS<-function(fileS){
  S<-read.csv(fileS)[,-1]
  nb<-ncol(S);na<-nrow(S)-nb
  for (i in 1:na){
    rownames(S)[i]<-paste0('U',as.character(i))
  }
  for (i in 1:nb){
    rownames(S)[na+i]<-paste0('B',as.character(i))
    colnames(S)[i]<-paste0('B',as.character(i))
  }
  S<-as.matrix(S)
  return (list(S,na+nb,na))
}

#Find a linearly independent C matrix
inde_C<-function(C,J,C_n,n,na){
  if (nrow(C)==0){
    C[1,]<-C_n
    J[which(C_n != 0)[1],]<-C_n
  }else{
    J_n<-C_n
    while (!all(J_n==0)){
      j<-which(J_n != 0)[1]
      if (is.na(J[j,j]) | J[j,j]==0){
        J[j,]<-J_n
        C[nrow(C)+1,]<-C_n
        break
      }else{
        J_n<-J_n-(sign(J_n[j])*sign(J[j,j]))[[1]]* J[j,]
      }
    }
  }
  if (nrow(C)==na+1) {return(list(C,J,TRUE))}else{return (list(C,J,FALSE))}
  
}

#Find all candidates
find_C<-function(S,bti,n,na){
  #Preparatory work
  C<-data.frame(matrix(numeric(0), ncol = n, nrow = 0));#used for store
  J<-data.frame(matrix(numeric(0), ncol = n, nrow = 0));#used for elimination
  Cn<-data.frame(matrix(numeric(0), ncol = nrow(S), nrow = 0));colnames(Cn)<-rownames(S)
  bottom<-colnames(S)
  #Start looking for candidates for the bottom node bti
  #direct candidate
  newCn_index<-nrow(Cn)+1
  Cn[newCn_index,]<-0
  Cn[newCn_index,bti]<-1
  
  C_n<-Cn[newCn_index,]
  results<-inde_C(C,J,C_n,n,na)
  C<-results[[1]];J<-results[[2]];signal<-results[[3]]
  if (signal==TRUE) return(C)#If the rank is full then we can find the next bottom node
  
  ##Other candidates
  ancestor <- S[ S%*%t(S[bti,,drop=FALSE])==1 &
                   rownames(S)!=bti, , drop=FALSE] #Fetch the ancestor node of bti (excluding itself)
  
  displace<-data.frame(matrix(0, ncol = n, nrow = na,dimnames = list(rownames(S)[1:na],rownames(S))))
  for (i in 1:na){
    displace[i,i]<- -1
    displace[i,(na+1):n]<- S[i,]
  }
  record<-c(rep(FALSE,times=na));names(record)<-rownames(S)[1:na]
  
  for (aci in rev(rownames(ancestor))){ #For each ancestor node (placed first to be subtracted)
    record[aci]<-TRUE
    
    minus  <- S[ S%*%S[bti,]==0 &
                   S%*%ancestor[aci,]!=0, , drop=FALSE] #Take the bti side nodes (subtracted from the combination)
    minus_upper<-minus[rowSums(minus)> 1,,drop=FALSE]#Upper nodes to be subtracted/used for composition
    minus_lower<-minus[rowSums(minus)<=1,,drop=FALSE]#bottom node to be subtracted/used for composition
    tmp<-ancestor[aci,,drop=FALSE];stack<-c();pointer<-1
    
    newCn_index<-nrow(Cn)+1
    Cn[newCn_index,]=0 #other nodes is 0
    Cn[newCn_index,aci]<- 1 #this ancestor is 1
    Cn[newCn_index,c(rownames(minus_lower)[tmp%*%t(minus_lower)==1])]<- -1 #Subtract the bottom node (set to -1)
    baseline<-Cn[newCn_index,]
    
    C_n<-Cn[newCn_index,]
    results<-inde_C(C,J,C_n,n,na)
    C<-results[[1]];J<-results[[2]];signal<-results[[3]]
    
    for (disp in rownames(minus_upper)){
      if (record[disp]==TRUE){
        next
      }else{
        record[disp]<-TRUE
      }
      newCn_index<-nrow(Cn)+1
      Cn[newCn_index,]<-baseline+displace[disp,]
      
      C_n<-Cn[newCn_index,]
      results<-inde_C(C,J,C_n,n,na)
      C<-results[[1]];J<-results[[2]];signal<-results[[3]]
      if (signal==TRUE) return(C)
    }
  }
  return (C)
}