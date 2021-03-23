#index function is used to index in cases where from > to situation might arise
#in this case an empty index is returned
index<-function(from,to) {
  if ( from > to ) {
    R<-c()
  } else {
    R<-from:to
  }
  R
}

negindex<-function(from,to) {
  if ( from < to ) {
    R<-c()
  } else {
    R<-from:to
  }
  R
}