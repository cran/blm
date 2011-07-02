remove.missing.blm <- function(f,data){

  var <- all.vars(f)
  missing.cases <- lapply(var,function(name){which(is.na(data[name]))})
  missing.cases <- unique(unlist(missing.cases))

  if(length(missing.cases)!=0){
    data = data[-missing.cases,]
    missing = length(missing.cases)
  }
  else{
    missing = NA
  }

  list(data=data,missing=missing)
}


remove.missing.lexpit <- function(f.linear,f.expit,data){

  var <- all.vars(f.linear)
  var <- unique(c(var,all.vars(f.expit)))
  
  missing.cases <- lapply(var,function(name){which(is.na(data[name]))})
  missing.cases <- unique(unlist(missing.cases))

  if(length(missing.cases)!=0){
    data = data[-missing.cases,]
    missing = length(missing.cases)
  }
  else{
    missing = NA
  }

  list(data=data,missing=missing)
}
