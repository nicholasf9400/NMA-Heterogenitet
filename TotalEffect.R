TotalEffect <- function(t1,t2,net,effect){
  
  tot.TE <- get(paste0('TE.', effect), net)[t1,t2]
  tot.seTE <- get(paste0('seTE.',effect), net)[t1,t2]
  tot.lower <- get(paste0('lower.', effect), net)[t1,t2]
  tot.upper <- get(paste0('upper.', effect), net)[t1,t2]

  sm <- net$sm
  
  if(sm == 'HR' | sm == 'OR' | sm == 'RR'){
    out <- data.frame(t1,
                      t2, 
                      '',
                      paste0(round(exp(tot.TE), 2), 
                             ' (', 
                             round(exp(tot.lower), 2), 
                             '-', 
                             round(exp(tot.upper), 2), 
                             ')'), 
                      tot.seTE, 
                      exp(tot.lower), 
                      exp(tot.upper), 
                      exp(tot.TE))
    
  }else{
    out <- data.frame(t1,
                      t2, 
                      '',
                      paste0(round(tot.TE, 2), 
                             ' (', 
                             round(tot.lower, 2), 
                             '-', 
                             round(tot.upper, 2), 
                             ')'), 
                      tot.seTE, 
                      tot.lower, 
                      tot.upper, 
                      tot.TE)
  }
  return(out)
}