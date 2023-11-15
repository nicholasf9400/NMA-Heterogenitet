TotalEffect <- function(t1,t2,net,effect){
  
  tot.TE <- get(paste0('TE.', effect), net)[t1,t2]
  tot.seTE <- get(paste0('seTE.',effect), net)[t1,t2]
  tot.lower <- get(paste0('lower.', effect), net)[t1,t2]
  tot.upper <- get(paste0('upper.', effect), net)[t1,t2]

  out <- data.frame(t1,t2, '',
                    paste0(tot.TE, ' (', tot.lower, '-', tot.upper, ')'), tot.seTE, tot.lower, tot.upper, tot.TE)
  return(out)
}