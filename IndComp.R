# FIRST ORDER INDRECT EFFECT FUNCTION

IndComb <- function(net, t1, t2, effect){
  require(igraph)
  
  #T1:T2 may be switched depends on alphabetical order
  t_list <- c(t1,t2)
  t_list <- t_list[order(t_list)]
  comp <- paste0(t_list[1],':',t_list[2]); rm(t_list)
  
  graph <- igraph::graph_from_adjacency_matrix(net$A.matrix, mode = 'undirected')
  ind.paths <- igraph::all_simple_paths(graph, from = t1, to = t2, cutoff = 2)
  
  ind <- sapply(ind.paths, 
                function(x){
                  if(length(x)==3){
                    return(1)
                  }else{return(0)}
                }
  )
  ind.res <- matrix(NA, length(sum(ind==1)), 3)
  k <- 1
  for (i in 1:length(ind)){
    if(ind[i]==1){
      path <- names(ind.paths[[i]])
      intermediate.int <- path[which(path != t1 & path != t2)]
      
      dir.effect <- get(paste0('TE.direct.', effect), net)
      se.effect <- get(paste0('seTE.direct.', effect), net)
      
      
      dr1 <- dir.effect[t1, intermediate.int]
      dr1.se <- se.effect[t1, intermediate.int]
      
      dr2 <- dir.effect[intermediate.int, t2]
      dr2.se <- se.effect[intermediate.int, t2]
      
      ind.effect <- round(dr1 - dr2, 2)
      seind.effect <- sqrt(dr1.se^2 + dr2.se^2)
      
      ind.res[k,] <- c(intermediate.int, ind.effect, seind.effect)
      k <- k+1
    }
  }
  
  #TOTAL INDIRECT EFFECT
  ind.TE <- get(paste0('TE.indirect.', effect), net)[t1,t2]
  ind.seTE <- get(paste0('seTE.indirect.', effect), net)[t1,t2]
  
  out <- data.frame(
    t1,
    ind.res[,1],
    t2,
    paste0(ind.res[,2], ' (', round(as.numeric(ind.res[,2]) -1.96*as.numeric(ind.res[,3]),2), '-', round(as.numeric(ind.res[,2]) +1.96*as.numeric(ind.res[,3]),2), ')'),
    ind.res[,3],
    round(as.numeric(ind.res[,2]) -1.96*as.numeric(ind.res[,3]), 2),
    round(as.numeric(ind.res[,2]) +1.96*as.numeric(ind.res[,3]), 2),
    ind.res[,2]
  )
  colnames(out) <- c('Treatment 1', "via", 'Treatment 2', paste0(net$sm, ' (95% CI)'), 'se', 'lower', 'upper', 'est')
  
  tot.ind <- data.frame('Total Indrect Effect (95% CI)','','', paste0(ind.TE, ' (',ind.TE-1.96*ind.seTE,ind.TE+1.96*ind.seTE, ')'), ind.seTE
                        ,ind.TE-1.96*ind.seTE,ind.TE+1.96*ind.seTE, ind.TE)
  colnames(tot.ind) <- c('Treatment 1', 'Treatment 2', "via", paste0(net$sm, ' (95% CI)'), 'se', 'lower', 'upper', 'est')
  
  out <- rbind(out, tot.ind)
  
  return(out)
}
