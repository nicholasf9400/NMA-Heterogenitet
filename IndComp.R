IndComb <- function(net, t1, t2, effect){
  
  require(igraph)
  require(dplyr)
  
  # Check class of input
  if (class(net) != 'netmeta') {
    stop('net must be netmeta object')
  }
  
  # Check if measure of effect requires log transform
  log_sm <- net$sm == 'HR' | net$sm == 'OR' | net$sm == 'RR'
  
  # t1:t2 may be switched depending on alphabetical order
  t_list <- c(t1,t2)
  t_list <- t_list[order(t_list)]
  comp <- paste0(t_list[1],':',t_list[2]); rm(t_list)

  # Use igraph to determine all first order loops in network
  graph <- igraph::graph_from_adjacency_matrix(net$A.matrix, mode = 'undirected')
  ind.paths <- igraph::all_simple_paths(graph, from = t1, to = t2, cutoff = 2)
  
  # Process output from the igraph package
  ind <- sapply(ind.paths, 
                function(x){
                  if(length(x)==3){
                    return(1)
                  }else{return(0)}
                }
  )
  
  # Extract and compute indirect effect for all first order loops in network
  ind.res <- matrix(NA, sum(ind==1), 4)
  k <- 1
  for (i in 1:length(ind)){
    if(ind[i] == 1){
      
      path <- names(ind.paths[[i]])
      
      # Determine intermediate treatment
      intermediate.int <- path[which(path != t1 & path != t2)]
      
      # Extract effect measures from input net
      dir.effect <- get(paste0('TE.direct.', effect), net)
      se.effect <- get(paste0('seTE.direct.', effect), net)
      
      # Extract direct estimates against intermediate treatment
      dr1 <- dir.effect[t1, intermediate.int]
      dr1.se <- se.effect[t1, intermediate.int]
      
      dr2 <- dir.effect[intermediate.int, t2]
      dr2.se <- se.effect[intermediate.int, t2]
      
      # Bucher adjusted indirect comparison as estimate for first order indirect
      # effect
      ind.effect <- round(dr1 - dr2, 2)
      seind.effect <- sqrt(dr1.se^2 + dr2.se^2)
      
      # Format data to perform truncated NMA for I2 estimation
      ind_data <- net$data %>% filter((treat1 == t1 | treat1 == t2 | treat1 == intermediate.int) & (treat2 == t1 | treat2 == t2 | treat2 == intermediate.int)) %>%
        filter(!(treat1 == t1 & treat2 == t2) | !(treat1 == t2 & treat2 == t2))
      
      # Correct class
      class(ind_data) <- c("pairwise", "data.frame")
      
      # Do NMA
      net_temp <- netmeta(TE = TE,
                          seTE = seTE,
                          treat1 = treat1, 
                          treat2 = treat2,
                          studlab = studlab,
                          data = ind_data)
      
      # Collect data for output
      ind.res[k,] <- c(intermediate.int, ind.effect, seind.effect, round(net_temp$I2, 2))
      k <- k+1
    }
  }
  
  ind.paths <- igraph::all_simple_paths(graph, from = t1, to = t2, cutoff = -1)
  
  paths <- list()
  k <- 1
  for (i in 1:length(ind.paths)) {
    if (length(ind.paths[[i]]) > 3) {
      paths[[k]] <- ind.paths[[i]]
      k <- k + 1
    }
  }
  rm(ind.paths)
  
  # Get relevant edges
  
  edges <- c()
  for (i in 1:length(paths)) {
    for (j in 1:(length(paths[[i]])-1)) {
      edges <- rbind(edges, c(names(paths[[i]][j]), names(paths[[i]][j+1])))
    }
  }
  

  
  
  
  # TOTAL INDIRECT EFFECT
  ind.TE <- get(paste0('TE.indirect.', effect), net)[t1,t2]
  ind.seTE <- get(paste0('seTE.indirect.', effect), net)[t1,t2]
  
  # Format output data
  if(log_sm){
    out <- data.frame(
      t1,
      t2,
      ind.res[,1],
      paste0(round(exp(as.numeric(ind.res[,2])), 2), ' (', round(exp(as.numeric(ind.res[,2]) -1.96*as.numeric(ind.res[,3])),2), '-', round(exp(as.numeric(ind.res[,2]) +1.96*as.numeric(ind.res[,3])),2), ')'),
      ind.res[,3],
      round(exp(as.numeric(ind.res[,2]) - 1.96*as.numeric(ind.res[,3])), 2),
      round(exp(as.numeric(ind.res[,2]) + 1.96*as.numeric(ind.res[,3])), 2),
      exp(as.numeric(ind.res[,2])),
      paste0(ind.res[,4],'%')
    )
    
    tot.ind <- data.frame('Total Indrect Effect (95% CI)','','', 
                          paste0(round(exp(ind.TE), 2), ' (', round(exp(ind.TE-1.96*ind.seTE), 2), '-', round(exp(ind.TE+1.96*ind.seTE), 2), ')'), 
                          ind.seTE,
                          exp(ind.TE - 1.96*ind.seTE), 
                          exp(ind.TE + 1.96*ind.seTE), 
                          exp(ind.TE), 
                          ' ')
  }else{
    # Output if effect measure is not on log scale
    out <- data.frame(
      t1,
      t2,
      ind.res[,1],
      paste0(round(ind.res[,2], 2), ' (', round(as.numeric(ind.res[,2]) -1.96*as.numeric(ind.res[,3]),2), '-', round(as.numeric(ind.res[,2]) +1.96*as.numeric(ind.res[,3]),2), ')'),
      ind.res[,3],
      round(as.numeric(ind.res[,2]) - 1.96*as.numeric(ind.res[,3]), 2),
      round(as.numeric(ind.res[,2]) + 1.96*as.numeric(ind.res[,3]), 2),
      ind.res[,2],
      paste0(ind.res[,4],'%')
    )
    
    tot.ind <- data.frame('Total Indrect Effect (95% CI)','','', 
                          paste0(round(ind.TE, 2), ' (',ind.TE-1.96*ind.seTE,ind.TE+1.96*ind.seTE, ')'), 
                          ind.seTE,
                          ind.TE - 1.96*ind.seTE,
                          ind.TE + 1.96*ind.seTE, 
                          ind.TE, ' ')
  }

  # Name output for compiling results
  colnames(out) <- c('Treatment 1', 'Treatment 2', "via", paste0(net$sm, ' (95% CI)'), 'se', 'lower', 'upper', 'est', "I2")
  colnames(tot.ind) <- colnames(out)
  
  # Output
  out <- rbind(out, tot.ind)
  return(out)
}
