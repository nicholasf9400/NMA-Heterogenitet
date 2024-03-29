IndComp <- function(net, t1, t2, effect){
  
  require(igraph, quietly = T)
  require(dplyr, quietly = T)
  
  # Check class of input
  if (class(net) != 'netmeta') {
    stop('net must be netmeta object')
  }
  
  # Check if measure of effect requires log transform
  log_sm <- (net$sm == 'HR' | net$sm == 'OR' | net$sm == 'RR')
  
  # t1:t2 may be switched depending on alphabetical order
  t_list <- c(t1, t2)
  t_list <- t_list[order(t_list)]
  comp <- paste0(t_list[1],':',t_list[2]); rm(t_list)

  # Use igraph to determine all first order loops in network
  graph <- igraph::graph_from_adjacency_matrix(net$A.matrix, mode = 'undirected')
  ind.paths <- igraph::all_simple_paths(graph, 
                                        from = t1, 
                                        to = t2, 
                                        cutoff = 2)
  
  # Process output from the igraph package
  ind <- sapply(ind.paths, 
                function(x){
                  as.numeric(length(x) == 3)
                }
  )
  
  
  # Extract and compute indirect effect for all first order loops in network
  ind.res <- matrix(NA, sum(ind), 4)
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
      # effect. NB: Does not take multiarm studies into account. 
      # ind.effect <- round(dr1 - dr2, 2)
      # seind.effect <- sqrt(dr1.se^2 + dr2.se^2)
      
      # Format data to perform truncated NMA for I2 estimation
      ind_data <- net$data %>% 
        filter(((treat1 == t1) & (treat2 == intermediate.int)) | 
                 (treat1 == intermediate.int & treat2 == t2))
      
      
      ind_data$id <- 1:nrow(ind_data)

      # Correct class
      class(ind_data) <- c("pairwise", "data.frame")
      
      # Do NMA to get I
      net_temp <- netmeta(TE = TE,
                          seTE = seTE,
                          treat1 = treat1, 
                          treat2 = treat2,
                          studlab = id,
                          data = ind_data)
      
      ind.effect <- get(paste0('TE.', effect), net_temp)[t1, t2]
      seind.effect <- get(paste0('seTE.', effect), net_temp)[t1, t2]
      I2.temp <- get('I2', net_temp)

      
      # Collect data for output
      ind.res[k,] <- c(intermediate.int, ind.effect, seind.effect, round(I2.temp, 2))
      k <- k+1
    }
  }

  # Remove edges representing direct evidence
  dat <- net$data %>% filter(!((treat1 == t1) & (treat2 == t2))) %>% 
    filter(!(treat1 == t2 & treat2 == t1))
  
  for (i in 1:nrow(ind.res)){
    
    dat <- dat %>% filter(!(treat1 == t1 & treat2 == ind.res[i, 1])) %>%
      filter(!(treat1 == ind.res[i, 1] & treat2 == t1)) %>%
      filter(!(treat1 == t2 & treat2 == ind.res[i, 1])) %>%
      filter(!(treat1 == ind.res[i, 1] & treat2 == t2))
  }
  
  # Do nma ignoring multiarms studies
  dat$new_id <- 1:nrow(dat)
  net_temp <- try(netmeta(TE, seTE, treat1, treat2, studlab = new_id, data = dat))
  
  # If NMA successful, continue
  if (class(net_temp) == 'netmeta') {
    if (any(net_temp$trts == t1) & any(net_temp$trts == t2)) {
      
      high_order_exists <- T
      
      nma_est <- get(paste0('TE.', effect), net_temp)[t1, t2]
      nma_se <- get(paste0('seTE.', effect), net_temp)[t1, t2]
      nma_lower <- get(paste0('lower.', effect), net_temp)[t1, t2]
      nma_upper <- get(paste0('upper.', effect), net_temp)[t1, t2]
      
      if (log_sm) {
        
        higher_ind_out <- c('Higher order indirect effect (95% CI)', '', '',
                            paste0(round(exp(nma_est), 2), '(', round(exp(nma_lower), 2), '-', round(exp(nma_upper), 2), ')'),
                            nma_se, 
                            exp(nma_lower), 
                            exp(nma_upper), 
                            exp(nma_est), 
                            '')
        
        
      }else{
        
        higher_ind_out <- c('Higher order indirect effect (95% CI)', '', '',
                            paste0(round(nma_est, 2), '(', round(nma_lower, 2), '-', round(nma_upper, 2), ')'),
                            nma_se, nma_lower, nma_upper, nma_est, '')
      }
      
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
      paste0(round(exp(as.numeric(ind.res[,2])), 2), ' (', 
             round(exp(as.numeric(ind.res[,2]) - 1.96 * as.numeric(ind.res[,3])), 2), '-', 
             round(exp(as.numeric(ind.res[,2]) + 1.96 * as.numeric(ind.res[,3])), 2), ')'),
      ind.res[,3],
      round(exp(as.numeric(ind.res[,2]) - 1.96 * as.numeric(ind.res[,3])), 2),
      round(exp(as.numeric(ind.res[,2]) + 1.96 * as.numeric(ind.res[,3])), 2),
      exp(as.numeric(ind.res[,2])),
      paste0(ind.res[,4],'%')
    )
    
    tot.ind <- data.frame('Total Indrect Effect (95% CI)','','', 
                          paste0(round(exp(ind.TE), 2), ' (', 
                                 round(exp(ind.TE - 1.96 * ind.seTE), 2), '-', 
                                 round(exp(ind.TE + 1.96 * ind.seTE), 2), ')'), 
                          ind.seTE,
                          exp(ind.TE - 1.96 * ind.seTE), 
                          exp(ind.TE + 1.96 * ind.seTE), 
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
  
  # If higher order indirect effects are found, add them to output.
  if (exists('high_order_exists')) {
    colnames(higher_ind_out) <- colnames(out)
    
    out <- rbind(out, higher_ind_out)
  }
  
  return(out)
}
