# FUNCTION FOR INVESTIGATING HETEROGENITY

## PRELIMENARIES

library(netmeta)
library(forestplot)
library(dplyr)
library(forestploter)

source('DirectComp.R')
source('IndComp.R')


data("Gurusamy2011")

d <- pairwise(treatment, death, n, data = Gurusamy2011, studlab = study, incr = 0.5, allstudies = T)

net <- netmeta(d, sm = 'RR', reference.group = 'Control/Placebo')


#Which direct effect contribute to a arbitrary mixed estimate
#Inputs
t1 <- net$treat1[2]
t2 <- net$treat2[3]
effect <- 'random'
net <- net

a <- decomp.design(net)
a$


#___________________________________________________________________________________________________




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

#TOTAL EFFECT



#_______________________________________________________________________________

## COMPLETE FUNCTION

#T1:T2 may be switched depends on alphabetical order
t_list <- c(t1,t2)
t_list <- t_list[order(t_list)]
comp <- paste0(t_list[1],':',t_list[2]); rm(t_list)

#Check if comparison is valid
if (!(comp %in% net$comparisons)) {
  stop('Comparison is not valid! Check spelling and try again.')
}

#Check if comparison is only direct evidence

prop.dir <- get(paste0('prop.direct.', effect), net)

if (prop.dir[comp]==1) {
  cat('Comparison is made using direct evidence only...\n')
  dir.effect <- DirectComp(net, t1, t2, effect)
}else{
  if (prop.dir[comp]==0) {
    cat('Comparison only consists of indirect evidence')
    ind.effect <- IndComb(net,t1,t2,effect)
  }else{
    
    a <- DirectComp(net,t1,t2,effect)
    dir.effect <- a$ef
    ind.effect <- IndComb(net,t1,t2,effect)
    tot.effect <- TotalEffect(t1, t2, net, effect)
    
    min.scale <- min(as.numeric(dir.effect$lower), as.numeric(ind.effect$lower), as.numeric(tot.effect$tot.lower))
    max.scale <- max(as.numeric(dir.effect$upper), as.numeric(ind.effect$upper), as.numeric(tot.effect$tot.upper))
    
    #Plotting
    theme <- forest_theme(
      base_size = 10,
      # Confidence interval point shape, line type/color/width
      ci_pch = 15,
      ci_col = "#762a83",
      ci_fill = "black",
      ci_alpha = 0.8,
      ci_lty = 1,
      ci_lwd = 1.5,
      ci_Theight = 0.2, # Set an T end at the end of CI 
      # Reference line width/type/color
      refline_lwd = 1,
      refline_lty = "dashed",
      refline_col = "grey20",
      # Vertical line width/type/color
      vertline_lwd = 1,
      vertline_lty = "dashed",
      vertline_col = "grey20",
      # Change summary color for filling and borders
      summary_fill = "#4575b4",
      summary_col = "#4575b4",
      # Footnote font size/face/color
      footnote_cex = 0.6,
      footnote_fontface = "italic",
      footnote_col = "blue")
    
    
    dt <- data.frame(cbind(dir.effect[,1:4],rep(' ', nrow(dir.effect)), dir.effect[,5:9]))
    colnames(dt) <- c('Study', 'Treatment 1', 'Treatment 2', 'Weight', 'Outcome', 'RR (95% CI)', 'se', 'lower', 'upper', 'est')
    
    p.dir <- forestploter::forest(dt[,1:6], 
                                  est = as.numeric(dt$est),
                                  lower=as.numeric(dt$lower),
                                  upper=as.numeric(dt$upper),
                                  ci_column = 5,
                                  title = 'Direct effects',
                                  theme = theme,
                                  xlim = c(min.scale, max.scale),
                                  is_summary = c(rep(F, nrow(dt)-1), T))
    
    p.dir <- edit_plot(p.dir, 
                       row = nrow(dt),
                       gp = gpar(fontface = 'bold'))
    
    p.dir <- add_grob(p.dir, row = nrow(dt)+1, col = 1:6, order = "backgroud", 
                      gb_fn = gridtext::richtext_grob,
                      text = a$fn,
                      gp = gpar(fontsize = 8),
                      hjust = 0, vjust = 1,
                      halign = 0, valign = 1,
                      x = unit(0,'npc'), y = unit(1, "npc"))
    
    #Indirect effect
    it <- data.frame(cbind(ind.effect[,1:3], rep(' ', nrow(ind.effect)), rep(' ', nrow(ind.effect)), ind.effect[,4:8]))
    colnames(it) <- c('Treatment 1', 'Treatment 2', 'via', ' ', 'Outcome', 'RR (95% CI)', 'se', 'lower', 'upper', 'est')
    
    p.ind <- forestploter::forest(it[,1:6], 
                                  est = as.numeric(it$est),
                                  lower = as.numeric(it$lower),
                                  upper = as.numeric(it$upper),
                                  ci_column = 5, 
                                  title = 'First Order Indirect effects',
                                  theme = theme,
                                  xlim = c(min.scale, max.scale),
                                  is_summary = c(rep(F, nrow(it)-1), T))
    
    p.ind <- edit_plot(p.ind,
                       row = nrow(it),
                       gp = gpar(fontface = 'bold'))
    
    tt <- data.frame(cbind(tot.effect[,1:2], '', '', '', tot.effect[,4:8]))
    colnames(tt) <- c('Treatment 1', 'Treatment 2', '','', 'Outcome', 'RR (95% CI)', 'se', 'lower', 'upper', 'est')
    
    p.tot <- forestploter::forest(tt[,1:6],
                                  est = as.numeric(tt$est),
                                  lower = as.numeric(tt$lower),
                                  upper = as.numeric(tt$upper),
                                  ci_column = 5, 
                                  title = 'Total effects',
                                  theme = theme,
                                  xlim = c(min.scale, max.scale),
                                  is_summary = T)
    
    
    grid:grid.newpage()
    grid:grid.draw(rbind(p.dir, p.ind, p.tot))
  }
}







