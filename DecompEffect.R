source('DirectComp.R')
source('IndComp.R')
source('TotalEffect.R')

DecompEffect <- function(net, t1, t2, effect){
  
  require(forestplot)
  require(forestploter)
  
  
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
      
      p.dir <- add_grob(p.dir, row = nrow(dt)+1, col = 1:6, order = "background", 
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
  
}
