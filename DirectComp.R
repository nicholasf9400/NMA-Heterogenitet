DirectComp <- function(net, comp, effect){
  require(dplyr)
  require(netmeta)
  
  #T1:T2 may be switched depends on alphabetical order
  t_list <- c(t1, t2)
  t_list <- t_list[order(t_list)]
  comp <- paste0(t_list[1],':',t_list[2]); rm(t_list)
  
  #What proportion is made from direct evidence?
  prop.dir <- get(paste0('prop.direct.', effect), net)
  
  #Check if comparison is in the network
  if (!(comp %in% net$comparisons)) {
    stop('Comparison is not valid! Check spelling and try again.')
  }
  if (prop.dir[comp]==0) {
    stop('Comparison does not contain direct evidence')
  }
  
  
  #Extract data from netmeta object
  dir.TE <- get(paste0('TE.direct.', effect), net)[ t1, t2 ]
  sm <- net$sm
  d <- net$data %>% filter((treat1==t1 & treat2==t2) | (treat1==t2 & treat2==t1))
  
  
  MA <- metagen(TE, seTE, studlab, data=d, sm=sm)
  MA.dir <- get(paste0('TE.', effect), MA)
  MA.SE.dir <- get(paste0('seTE.', effect), MA)
  
  #Compute weight of each study in direct effect estimate
  dir.w <- prop.dir[comp]
  s.MA <- get(paste0('w.', effect), summary(MA))
  weights <- sapply(1:length(s.MA), 
                    function(x){
                      s.MA[x]/sum(s.MA)*dir.w
                      }
                    )
  
  
  #Gather heterogeneity statistics for footnote
  fn <- paste0("Heterogeneity: &tau;<sup>2</sup> = ", round(MA$tau2,2), "; &chi;<sup>2</sup> = ", round(MA$Q,2), ", df =", round(MA$df.Q,2), " (P=", round(MA$pval.Q,2), "); I<sup>2</sup> = ", round(MA$I2,2), "%",
               "<br><span style ='color:blue'>**Total for overall effect:**</span> Z =", round(get(paste0('zval.', effect), MA),2), " (P =", round(get(paste0('pval.', effect), MA),2), ")")
  
  total <- c(MA.dir, MA.SE.dir)
  
  #Compile output for direct effect
  if((sm == "RR")|(sm =="OR")|(sm == 'HR')){
    ef <- data.frame(
      d$studlab,
      d$treat1,
      d$treat2,
      paste0(round(weights*100,2),'%'),
      paste0(exp(round(d$TE,2)), ' (', exp(round(d$TE -1.96*d$seTE,2)), '-', exp(round(d$TE +1.96*d$seTE,2)), ')'),
      exp(round(d$seTE,2)),
      exp(round(d$TE -1.96*d$seTE, 2)),
      exp(round(d$TE +1.96*d$seTE, 2)),
      exp(round(d$TE, 2))
    )
    ef <- rbind(ef, c('Total (95% CI)', '', '', paste0(round(dir.w*100,2),'%'), 
                      paste0(exp(round(total[1],2)), ' (', exp(round(total[1] - 1.96 * total[2],2)),'-', exp(round(total[1] + 1.96 * total[2],2)), ')'), total[2], exp(total[1] - 1.96 * total[2]), exp(total[1] + 1.96 * total[2]), exp(total[1])))
    
  }else{
    ef <- data.frame(
      d$studlab,
      d$treat1,
      d$treat2,
      paste0(round(weights*100,2),'%'),
      paste0(round(d$TE,2), ' (', round(d$TE -1.96*d$seTE,2), '-', round(d$TE +1.96*d$seTE,2), ')'),
      round(d$seTE,2),
      round(d$TE -1.96*d$seTE, 2),
      round(d$TE +1.96*d$seTE, 2),
      round(d$TE, 2)
      )
  }
  colnames(ef) <- c('Study', 'Treatment 1', 'Treatment 2', 'Weight', paste0(sm, ' (95% CI)'), 'se', 'lower', 'upper', 'est')
  
  #Reformat outputs as lists
  out <- list(ef= ef, fn = fn)
  
  return(out)
}
#-------------------------------------------------------------------------------

#Testing
library(netmeta)
data("Gurusamy2011")

d <- pairwise(treatment, death, n, data = Gurusamy2011, studlab = study, incr = 0.5, allstudies = T)

net <- netmeta(d, sm = 'RR', reference.group = 'Control/Placebo')


#Which direct effect contribute to a arbitrary mixed estimate
#Inputs
t1 <- net$treat1[2]
t2 <- net$treat2[3]
effect <- 'random'
net <- net

test <- DirectComp(net, t1, t2, 'random')

