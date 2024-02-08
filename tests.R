# Test fil for funktioner

library(netmeta)
data("Linde2016")

d <- pairwise(treat = list(treat1, treat2),
              event = list(resp1, resp2),
              n = list(n1, n2),
              data = Linde2016, 
              studlab = id,
              sm = 'RR')

net <- netmeta(d, sm = 'RR', reference.group = 'Placebo')

#Which direct effect contribute to a arbitrary mixed estimate
#Inputs
t1 <- "TCA"
t2 <- "Placebo"
effect <- 'random'
net <- net

source('DirectComp.R')

DirectComp(net, t1, t2, 'random')

source('IndComp.R')

IndComp(net, t1, t2, 'random')

source('TotalEffect.R')

TotalEffect(t1, t2, net, 'random')

source('DecompEffect.R')

DecompEffect(net, t1, t2, 'random')
