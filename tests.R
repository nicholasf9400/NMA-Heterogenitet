# Test fil for funktioner

library(netmeta)
data("Gurusamy2011")

d <- pairwise(treatment, death, n, data = Gurusamy2011, studlab = study)

net <- netmeta(d, sm = 'RR', reference.group = 'Control/Placebo')


#Which direct effect contribute to a arbitrary mixed estimate
#Inputs
t1 <- net$treat1[2]
t2 <- net$treat2[3]
effect <- 'random'
net <- net

source('DirectComp.R')

DirectComp(net, t1, t2, 'random')

source('IndComp.R')

IndComb(net, t1, t2, 'random')

source('TotalEffect.R')

TotalEffect(t1, t2, net, 'random')

source('DecompEffect.R')

DecompEffect(net, t1, t2, 'random')
