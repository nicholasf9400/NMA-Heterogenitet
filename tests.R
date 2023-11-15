# Test fil for funktioner

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

