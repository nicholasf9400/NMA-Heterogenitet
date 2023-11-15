setwd("C:/Users/NicholasFitzhugh/OneDrive - Behandlingsraadet/DACS/Projekter/NMA pakke/Code")

## PRELIMENARIES

library(netmeta)
library(forestplot)
library(dplyr)
library(forestploter)



data("Gurusamy2011")

d <- pairwise(treatment, death, n, data = Gurusamy2011, studlab = study, incr = 0.5, allstudies = T)


net <- netmeta(d, sm = 'RR', reference.group = 'Control/Placebo')

netgraph(net)


d <- d %>% filter((treat1 == "Control/Placebo" & treat2 == "Aprotonin")|(treat1 == "Aprotonin" & treat2 == "Control/Placebo")
                  |(treat1 == "Control/Placebo" & treat2 == "Tranexamic acid")|(treat1 == "Tranexamic acid" & treat2 == "Control/Placebo")|
                    (treat1 == "Aprotonin" & treat2 == "Tranexamic acid")|(treat1 == "Tranexamic acid" & treat2 == "Aprotonin"))

d <- d %>% mutate(indep = ifelse(
  ((treat1 == "Control/Placebo" & treat2 == "Tranexamic acid")|(treat1 == "Tranexamic acid" & treat2 == "Control/Placebo")),
  1, ifelse(
    (treat1 == "Aprotonin" & treat2 == "Tranexamic acid")|(treat1 == "Tranexamic acid" & treat2 == "Aprotonin"), -1, 0
  )
))

ma_bin <- metabin(event1, n1, event2, n2, studlab, d, subgroup = indep)

summary(ma_bin)


metareg(ma_bin, formula = ~-1 + .subgroup)
