rm(list = ls())
library(latex2exp)
library(tidyverse)
load('all_methods_uniform_2000.rdata')


res2 = final[[1]]
res2 = res2[,-c(5,7,8,9)]
res2 = res2[,-c(6,7,8,9)]

tmp = res2
#tmp = res2/res2[,6]
colnames(tmp)[1] = "alpha"
tmp = data.frame(tmp)
tmp[,1] = res2[,1]
tmp$ytype = "n*MSE"

df1 = tmp %>%
  #select(alpha, forster_poly, ls_poly,forster_ns,ls_ns, forster_bs, ls_bs) %>%
  gather(key = "Type", value = "value", -alpha, -ytype)

tmp = res2/res2[,6]
colnames(tmp)[1] = "alpha"
tmp = data.frame(tmp)
tmp[,1] = res2[,1]
tmp$ytype = "Ratio of MSE"

df2 = tmp %>%
  #select(alpha, forster_poly, ls_poly,forster_ns,ls_ns, forster_bs, ls_bs) %>%
  gather(key = "Type", value = "value", -alpha, -ytype)

df = rbind(df1, df2)
df[df=="forster_bs"] = "FW_bs"


ggplot(df, aes(x = alpha, y = value)) + 
  geom_line(aes(color = Type, linetype = Type))+
  geom_point(aes(color = Type, shape = Type))


final_plot = ggplot(df, aes(x = alpha, y = value)) +
  geom_line(aes(color = Type, linetype = Type), linewidth = 2)+
  geom_point(aes(color = Type, shape = Type), size = 5)+
  facet_grid(ytype~.,,scales = "free",switch='y')+
  xlab(TeX(r'($\hat{\pi}$ convergence rate ($\alpha$ in RMSE $n^{-\alpha})$)'))+
  theme(strip.text = element_text(size = 50, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        legend.position="top",
        legend.title = element_text(size=50,face="bold"),
        legend.text = element_text(size=50),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title=element_text(size=60,face="bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face="bold",size=50),
        axis.text.y = element_text(face="bold",size=50)
        #axis.title.y =  element_text(face="bold",size=40),
        #strip.text.x = element_blank(),
  )

ggsave(("ggplot_uniform_final.pdf"),final_plot, width = 16, height = 16, units = "in")





