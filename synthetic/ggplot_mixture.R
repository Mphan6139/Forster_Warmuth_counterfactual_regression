library(latex2exp)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot, scales)
load('plotdata_mixture_2000.rdata')


tmp = res2/res2[,6]
colnames(tmp)[1] = "alpha"
tmp = data.frame(tmp)
tmp[,1] = res2[,1]

df1 = tmp %>%
  #select(alpha, forster_poly, ls_poly,forster_ns,ls_ns, forster_bs, ls_bs) %>%
  gather(key = "Type", value = "value", -alpha)
df1[df1=="forster_poly"] = "FW_poly"
df1[df1=="forster_bs"] = "FW_bs"
df1[df1=="forster_ns"] = "FW_ns"

final1 = ggplot(df1, aes(x = alpha, y = value)) +
  geom_line(aes(color = Type, linetype = Type), linewidth = 2)+
  scale_linetype_manual(values = c(rep("solid",3), rep("dashed", 3)))+
  scale_color_manual(values=c(rep(hue_pal()(6)[1:3],2)))+
  geom_point(aes(color = Type, shape = Type), size = 5)+
  #ylab("n*MSE")+
  ylab("Ratio of MSE")+
  xlab(TeX(r'( $n=2000$)'))+
  theme(strip.text = element_text(size = 50, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        legend.position="top",
        legend.title = element_text(size=50,face="bold"),
        legend.text = element_text(size=40),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x=element_text(size=60,face="bold", vjust = -0.6),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=50),
        axis.text.y = element_text(face="bold",size=50),
        axis.title.y =  element_text(face="bold",size=50)
        #strip.text.x = element_blank(),
  )


load('plotdata_mixture_400.rdata')


tmp = res2/res2[,6]
colnames(tmp)[1] = "alpha"
tmp = data.frame(tmp)
tmp[,1] = res2[,1]

df2 = tmp %>%
  #select(alpha, forster_poly, ls_poly,forster_ns,ls_ns, forster_bs, ls_bs) %>%
  gather(key = "Type", value = "value", -alpha)
df2[df2=="forster_poly"] = "FW_poly"
df2[df2=="forster_bs"] = "FW_bs"
df2[df2=="forster_ns"] = "FW_ns"

final2 = ggplot(df2, aes(x = alpha, y = value)) +
  geom_line(aes(color = Type, linetype = Type), linewidth = 2)+
  scale_linetype_manual(values = c(rep("solid",3), rep("dashed", 3)))+
  scale_color_manual(values=c(rep(hue_pal()(6)[1:3],2)))+
  geom_point(aes(color = Type, shape = Type), size = 5)+
  #ylab("n*MSE")+
  ylab("Ratio of MSE")+
  xlab(TeX(r'($n=400$$)'))+
  theme(strip.text = element_text(size = 50, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        legend.position="top",
        legend.title = element_text(size=50,face="bold"),
        legend.text = element_text(size=50),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x=element_text(size=60,face="bold", vjust = 0.6),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=50),
        axis.text.y = element_text(face="bold",size=50),
        axis.title.y =  element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(4,"line"),
  )

#ggsave(("ggplot_mixture_forster_ls_400.pdf"),final, width = 16, height = 12, units = "in")

legend = get_legend(final2)
first = plot_grid(final1+theme(legend.position = "none"), final2+theme(legend.position = "none"))
plot = plot_grid(legend, first, ncol = 1, rel_heights = c(.15, 1))
x.grob <- textGrob(TeX(r'($\hat{\pi}$ convergence rate ($\alpha$ in RMSE $n^{-\alpha})$)'), gp=gpar(fontface="bold", fontsize=60))
final = grid.arrange(arrangeGrob(plot, bottom = x.grob))
ggsave(("ggplot_mixture.pdf"),final, width = 20, height = 16, units = "in")

