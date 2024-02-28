library(ggplot2)
load('unconfound_plotdata.RDATA')
est_forster_poly = unconfoundedness_data[[1]]
sd_forster_poly = unconfoundedness_data[[2]]
est_forster_bs = unconfoundedness_data[[3]]
sd_forster_bs = unconfoundedness_data[[4]]
x1_list = unconfoundedness_data[[5]]
plotdata_poly_unconfoundedness = data.frame(x=x1_list, y = est_forster_poly, lower=est_forster_poly - 1.96*sd_forster_poly, upper =est_forster_poly + 1.96*sd_forster_poly, Basis = "polynomials" )
plotdata_bs_unconfoundedness = data.frame(x=x1_list, y = est_forster_bs, lower=est_forster_bs - 1.96*sd_forster_bs, upper =est_forster_bs + 1.96*sd_forster_bs,Basis = "splines" )
plotdata_unconfoundedness = rbind(plotdata_poly_unconfoundedness,plotdata_bs_unconfoundedness)
plotdata_unconfoundedness$Setting = "Unconfoundedness"

load('proxy_plotdata.RDATA')
est_forster_poly = proxy_data[[1]]
sd_forster_poly = proxy_data[[2]]
est_forster_bs = proxy_data[[3]]
sd_forster_bs = proxy_data[[4]]
x1_list = proxy_data[[5]]
plotdata_poly_proxy = data.frame(x=x1_list, y = est_forster_poly, lower=est_forster_poly - 1.96*sd_forster_poly, upper =est_forster_poly + 1.96*sd_forster_poly, Basis = "polynomials" )
plotdata_bs_proxy = data.frame(x=x1_list, y = est_forster_bs, lower=est_forster_bs - 1.96*sd_forster_bs, upper =est_forster_bs + 1.96*sd_forster_bs,Basis = "splines" )
plotdata_proxy = rbind(plotdata_poly_proxy,plotdata_bs_proxy)
plotdata_proxy$Setting = "Proximal"
plotdata=rbind(plotdata_unconfoundedness, plotdata_proxy)

final = ggplot(data = transform(plotdata, Setting=factor(Setting, levels=c("Unconfoundedness", "Proximal"))), aes(y=y, x=x, fill = Basis, color = Basis))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=x), alpha = 0.3)+
  facet_grid(~Setting,scales = "free",switch='x')+
  ylab("CATE")+
  xlab("surv2md1")+
  theme(strip.text = element_text(size = 40, face = "bold"), 
        axis.ticks.length=unit(.4, "cm"),
        legend.position="top",
        legend.title = element_text(size=30,face="bold"),
        legend.text = element_text(size=30),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title=element_text(size=25,face="bold"),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold",size=35),
        axis.text.y = element_text(face="bold",size=35),
        axis.title.y =  element_text(face="bold",size=40),
        #strip.text.x = element_blank(),
        )+
  geom_hline(aes(yintercept=0), linetype="dashed")    
        
ggsave(("ggplot_poly_bs.pdf"),final, width = 16, height = 12, units = "in")

