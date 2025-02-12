source("cust_fun.R")
source("cust_MySL.R")

# This true psi is based on 1M sample size simulation. 
# If you change single_gen function, make sure you change psi as well for the bias evaluation.

psi = 3.164386

single_run_simulation = function(df,N,nc=3, fold_num = 1, CV,seed = NULL, k=4, bs= 'bs',choice=c(1,5,7), raw = F){
  if(is.null(seed)){
    set.seed(seed)
  }
  
  current_df = single_gen(N)
   
  y = current_df$y
  a = current_df$a
  z = current_df$z
  x1 = current_df$x1
  x2 = current_df$x2
  s <- rep(1:nc, each=N/nc)
  
  y1 = current_df$y1
  y0 = current_df$y0
  
  data_input = data.frame(x1,x2,a,z,y,s)
  
  psi_mini = mean((current_df$y1 - current_df$y0)[a==1])
  psi_full = psi
  
  res_list = rep(list(NA), fold_num)
  bias_list = rep(list(NA), fold_num)
  var_list = rep(list(NA), fold_num)
  
  for(j in 1:fold_num){
    data_input <- data_input[sample(N), ]
    data_input$s = s
    result_single = three_fold_eval(data_input, CV=CV, k=k, bs = bs, choice=choice)
    res_list[[j]] = result_single
    bias_list[[j]] = as.data.frame(lapply(result_single, function(x) mean(x)))
    var_list[[j]] = as.data.frame(lapply(result_single, function(x) var_compute(x,data_input$a,data_input$s,nc)))
  }
  
  

  bias_single = apply(do.call(rbind, bias_list), 2, median, na.rm = TRUE) - psi_full
  bias_single2 = apply(do.call(rbind, bias_list), 2, median, na.rm = TRUE) - psi_mini
  var_single = var_compute_m(bias_list, var_list)
  
  # Save all these results
  tmp = data.frame(N=N,
                   bias_wald = bias_single[1],
                   bias_if = bias_single[2],
                   bias_if_fw = bias_single[3],
                   bias_if_fw2 = bias_single[4],
                   bias2_wald = bias_single2[1],
                   bias2_if = bias_single2[2],
                   bias2_if_fw = bias_single2[3],
                   bias2_if_fw2 = bias_single2[4],
                   var_wald = var_single[1],
                   var_if = var_single[2],
                   var_if_fw = var_single[3],
                   var_if_fw2 = var_single[4]
  )
  return(tmp)
}


n_vec = c(300,600,1200,2400,3600)
SL_choice = c(1,2,3,4,5,6,7)

i = 1
p1= paste0("./R_157_300/",i,".Rda")
p2= paste0("./R_157_600/",i,".Rda")
p3= paste0("./R_157_1200/",i,".Rda")
p4= paste0("./R_157_2400/",i,".Rda")
p5= paste0("./R_157_3600/",i,".Rda")
r1 = single_run_simulation(df=df_pool,N=n_vec[1],CV=T,k=3,choice=SL_choice)
save(r1,file=p1)
r2 = single_run_simulation(df=df_pool,N=n_vec[2],CV=T,k=4,choice=SL_choice)
save(r2,file=p2)
r3 = single_run_simulation(df=df_pool,N=n_vec[3],CV=T,k=5,choice=SL_choice)
save(r3,file=p3)
r4 = single_run_simulation(df=df_pool,N=n_vec[4],CV=T,k=6,choice=SL_choice)
save(r4,file=p4)
r5 = single_run_simulation(df=df_pool,N=n_vec[5],CV=T,k=7,choice=SL_choice)
save(r5,file=p5)