source("QIV_utils.R")
N = 3000
current_df = single_gen_QIV(N)
y = current_df$y
a = current_df$a
z = current_df$z
x1 = current_df$x1
x2 = current_df$x2
nc=3
s <- rep(1:nc, each=N/nc)

y1 = current_df$y1
y0 = current_df$y0

Data = data.frame("x1"=x1,"x2"=x2,"a"=a,"z"=z,"y"=y,s)
psi_mini = mean((y1-y0)[a==1])

single_run_simulation_QIV = function(df,N,nc=3, fold_num = 1, CV,seed = NULL, k=4, bs= 'bs',choice=c(1,5,7), raw = F){
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
  
  psi_mini = mean((y1 - y0)[a==1])
  
  res_list = rep(list(NA), fold_num)
  bias_list = rep(list(NA), fold_num)
  var_list = rep(list(NA), fold_num)
  
  for(j in 1:fold_num){
    data_input <- data_input[sample(N), ]
    data_input$s = s
    result_single = three_fold_eval_QIV(data_input, CV=CV, k=k, bs = bs, choice=choice)
    res_list[[j]] = result_single
    bias_list[[j]] = as.data.frame(lapply(result_single, function(x) mean(x)))
    var_list[[j]] = as.data.frame(lapply(result_single, function(x) var_compute(x,data_input$a,data_input$s,nc)))
  }
  
  
  
  bias_single2 = apply(do.call(rbind, bias_list), 2, median, na.rm = TRUE) - psi_mini
  var_single = var_compute_m(bias_list, var_list)
  
  # Save all these results
  tmp = data.frame(N=N,
                   bias2_wald = bias_single2[1],
                   bias2_if = bias_single2[2],
                   var_wald = var_single[1],
                   var_if = var_single[2])
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
r1 = single_run_simulation_QIV(df=df_pool,N=n_vec[1],CV=T,k=3,choice=SL_choice)