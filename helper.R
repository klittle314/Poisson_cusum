
set.seed(1234)
##########################################################################################
# function to calculate test statistic
S_update <- function(S_current, k_value, current_count){
  
  S <- max(0,current_count - k_value + S_current, na.rm=TRUE)
  
  return(S)
}

#########################################################################################
#function to make nrep replications of Poisson series of length n, calculate cusum statistic using k_test
# N_day is the days of interest--e.g. what is the cusum statistic at 5 days, then N_day = 5
make_series <- function(n,nrep,k_test,Po_mean,N_day) {
  df_list <- list()
  
  df1 <- data.frame(matrix(ncol=4,nrow=0))
  
  df2 <- data.frame(matrix(ncol=3,nrow=0))
  
  names(df1) <- c("j_index","i_index","Poisson_values","cusum")
  
  names(df2) <- c("j_index","N_day_cusum", "max_N_day_cusum")
  
  for(j in 1:nrep) {
    Poisson_values <- rpois(n,Po_mean)
    
    cusum <- vector("double", n)
    
    i_index <-  vector("integer", n)
    
    S_current <- 0
    
    j_index <- j
    
    for(i in seq_along(Poisson_values)) {
      cusum[i] <- S_update(S_current,k_value = k_test,current_count = Poisson_values[i])
      
      S_current <- cusum[i]
      
      i_index[i] <- i
    }
    
    N_day_cusum <- cusum[N_day]
    
    max_N_day_cusum <- max(cusum[1:N_day], na.rm=TRUE)
    
    df1 <- rbind.data.frame(df1,data.frame(j_index,i_index,Poisson_values,cusum))
    
    df2 <- rbind.data.frame(df2,data.frame(j_index,N_day_cusum,max_N_day_cusum))
    
  }
  
  df_list$df1 <- df1
  
  df_list$df2 <- df2
  
  return(df_list)
}

###################################################################
#function to find Run Length of a series produced by the function make_series:   first value in vector that exceeds a given value
# 
find_RL <- function(j,data,h){
  
  df <- data %>% filter(j_index == j)
  #if there is no cusum value greater than h, function call returns Inf.  We convert to NA for further analysis.
  RL_test <- min(which(df$cusum >= h))
  
  RL <- ifelse(is.infinite(RL_test),NA,RL_test)
  
  return(RL)
}

################################################################################
#function to check h, k in Lucas 1985:  produce a vector of Run Lengths, starting with nrep_use repetitions of vectors of 
#length n_use, a random realization of Poission with mean = Po_mean_use.   The cusum statistic needs the value k_test.
#The function find_RL looks at each repetition of Poisson values and finds the first time the value of the cusum is at least h

check_hk_RL <- function(n_use,nrep_use,h_use,k_use,Po_mean_use) {
  
  df_list <- make_series(n=n_use,nrep=nrep_use,k_test=k_use,Po_mean=Po_mean_use, N_day = 7)
  
  index_RL <- c(1:max(df_list$df1$j_index))
  
  RL_out <- unlist(lapply(index_RL,find_RL,data=df_list$df1,h=h_use))
  
  return(RL_out)
}

##################h grid test set up function:  pair with lapply to get a summary table in function make_table_RL
RL_summary <- function(h_use0,n_use0,nrep_use0,k_use0,Po_mean_use0) {
  
  df_list <- make_series(n=n_use0,nrep=nrep_use0,k_test=k_use0,Po_mean=Po_mean_use0, N_day = 7)
  
  index_RL <- c(1:max(df_list$df1$j_index)) #or just use nrep?  should save a computation
  
  RL_vector <- unlist(lapply(index_RL,find_RL,data=df_list$df1,h=h_use0))
  
  #for large values of h, the simulated series may never have a cusum statistic that reaches or crosses the h boundary.  In find_RL, Inf is converted to NA.
  #here we detect NAs and insert NA in the named vector resulting from summary function, to align summaries in the make_table_RL function call to rbind.
  if(!anyNA(RL_vector)) {
    RL_summary <- c(summary(RL_vector),"NA's"=0)
  } else RL_summary <- summary(RL_vector)
  
  summary_out <- c(RL_summary,quantile(RL_vector,probs = c(.1,.5,1,2.5,5,10,50)/100, na.rm=TRUE))
  
  return(summary_out)
}

##################function to make a summary table of RL distribution info for given values of the parameters
#could use a grid for h based on seq function or a vector, e.g. to reproduce the levels of h in Lucan table 2.
make_table_RL <- function(h_vals_use,
                          n_use,
                          nrep_use,
                          k_use,
                          Po_mean_use){
  
                        
  
  h_grid_out <- lapply(h_vals_use,
                       RL_summary,
                       n_use0 = n_use, 
                       nrep_use0 = nrep_use,    
                       k_use0 = k_use, 
                       Po_mean_use0 = Po_mean_use)
  
  summary_table <- cbind.data.frame(h_vals_use,do.call(rbind,h_grid_out))
  
  names(summary_table)[1]  <- "h"
  
  return(summary_table)
}