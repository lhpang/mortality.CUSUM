O_E.limit=function(nloop,N_list,delta, enrl_t, cho_time, xbeta,Lambda0,theta1,theta0,p,tau=4,merge=4){
  h_val=NULL
  c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1
  sign_c1=sign(c1)
  abs_c1=abs(c1)
  N=length(xbeta)
  n_days=ceiling(366*(tau/merge))
  # if(n_days>length(Lambda0)) warning("Error:n_days>length(Lambda0).")
  E_t_pre=O_t_pre=matrix(0,nrow=N,ncol=n_days)
  for (i in 1:N){
    if(delta[i]){
      O_t_pre[i,(cho_time[i]+1):n_days]=1
    }
    E_t_pre[i,(cho_time[i]+1):n_days]=Lambda0[cho_time[i]-enrl_t[i]+1]*exp(xbeta[i])
    E_t_pre[i,(enrl_t[i]:cho_time[i]+1)]=Lambda0[1:(cho_time[i]-enrl_t[i]+1)]*exp(xbeta[i])
  }

  for(N3 in N_list){
    print(paste0("size=",N3))
    pb <- txtProgressBar(min = 1, max = nloop, style = 3)
    h_temp=matrix(NA,nrow=nloop,ncol=length(c1))
    for(loop in 1:nloop){
      setTxtProgressBar(pb, loop)
      E_t=O_t=0
      for(m in 1:merge){
        select=sample.int(N,size=N3,replace=T) #error! correctee on June 26 2018 size=round(N3/merge) should be changed to size=N3
        O_t=c(O_t,O_t[length(O_t)]+colSums(O_t_pre[select,]))
        E_t=c(E_t,E_t[length(E_t)]+colSums(E_t_pre[select,]))
      }

      O_E_t=O_t-E_t

      for(k in 1:length(c1)){
        M_t_pre=sign_c1[k]*O_E_t-abs_c1[k]*E_t
        M_t_pre2 = cummin(M_t_pre) #longhao:select the min among all subjects and time
        h_temp[loop,k]=max(M_t_pre-M_t_pre2)
      }
    }
    h_size=apply(h_temp,2,quantile,probs=1-p)
    print(h_size)
    h_val=rbind(h_val, h_size)
    close(pb)
  }
  return( data.frame(N_list,h_val))
}
