#include <Rcpp.h>
//#include <Rcpp/vector.h>
using namespace Rcpp;
//Let batch enrollment at the beginning
// [[Rcpp::export]]
NumericVector O_E_CUSUM_hval(int nloop,int N3,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double p,double yr_int=1,double start_yr=0){
/*xmding:
  nloop:   number of loops
  N3:      sample size
  theta1:  alternative hypothesis (logRA1,logRA2), could be any length, as long as consisting with number of theta0
  theta0:  null hypothesis (logR0,-logR0)
  mu:      true relative risk logR0
  tau:     number of follow-up years
  yr_er:   the control population??s yearly motality rate ER=(total new events in one-year follow-up)/(total patient-years)
  p:       type I error
  yr_int:  follow-up years of interest (usually the whole period)
  start_yr:baseline year (usually year 0)
*/
  if(yr_int>tau) warning("Error: your yr_int>tau, this may be invalid!");
  int start_day=trunc(start_yr*365)+1;                         ///xmding: first day of interest
  NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;///xmding: k1 in 20190508 notes
  IntegerVector sign_c1=sign(c1); ///xmding: c1*theta1>0 (theta0=0 case)
  NumericVector abs_c1=abs(c1);
  NumericVector h_sel(c1.length()); 
  double yr_size=N3/tau;                              //////xmding: seems yr_size is unused
  double h_temp;
  int sN=floor(nloop*p);                              ///xmding: pivot of the quantile; should use floor to control type I error
  //////xmding: would it be better to add a warning in case sN=0? Though it usually would not happen.
  NumericMatrix h_quantiles(sN,c1.length());          ///xmding: fill with 0 by default
  double gamma_subject=-log(1-yr_er)/365;             ///xmding: \lambda_i(t); constant over at risk period
  int ndays=trunc(tau*365);                           ///xmding:trunc=floor when positive
  IntegerVector time_list=seq_len(ndays-start_day+1); ///xmding: length of intereting days
  time_list.push_front(0);                            ///xmding: add a zero to the front
  for(int loop=1;loop<=nloop; loop++){                //////xmding: consider moving variable-generating syntax out of this for loop
    //Rcout<<loop<<std::endl;                         //////xmding: annotated the syntax for print out loop
    //srand (loop);
    NumericVector enrl_t(N3);                         ///xmding: transplant time S_i
    enrl_t.fill(0);                                   //////xmding: E, O, and M need to fill with 0 too?
    //Rcout<<N3<<std::endl;
    NumericVector pre_time=trunc(rexp(N3,(gamma_subject)*exp(mu))); ///xmding: rexp generates random numbers from exponential distribution
    ///xmding: simulate failure time after S_i: X_i
    LogicalVector delta_list=((pre_time<=365*yr_int) & ((enrl_t+pre_time)<(tau*365)));
    ///xmding: qualified failure: failure happened within years of interest after transplant,
    ///xmding: and failure happened within follow up time
    NumericVector time=trunc(pmin(pre_time,pmin(365*yr_int,tau*365-enrl_t)));
    ///longhao: min(X_i-S_i,1,C_i-S_i)
    ///xmding: the shortest time among failure after transplant, years of interest, and follow-up period after transplant
    ///xmding: if delta_list=1, then time=pre_time
    NumericVector cho_time=enrl_t+time;
    ///xmding: the earliest time among failure, out of risk, stop follow up
    NumericMatrix E_t_pre(N3,ndays-start_day+1);
    NumericMatrix O_t_pre(N3,ndays-start_day+1);
    NumericVector M_t_pre(ndays-start_day+1);
    NumericVector M_t_pre2(ndays-start_day+1);
    for (int i=0;i<N3;i++){
      if(cho_time[i]>=start_day){
      ///xmding: failure, out of risk, stop follow up, all happened after start_day of interest
        //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
        int true_in=(start_day>enrl_t[i])?start_day:enrl_t[i];//enrl_t[i]
        ///xmding: true start day is the later of start_day of interest and the transplant day
        //int true_in=enrl_t[i];
        for(int j=0;j<ndays-start_day+1;j++){                       ///xmding: i for patient, j for time
          if(j>=cho_time[i]-start_day && delta_list[i]){
            O_t_pre(i,j)=1;                                         ///xmding: observed qualified failure; O=0 otherwise
          }
          if(j>=cho_time[i]-start_day) {///xmding: Lambda=lambda*t before failure and stays a constant after failure? simulation scenario
            E_t_pre(i,j)=gamma_subject*(cho_time[i]-true_in);// start_day<=cho_time && enrl_time< cho_Time
            ///xmding: (cho-start)-(true-start)
          } else if(j>=true_in-start_day && j<=cho_time[i]-start_day){
            E_t_pre(i,j)=gamma_subject*(j-(true_in-start_day));
          }
        }
        
      }
    }
    NumericVector O_t(ndays-start_day+1);
    NumericVector E_t (ndays-start_day+1);
    NumericVector O_E_t (ndays-start_day+1);
    for(int j=0; j<ndays-start_day+1; j++){
      O_t[j]=sum(O_t_pre.column(j));
      E_t[j]=sum(E_t_pre.column(j));
    }
    O_E_t=O_t-E_t;
    //print(O_t);
    //print(E_t);
    for(int k=0;k<c1.length();k++){
      M_t_pre=sign_c1[k]*O_E_t-abs_c1[k]*E_t;        //////xmding: sign_c*(O_E-c*E); no need to have abs_c;
      M_t_pre.push_front(0);
      M_t_pre2 = NumericVector(cummin(M_t_pre));     ///xmding: inf(C-kA)
      h_temp=max(M_t_pre-M_t_pre2);
      //Rcout<<h_temp<<std::endl;
      bool mainloop=false;
      for(int sn=0; sn<sN && mainloop==false; sn++){ ///xmding: order the largest sN h_temps decreasingly
        if(h_temp>h_quantiles(sn,k)){mainloop=true;
          for(int ssn=sN-1;ssn>sn;ssn--){h_quantiles(ssn,k)=h_quantiles(ssn-1,k);}
          h_quantiles(sn,k)=h_temp;
        }
      }
      // Rcout<<h_temp<<std::endl;
      // print(h_quantiles);
    }
    
  }
  //print(h_quantiles);
  return( h_quantiles.row(sN-1)); ///xmding: row(sN-1) is the sN_th row
  // return(0);
}

/***R
O_E_CUSUM_hval_resample=function(nloop,N_list,delta, enrl_t, cho_time, xbeta,Lambda0,theta1,theta0,p,tau=4,merge=4){
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

CUSUM_data_gen=function(mu,yr_size,tau,yr_er,yr_int=1,start_yr=1,seed=1,beta=0,change_yr=F,change_rate=0){
  set.seed(seed)
  size=trunc(yr_size*tau)
  gamma_subject=-log(1-yr_er)/365 #longhao: hazard rate
  ndays=trunc(tau*365)
  start_day=trunc(start_yr*365)+1
  time_list=0:(ndays-start_day+1)
  Lambda0=(0:ndays)*gamma_subject #longhao: cumulative hazard
  enrl_gp=rexp(size*2, yr_size)
  enrl_t=cumsum(enrl_gp)
  enrl_t=floor((enrl_t[enrl_t<tau])*365)
  N3=length(enrl_t)
  pre_time=trunc(rexp(N3,(gamma_subject)*exp(mu)))
  if(change_yr!=F){
    ind_change <- which((pre_time+enrl_t)>(change_yr*365))
    pre_time[ind_change]<- ceiling((pre_time[ind_change]+enrl_t[ind_change]-pmax(change_yr*365,enrl_t)[ind_change])*exp(mu)/exp(change_rate)+pmax(change_yr*365,enrl_t)[ind_change]-enrl_t[ind_change])
  }
  delta_list=((pre_time<=365*yr_int) & ((enrl_t+pre_time)<(tau*365)))
    time=trunc(pmin(pre_time,pmin(365*yr_int,tau*365-enrl_t)));
  cho_time=enrl_t+time;
  x=rnorm(N3)
    xbeta=x*beta
    return(list(time_list=time_list,delta_list=delta_list,cho_time=cho_time,enrl_t=enrl_t,xbeta=xbeta,Lambda0=Lambda0,total_size=N3))
    
}

*/

