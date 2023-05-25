### Mortality CUSUM
Mortality_CUSUM<-function(date,ESCO_ord,data_path,facility_id_name,start_name,end_name,
                          dial_drd_name,xbeta_name,patient_id_name,start_date,end_date,
                          HPT1,HPT2,Lambda_list_prior,year_event_rate,pre,
                          ESCO_de_identified_path,y_axis_path,h_val_path,rho){

  ### Setting check
  if (!(HPT1==1.2&& HPT2==0.8) && is.null(h_val_path)){
    stop("If the hypothesis is not HPT1=1.2 and HPT2=0.8, please input the h_val_path accordingly,the default h_val_path is set to handle HPT1=1.2 and HPT2=0.8 hypothesis tests. ")
  }
  if (HPT1<=1||HPT2>=1||HPT2<=0){
    stop("HPT1 must be >1 and HPT2 must be in (0,1)")
  }
  if(rho>1||rho<0){
    stop("Invalid input of rho, which should be in range [0,1]")
  }
  if((year_event_rate>0.2||year_event_rate<0.1)&&is.null(h_val_path)){
    stop("Need to import your h value table through h_Val_path, because your year_event_rate is not in [0.1,0.2]")
  }
  #longhao: why?

  #parameters
  theta1=log(HPT1)
  theta2=log(HPT2)
  c1=(exp(theta1)-1)/theta1-1
  c2=(exp(theta2)-1)/theta2-1

  #whether and where to restart
  restart1=rho<1
  restart2=rho<1
  rho1=rho2=rho
  #what's the use of rho?

  # read in all the data
  message("Read in data...")
  data_all=read.csv(data_path)
  #longhao: prevent warnings
  #suppressWarnings(data_all <- fread(data_path, header=T, na.strings=c(".", "na", " "), sep=",", data.table=F))
  #longhao:y_axis is for future plotting
  #y_axis<-read.csv(y_axis_path)
  #bsline=read.csv(bsline_path)
  #Lambda_list_prior=-log(bsline[,bsline_name]) #cumulative hazards, -log(S(t))

  # Read in the hvalues
  message("Read in h values...")
  yer=round(year_event_rate,3)*100
  if(is.null(h_val_path)) {
    h_tab=read.csv("selected_h_test_1_2.csv")
    colnames(h_tab)<-c("N_all","flag_up_h","flag_low_h")
  }else {h_tab<-read.csv(h_val_path)
  colnames(h_tab)<-c("N_all","flag_up_h","flag_low_h")} # must have three columns
  h_val_size_max<-max(h_tab$N_all) #used to check whether the sample size is in the safe range: n<= h_val_size_max

  ### date and years
  startdate=as.Date(start_date, format = "%m/%d/%Y")
  enddate=as.Date(end_date, format = "%m/%d/%Y",origin=startdate)
  years=as.numeric(format(startdate,'%Y')):as.numeric(format(enddate,'%Y')) ## you need to change the follow-up years

  ###ordering
  if((!is.null(ESCO_de_identified_path))){
    #ESCO_de_identified<-read.csv(ESCO_de_identified_path)
    #longhao: find the ESCO_ord elements, for future plotting
    #ESCO_deorder<-match(ESCO_ord,ESCO_de_identified$ESCO)
  }
  if(!is.null(y_axis_path)){
    y_range_ord<-y_axis[match(ESCO_ord,y_axis$ESCO_ID),] #reorder the setting accordingly
  }

  #longhao: compute missingness
  n.incomplete=sum(!complete.cases(data_all))
  if(n.incomplete>0){
    warning(paste0("Becuase of missingness, delete ", n.incomplete," incomplete observations..."))
    data_all=na.omit(data_all) #delete observations that has missing varibles
  }
  ###########################################################
  ###Build the loop to calculate for each year
  if(length(years)>=2){
    startdatelist = startdate
    enddatelist = enddate
    for(i in 1:(length(years)-1)){
      startdatelist<-c(startdatelist,as.Date(paste("01/01",years[i+1],sep = "/"), format = "%m/%d/%Y",origin=startdate))
      enddatelist<-c(enddatelist,as.Date(paste("12/31",years[i],sep = "/"), format = "%m/%d/%Y",origin=startdate))
    }
    enddatelist = c(enddatelist[-1],enddate)
  }else if(length(years)==1){
    startdatelist=startdate
    enddatelist<-enddate
  }

  message("Analysis in progress...")
  size=NULL

  covariates_ls<-c(facility_id_name,start_name,end_name,dial_drd_name,xbeta_name,patient_id_name)
  for (year_ind in years){
    #year_ind=2015
    year_start=as.Date(paste0("01/01/",year_ind),format="%m/%d/%Y",origin=startdate)
    data=data_all[which(data_all$year==year_ind),]
    if(sum(data_all$year==year_ind)==0){
      break
    }
    # colnames(data)
    dataforuse=data[,covariates_ls]
    #start date
    X<-dataforuse[,start_name]<-as.Date(as.numeric(dataforuse[,start_name]),format = "%m/%d/%Y",origin="01/01/1970")
    delta=dataforuse[,dial_drd_name]
    time=as.Date(as.numeric(dataforuse[,end_name]),format="%m/%d/%Y",origin="01/01/1970")
    time_X=time-X+1 #time_X is the Follow-Up time in certain year

    #warn<-sum(abs(time_X-dataforuse[,dial_dar_name])*365.25>5) # use 5 days as an cutoff    ######xmding: why times 365.25, already in days
    #if(warn!=0) message(paste("Year",year_ind," Warning: different fu calculated from ", dial_dar," and start-end+1. ",
    #                          " The count for the problematic observations is  ", warn,". ",sep=""))
    facility=dataforuse[,facility_id_name] #ESCO_ID
    Xbeta=dataforuse[,xbeta_name]
    person_id<-dataforuse[,patient_id_name]
    N=nrow(dataforuse)
    N.uniq<-length(unique(person_id))
    if(N!=N.uniq) message(paste("Year ",year_ind,": different number of records and patients: ",N-N.uniq,sep=""))
    if(N>0&&sum(delta)/N<0.05) message(paste("Year ",year_ind," Small yearly event rate: ",sum(delta)/N,sep=""))
    delta_list=which(delta==1) #use 0.05 as a criteria a cutoff for a small event rate signal
    ####longhao: for each id, startdate and enddate is recorded year-wise
    num_start_dt<-as.numeric(startdatelist[which(year_ind==years)])
    num_end_dt<-as.numeric(enddatelist[which(year_ind==years)])
    ndays=num_end_dt- num_start_dt+1

    ########################################
    E_t_pre=matrix(0, nrow=N, ndays)
    O_t_pre=matrix(0, nrow=N, ndays)
    time_list=num_start_dt:num_end_dt
    for (i in 1:N){
      if (i %in% delta_list) {
        #longhao: qualifying failures
        O_t_pre[i,which(time_list>=time[i])]=1
      }
      #longhao:risk adjusted expected cumulative failures
      #longhao:the cumulative hazard stays the same after death?
      E_t_pre[i,which(time_list>=time[i])]=exp(Xbeta[i])*(Lambda_list_prior[time[i]-year_start+1]-Lambda_list_prior[X[i]-year_start+1])
      at.risk=which(time_list>=X[i]&time_list<time[i])
      E_t_pre[i,at.risk]=exp(Xbeta[i])*(Lambda_list_prior[time_list[at.risk]+1-as.numeric(year_start)]-Lambda_list_prior[X[i]-year_start+1])
    }

    #longhao: a function splits data into subsets according to ESCO_ord
    facility_sum=function(x){
      as.vector(sapply(split(x,factor(facility,levels=ESCO_ord)), sum)) # will be ordered by facility itself automatically
    }

    #longhao: O_t is a matrix with nrow=#facility, ncol=$nday
    O_t=apply(O_t_pre, 2, facility_sum)
    E_t=apply(E_t_pre, 2, facility_sum)
    O_E_t=O_t-E_t

    quartersize=function(X,time,person_id,facility,yr,ESCO=ESCO_ord){
      w1=which(X>as.Date(paste("3/31/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) |
                 time<as.Date(paste("1/1/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) )
      w2=which(X>as.Date(paste("6/30/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) |
                 time<as.Date(paste("4/1/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) )
      w3=which(X>as.Date(paste("9/30/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) |
                 time<as.Date(paste("7/1/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) )
      w4=which(X>as.Date(paste("12/31/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) |
                 time<as.Date(paste("10/1/",yr,sep=""), format = "%m/%d/%Y",origin=startdate) )
      if(length(w1)==0)  v1=sapply(split(person_id,factor(facility,levels=ESCO)),function(v)length(unique(v))) else v1=sapply(split(person_id[-w1],factor(facility[-w1],levels=ESCO)),function(v)length(unique(v)))
      #longhao: for each subset, compute the number of unique ESCO group of person_id for people in each quarter
      if(length(w2)==0)  v2=sapply(split(person_id,factor(facility,levels=ESCO)),function(v)length(unique(v))) else v2=sapply(split(person_id[-w2],factor(facility[-w2],levels=ESCO)),function(v)length(unique(v)))

      if(length(w3)==0)  v3=sapply(split(person_id,factor(facility,levels=ESCO)),function(v)length(unique(v))) else v3=sapply(split(person_id[-w3],factor(facility[-w3],levels=ESCO)),function(v)length(unique(v)))

      if(length(w4)==0)  v4=sapply(split(person_id,factor(facility,levels=ESCO)),function(v)length(unique(v))) else v4=sapply(split(person_id[-w4],factor(facility[-w4],levels=ESCO)),function(v)length(unique(v)))


      return(t(rbind(v1,v2,v3,v4)))
    }
    size_part1=quartersize(X,time,person_id,facility,year_ind) # quarterly size(#person_id) by facilities/ESCOs
    size_part2=sapply(split(X[!duplicated(person_id)],factor(facility[!duplicated(person_id)],levels=ESCO_ord)),length) # yearly size by facilities/ESCOS
    size_p<-cbind(size_part2,size_part1)
    colnames(size_p)=c(paste(year_ind,"_total",sep=""),paste(year_ind,"Q1",sep=""),paste(year_ind,"Q2",sep=""),paste(year_ind,"Q3",sep=""),paste(year_ind,"Q4",sep=""))
    size=cbind(size,size_p)
    ESCO=c("Date",ESCO_ord) #longhao: should be ESCO here?
    #longhao:for each day, the expected cumulative hazard for each strata
    E_t=cbind(ESCO,rbind(time_list,E_t))
    O_t=cbind(ESCO,rbind(time_list,O_t))
    O_E_t=cbind(ESCO,rbind(time_list,O_E_t))
    write.csv(E_t,paste("E_t_",year_ind,'_',date,".csv",sep=""),row.names = F)
    write.csv(O_t,paste("O_t_",year_ind,'_',date,".csv",sep=""),row.names = F)
    write.csv(O_E_t,paste("O_E_t_",year_ind,'_',date,".csv",sep=""),row.names = F)
  } #end of for loop in years
  #count the number of at risk patients for each quarters, combine with CUSUM_prep code
  size_max=max(size)
  if(size_max>h_val_size_max) {
    warning(paste0("
                   Your sample size exceeds the maximum size the current h values table is able to handle,
                   the largest size in your data for each facility is ",size_max,", but the current table is ",max(h_tab$N_all),".
                   Please contact lilywang@umich.edu to update the default form for you (include your largest size value shown above),
                   or conduct the simulations yourself using the simulation code and direct it through variable h_val_path"))
  }
  size[is.na(size)]<-0
  write.csv(size,paste("Inform_",date,".csv",sep=""))

  #############################################################################
  ###CUSUM prep part##########################################################
  #############################################################################

  ##########################################################################
  #read in the sized for each year
  inform=read.csv(paste("Inform_",date,".csv",sep=""),header=T)
  ##########################################################################
  #longhao: year size smaller than all the h_tab
  pick<-function(x,h_m){
    h1=h_m$flag_up_h[which(x<=h_m$N_all)][1]
    h2=h_m$flag_low_h[which(x<=h_m$N_all)][1]
    return (c(h1,h2))
  }
  ##########################################################################
  #Read in the O and E datasets we calculated and merge them
  H_values=NULL
  for (year_ind in years){
    if(sum(data_all$year==year_ind)==0){
      break
    }
    N=inform[,paste("X",year_ind,"_total",sep="")] #yearly size: update the hvalues by yearly sizes
    h=sapply(N,pick,h_tab)
    rownames(h)=c(paste("upper_",year_ind,sep=""),paste("lower_",year_ind,sep=""))
    h1=h[1,]
    h2=h[2,]
    H_values=cbind(H_values,t(h))
    E_t=read.csv(file=paste("E_t_",year_ind,'_',date,".csv",sep=""))
    O_t=read.csv(file=paste("O_t_",year_ind,'_',date,".csv",sep=""))
    O_E_t=read.csv(file=paste("O_E_t_",year_ind,'_',date,".csv",sep=""))
    d2=dim(E_t)
    time=as.vector(unlist(E_t[1,2:d2[2]])) #longhao: each day for a year
    ESCO=as.vector(unlist(E_t[2:d2[1],1]))
    E_t=E_t[2:d2[1],2:d2[2]]
    O_t=O_t[2:d2[1],2:d2[2]]
    O_E_t=O_E_t[2:d2[1],2:d2[2]]
    h1_m=matrix(rep(h1,d2[2]-1),nrow=d2[1]-1,ncol=d2[2]-1)
    h2_m=matrix(rep(h2,d2[2]-1),nrow=d2[1]-1,ncol=d2[2]-1)
    if (year_ind==years[1]) {
      Comb_E_t=E_t
      Comb_O_t=O_t
      Comb_O_E_t=O_E_t
      H1=h1_m
      H2=h2_m
      TIME=time
      i=0
    }
    else{
      i=i+1
      d1=dim(Comb_E_t)
      TIME=c(TIME,time)
      E_t=Comb_E_t[,d1[2]]+E_t #longhao: next year add the cumulative hazard of the last day of last year
      O_t=Comb_O_t[,d1[2]]+O_t
      O_E_t=Comb_O_E_t[,d1[2]]+O_E_t
      Comb_E_t=cbind(Comb_E_t,E_t)
      Comb_O_t=cbind(Comb_O_t,O_t)
      Comb_O_E_t=cbind(Comb_O_E_t,O_E_t) #use for M_t_restart's middle curve
      H1=cbind(H1,h1_m)
      H2=cbind(H2,h2_m)
    }
  }
  ############################################################################
  #Redefine M_t_restart to have the restarting function
  FU=1:ncol(Comb_E_t) #longhao: FU refers to the first day of first year till last day of last year
  N=ceiling(as.numeric(enddate-startdate)/91) #number of quarters  ######xmding: ceiling instead of round? 91.3125?
  ndays=length(FU)
  nesco=length(ESCO)

  M1_pre=Comb_O_E_t-c1*Comb_E_t #use for M_t_restart
  M2_pre=-Comb_O_E_t+c2*Comb_E_t

  M1_pt1=NULL #vectors
  M1_pt2=NULL
  M2_pt1=NULL
  M2_pt2=NULL
  cross1=rep(0,nesco)
  cross1_t=vector("list",nesco)
  cross2=rep(0,nesco)
  cross2_t=vector("list",nesco)

  M1_restart=matrix(nrow=nesco,ncol=ndays)
  M2_restart=matrix(nrow=nesco,ncol=ndays)
  #O_E_restart=matrix(nrow=nesco,ncol=ndays)

  M1_min=M1_pre[,1]
  M2_min=M2_pre[,1]
  start1=rep(1,nesco) #record the restart time for hypothesis 1
  start2=rep(1,nesco) #record the restart time for hypothesis 2

  for(i in 1:nesco){
    for(t in 1:ndays){
      M1_min[i]=min(M1_min[i],M1_pre[i,t])
      M2_min[i]=min(M2_min[i],M2_pre[i,t])
      M1_pt1[i]=M1_min[i]-M1_pre[i,t]+H1[i,t]          ###xmding: both M1 deduce M1_pre[i,start1[i]]
      M2_pt1[i]=M2_min[i]-M2_pre[i,t]+H2[i,t]
      if(restart1==T && cross1[i]>0) M1_pt2[i]=(1-rho1)*H1[i,t]-(M1_pre[i,t]-M1_pre[i,start1[i]]) #incoporate restart mechanism in the difference
      ###xmding: why -M1_pre[i,start1[i]]? restart ignore past information
      if(restart2==T && cross2[i]>0) M2_pt2[i]=(1-rho2)*H2[i,t]-(M2_pre[i,t]-M2_pre[i,start2[i]])
      M1_restart[i,t]=min(M1_pt1[i],M1_pt2[i],na.rm=T)
      M2_restart[i,t]=min(M2_pt1[i],M2_pt2[i],na.rm=T)
      #O_E_restart[i,t]=Comb_O_E_t[i,t]-Comb_O_E_t[i,start[i]] with shared or min start time
      if((M1_restart[i,t]<0)||(M2_restart[i,t]<0)){
        if(M1_restart[i,t]<0 && restart1==T) {
          start1[i]=t
          cross1[i]=cross1[i]+1 #longhao: one more cross if restart and signal, if signal but not restart, don't record
          cross1_t[[i]]=c(cross1_t[[i]],t)
        }else if(M1_restart[i,t]<0 && restart1==F && sum(cross1[i],cross2[i])==0){
          cross1[i]=1                                  ######xmding: if restart=F, only record the first cross? No
          cross1_t[[i]]=t
        }
        if(M2_restart[i,t]<0 && restart2==T) {
          start2[i]=t #longhao: reset the start date to t
          cross2[i]=cross2[i]+1
          cross2_t[[i]]=c(cross2_t[[i]],t)
        }else if(M2_restart[i,t]<0 && restart2==F && sum(cross1[i],cross2[i])==0){
          cross2[i]=1
          cross2_t[[i]]=t
        }

        if(restart1==T||restart2==T){
          M1_min[i]=M1_pre[i,start1[i]]                ###xmding: restart, remove information before cross
          M2_min[i]=M2_pre[i,start2[i]]
        }
      }
    }
  }


  ######################################################
  #longhao:two boundaries
  O_E_upper=Comb_O_E_t[,FU]+M1_restart
  O_E_lower=Comb_O_E_t[,FU]-M2_restart
  ###################
  #longhao: indicate which are restarted
  indicate_lab=matrix(0,dim(M1_restart)[1],dim(M1_restart)[2])
  for(i in 1:nrow(M1_restart)){
    for(j in 1:ncol(M1_restart)){
      if(M1_restart[i,j]<0) {indicate_lab[i,j]=1}
      else if(M2_restart[i,j]<0) {indicate_lab[i,j]=-1}
    }
  }


  ###################################################################################
  #Plot the data
  message("Producing plots...")

  N=round(as.numeric(enddate-startdate)/91) #number of quarters ######xmding: ceiling; 91.3125
  lab.dt<-vector("list",N)
  lab<-vector("character",N)
  i=1
  for(year_ind in years){
    if(sum(data_all$year==year_ind)==0){
      break
    }
    q1<-as.Date(paste("3/31/",year_ind,sep=""), format = "%m/%d/%Y",origin=startdate)
    q2<-as.Date(paste("6/30/",year_ind,sep=""), format = "%m/%d/%Y",origin=startdate)
    q3<-as.Date(paste("9/30/",year_ind,sep=""), format = "%m/%d/%Y",origin=startdate)
    q4<-as.Date(paste("12/31/",year_ind,sep=""), format = "%m/%d/%Y",origin=startdate)

    if(q1>startdate && q1<=enddate) {lab.dt[[i]]<-q1;lab[i]=paste(year_ind," Q1",sep="");i=i+1}
    if(q2>startdate && q2<=enddate) {lab.dt[[i]]<-q2;lab[i]=paste(year_ind," Q2",sep="");i=i+1}
    if(q3>startdate && q3<=enddate) {lab.dt[[i]]<-q3;lab[i]=paste(year_ind," Q3",sep="");i=i+1}
    if(q4>startdate && q4<=enddate) {lab.dt[[i]]<-q4;lab[i]=paste(year_ind," Q4",sep="");i=i+1}

  }


  par(mfrow = c(1, 1))
  O_E_upper_plot=O_E_upper
  O_E_lower_plot=O_E_lower
  FU_plot=FU
  for (i in 1:nesco) {
    Ylimit=c(min(c(min(O_E_lower[i,]),min(O_E_upper[i,]),min(Comb_O_E_t[i,])))-8, max(c(max(O_E_lower[i,]),max(O_E_upper[i,]),max(Comb_O_E_t[i,])))+15)
    need1=(restart1==F && cross1[i]>0)
    need2=(restart2==F && cross2[i]>0)
    need3=(restart1==T && cross1[i]>0)
    need4=(restart2==T && cross2[i]>0)
    if(need1) O_E_upper_plot[i,FU_plot>cross1_t[[i]][1]]<-NA ######xmding: not plot time after cross if restart=F
    if(need2) O_E_lower_plot[i,FU_plot>cross2_t[[i]][1]]<-NA

    png(paste(pre,ESCO[i],"_",date,".png",sep=""),height=800,width=800,res=150)
    if(is.null(y_axis_path)){
      plot(TIME, Comb_O_E_t[i,FU_plot], pch=1, lty=1, type="l",xaxt='n', ylim=Ylimit,main=paste("Observed-Expected Mortality CUSUM: ",ESCO[i],sep=""),xlab="",ylab='Excess Deaths (O-E)')#,sub="Weak Hyp")
    }else{
      plot(TIME, Comb_O_E_t[i,FU_plot], pch=1, lty=1, type="l",xaxt='n', ylim=unlist(y_range_ord[i,2:3]),main=paste("Observed-Expected Mortality CUSUM: ",ESCO[i],sep=""),xlab="",ylab='Excess Deaths (O-E)')#,sub="Weak Hyp")
    }
    lines(TIME,O_E_upper_plot[i,FU_plot], col='red', lty=2,cex=0.8)
    lines(TIME,O_E_lower_plot[i,FU_plot], col='blue', lty=2,cex=0.8)
    mtext("Continuous O-E")
    axis(1, at = as.numeric(lab.dt), las=2,label=lab,cex.lab=0.8,cex.axis=0.8)
    #legend("topleft",legend=c("upper limit","O-E","lower limit"),text.font=2,bty="n",lty=c(2,1,3),col=c("red","black","blue"),cex=0.7,horiz = T)
    if(need1||need3){
      abline(v=cross1_t[[i]]+as.numeric(startdate),col='red',lty=1) #longhao:signal points
      #arrows(x0=(cross1_t[[i]]-60),y0=as.vector(unlist(O_E_upper_plot[i,cross1_t[[i]]])),x1=cross1_t[[i]],y1=as.vector(unlist(O_E_upper_plot[i,cross1_t[[i]]])),col='red',lwd=1.5,length=0.1)
    }
    if(need2||need4){
      abline(v=cross2_t[[i]]+as.numeric(startdate),col='blue',lty=1)
      #arrows(x0=(cross2_t[[i]]-60),y0=as.vector(unlist(O_E_lower_plot[i,cross2_t[[i]]])),x1=cross2_t[[i]],y1=as.vector(unlist(O_E_lower_plot[i,cross2_t[[i]]])),col='blue',lwd=1.5,length=0.1)
    }
    dev.off()
    if((!is.null(ESCO_de_identified_path))){
      png(paste(pre,ESCO_deorder[i],"_",date,".png",sep=""),height=800,width=800,res=150)
      if(is.null(y_axis_path)){
        plot(TIME, Comb_O_E_t[i,FU_plot], pch=1, lty=1, type="l",xaxt='n', ylim=Ylimit,main=paste("Observed-Expected Mortality CUSUM: ESCO ",ESCO_deorder[i],sep=""),xlab="",ylab='Excess Deaths (O-E)')#,sub="Weak Hyp")
      }else{
        plot(TIME, Comb_O_E_t[i,FU_plot], pch=1, lty=1, type="l",xaxt='n', ylim=unlist(y_range_ord[i,2:3]),main=paste("Observed-Expected Mortality CUSUM: ",ESCO_deorder[i],sep=""),xlab="",ylab='Excess Deaths (O-E)')#,sub="Weak Hyp")
      }######xmding: difference in main
      lines(TIME,O_E_upper_plot[i,FU_plot], col='red', lty=2,cex=0.8)
      lines(TIME,O_E_lower_plot[i,FU_plot], col='blue', lty=2,cex=0.8)
      mtext("Continuous O-E")
      axis(1, at = as.numeric(lab.dt), las=2,label=lab,cex.lab=0.8,cex.axis=0.8)
      #legend("topleft",legend=c("upper limit","O-E","lower limit"),text.font=2,bty="n",lty=c(2,1,3),col=c("red","black","blue"),cex=0.7,horiz = T)
      if(need1||need3){
        abline(v=cross1_t[[i]]+as.numeric(startdate),col='red',lty=1)
        #arrows(x0=(cross1_t[[i]]-60),y0=as.vector(unlist(O_E_upper_plot[i,cross1_t[[i]]])),x1=cross1_t[[i]],y1=as.vector(unlist(O_E_upper_plot[i,cross1_t[[i]]])),col='red',lwd=1.5,length=0.1)
      }
      if(need2||need4){
        abline(v=cross2_t[[i]]+as.numeric(startdate),col='blue',lty=1)
        #arrows(x0=(cross2_t[[i]]-60),y0=as.vector(unlist(O_E_lower_plot[i,cross2_t[[i]]])),x1=cross2_t[[i]],y1=as.vector(unlist(O_E_lower_plot[i,cross2_t[[i]]])),col='blue',lwd=1.5,length=0.1)
      }
      dev.off()
    }
  }

  #############################################################
  ##################  Summary of the H values #################
  #############################################################
  #H_ind<-c(seq(1,2*length(years)-1,by=2),seq(2,2*length(years),by=2))
  #H_summary<-H_values[,H_ind]
  rownames(H_values)<-ESCO
  write.csv(H_values,paste(pre,"Hsummary_",date,".csv"))
  ##############################################################
  ### Summary of the Data's for plot and tabulization ##########
  ##############################################################
  quarter_ind=as.numeric(lab.dt)-as.numeric(startdate)+1
  N=round(as.numeric(enddate-startdate)/91) #number of quarters          ######xmding: ceiling; 91.3125
  Qmed<-sapply(lab,function(a) gsub(" ", "", a) ) #longhao: remove space for every element in lab
  names(Qmed)<-NULL
  Qind=match(sapply(Qmed,function(x) paste("X",x,sep="")),names(inform)) #longhao: why pasting "X" may be related to names(inform)
  Hit<-function(v1,v2,cp){
    v<-rep('N',length(cp))
    low_Q<-up_Q<-both_Q<-NULL
    up_Q=unlist(sapply(v1,function(x) which(x<=cp)[1]))
    low_Q=unlist(sapply(v2,function(x) which(x<=cp)[1]))
    v[intersect(up_Q,low_Q)]<-"UL"
    v[setdiff(up_Q,low_Q)]<-'U'
    v[setdiff(low_Q,up_Q)]<-'L'
    return(v)
  } #longhao: indicate the hitting situation for each quarter



  summary_table_nona=summary_table=NULL
  table=NULL
  for(i in 1:nesco){
    table=data.frame(ESCO=ESCO[i],Quarter=lab,sample_size=as.vector(t(unname(inform[i,Qind]))),Observe=as.vector(t(unname(Comb_O_t[i,quarter_ind]))),Expected=as.vector(t(unname(Comb_E_t[i,quarter_ind])))
                     ,Observed_Expected=as.vector(t(unname(Comb_O_E_t[i,quarter_ind]))),Lower_limit=as.vector(t(unname(O_E_lower_plot[i,quarter_ind]))),Upper_limit=as.vector(t(unname(O_E_upper_plot[i,quarter_ind])))
                     ,Crossing_Indicator=Hit(cross1_t[[i]],cross2_t[[i]],quarter_ind))
    summary_table=rbind(summary_table,table,rep(NA,ncol(table))) #longhao: with NA row as transition between each ESCO_ord
    summary_table_nona=rbind(summary_table_nona,table)
  }

  #summary_table[1:16,]

  write.csv(summary_table,paste(pre,"Summary_",date,".csv",sep=""))
  write.csv(summary_table_nona,paste(pre,"Summary_nona",date,".csv",sep=""))

  message("Analysis finished :-)")

}
#comment two axis
