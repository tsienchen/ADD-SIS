library(readr)
library(Icens)
library(interval)
library(MASS)
library(dplyr)
library(energy)


############################### Loading Data  ########################
############################### 
############################### 
############################### 
############################### 
dt <- read_csv("WordMatrix_Filter.csv")
REAL_TEST = subset(dt, select = -c(X1) )

l = as.numeric(REAL_TEST$l) 
r =  as.numeric(REAL_TEST$r) 
s = as.data.frame(subset(REAL_TEST, select = -c(l,r) ))

############################### 
############################### 
intervals <- initcomputeMLE(l,r)
max_tmp <- max(intervals$intmap[2,])
min_tmp <- min(intervals$intmap[1,])
if (is.finite(max_tmp)) {end <- max_tmp} 
if (is.infinite(max_tmp)) {end <-max(intervals$intmap[1,])+100 }
if (is.finite(min_tmp)) {start <- min_tmp} 
if (is.infinite(min_tmp)) {start <-min(intervals$intmap[2,])-100 }
X1<-c(as.integer(start/100):as.integer(end/100))
X2<-X1+1
basis<- data.frame(x1=100*X1,x2=100*X2)

int_split <- function(i,df){
  start <- as.numeric(df[i,][1])/100
  end <- as.numeric(df[i,][2])/100
  difs <- end-start-1
  x1 <- 100*(start+c(0:difs))
  x2 <- 100*(start+1+c(0:difs))
  x3 <- as.numeric(rep(df[i,][3]/(as.numeric(difs+1)),difs+1))
  return (data.frame(x1=x1,x2=x2,x3=x3))
  
}



screen_add <- function(no,s,l,r) {

  #word1 <- sapply(c(1:length(sub)[1]),function(i) grepl(dict[no],kk[i]))
  x1 <- s[,no]
  every = data.frame(ll=l,rr=r,x=x1)
  part1 = every[every$x==TRUE,c('ll','rr')]
  part2 = every[every$x==FALSE,c('ll','rr')]
  a <- initcomputeMLE(part1$ll,part1$rr)
  b <- initcomputeMLE(part2$ll,part2$rr)
  
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))
  if (length(a$pf)==1){
    if (max(df2$X2)>1) {kk2 = lapply(c(1:dim(df2)[1]),int_split,df=df2)
    pp2 = bind_rows(kk2)}
    if (max(df2$X2)==1) {ind_start <- (df2$X1/100)[1]
    ind_end<- (df2$X1/100)[2]
    pp2 <- data.frame(x1=100*c(ind_start:(ind_end-1)),
                      x2=100*c((ind_start+1):(ind_end)),
                      x3=rep(1/(end-start),as.numeric(ind_end-ind_start)))
    }
    tmp = merge(pp2,basis,all.y = TRUE)
    tmp[is.na(tmp$x3),c('x3')]<-0
    score_ = sum(100 * (cumsum(tmp$x3)))
  }
  else {
    kk1 = lapply(c(1:dim(df1)[1]),int_split,df=df1)
    pp1 = bind_rows(kk1)
    
    if (max(df2$X2)>1) {kk2 = lapply(c(1:dim(df2)[1]),int_split,df=df2)
    pp2 = bind_rows(kk2)}
    if (max(df2$X2)==1) {ind_start <- (df2$X1/100)[1]
    ind_end<- (df2$X1/100)[2]
    pp2 <- data.frame(x1=100*c(ind_start:(ind_end-1)),
                      x2=100*c((ind_start+1):(ind_end)),
                      x3=rep(1/(end-start),as.numeric(ind_end-ind_start)))
    }
    a1 = merge(pp1,basis,all.y = TRUE)
    a2 = merge(pp2,basis,all.y = TRUE)
    tmp <- merge(a1,a2,by='x1')
    tmp[is.na(tmp$x3.x),c('x3.x')]<-0
    tmp[is.na(tmp$x3.y),c('x3.y')]<-0
    score_ = sum(100 * abs(cumsum(tmp$x3.x)-cumsum(tmp$x3.y)))}
  
  return(score_)
  
  
}

#score <- sapply(c(1:as.numeric(dim(s)[2])),screen_add,s=s)
#res_ <- data.frame(screen=score)

#fitting_prepare <- function(res,sub){
#  totno <- dim(s)[2]+ 1
#  ks <- data.frame(rank=totno-rank(res$screen),index=c(1:dim(s)[2] ))
#  x <- as.matrix(s[,ks[ks$rank<sub+1,]$index])
#  data <- data.frame(l=l,r=r,sapply(c(1:dim(x)[2]),function(i) x[,i]))
#  return (data)
  
#}


#simulation_matrix<- function(res,sub,input_matrix){
#  totno <- dim(input_matrix)[2]+ 1
#  ks <- data.frame(rank=totno-rank(res$screen),index=c(1:dim(input_matrix)[2] ))
#  print (ks[ks$rank<sub+1,]$index)
#  x1 <- as.matrix(input_matrix[,ks[ks$rank<sub+1,]$index])
#  x2 <- as.matrix(input_matrix[,ks[ks$rank>=sub+1,]$index])
#  x <- cbind(x1,x2)
#  data <- data.frame(sapply(c(1:dim(x)[2]),function(i) x[,i]))
#  return (data)
  
#}



fitting_prepare <- function(res,index){
  #totno <- dim(s)[2]+ 1
  #ks <- data.frame(rank=totno-rank(res$screen),index=c(1:dim(s)[2] ))
  x <- as.matrix(s[,index])
  data <- data.frame(l=l,r=r,sapply(c(1:dim(x)[2]),function(i) x[,i]))
  return (data)
  
}


simulation_matrix<- function(res,index,input_matrix){
  #totno <- dim(input_matrix)[2]+ 1
  #ks <- data.frame(rank=totno-rank(res$screen),index=c(1:dim(input_matrix)[2] ))
  #print (ks[ks$rank<sub+1,]$index)
  x1 <- as.matrix(input_matrix[,index])
  x2 <- as.matrix(input_matrix[,c(1:677)[-index]])
  x <- cbind(x1,x2)
  data <- data.frame(sapply(c(1:dim(x)[2]),function(i) x[,i]))
  return (data)
  
}

model_fitting <- function(data,dis) {
  dt_tmp <- data[3:ncol(data)]
  if (dis == "weibull"){
    dt_tmp$l = data$l
    dt_tmp$r = data$r
    total <- dim(data)[2] - 2
    xkk = paste("X", 1:total, sep = "",collapse="+")
    model <- survreg(as.formula(paste("Surv(l,r,type='interval2')",
                                      xkk,sep="~")),data=dt_tmp,dist=dis)
    #scsc <- model$loglik
    #ress <- c(scsc[2],2*(scsc[2]-scsc[1]))
  }
  else {
    dt_tmp$l = log(data$l)
    dt_tmp$r = log(data$r)
    total <- dim(data)[2] - 2
    xkk = paste("X", 1:total, sep = "",collapse="+")
    model <- survreg(as.formula(paste("Surv(l,r,type='interval2')",
                                      xkk,sep="~")),data=dt_tmp,dist=dis)
    #scsc <- model$loglik
    #ress <- c(scsc[2],2*(scsc[2]-scsc[1]))}
  
  }
  return (model)
  }
  
screen_ks<- function(i,y,input)
{every = data.frame(cbind(y,input[,i]))
part1 = every[every$V2==0,c('y')]
part2 = every[every$V2==1,c('y')]
return (as.numeric(ks.test(part1,part2)$statistic))
}




screen_dc <- function(i,y,input)
{return (dcor(y,input[,i]))}




Fk<-function(X0,x) {
  Fk=c()
  for (i in 1:length(x))
  { Fk[i]=sum(X0<=x[i])/length(X0) }
  return(Fk)
}

Fkr<-function(X0,Y,yr,x) {
  Fkr=c()
  ind_yr=(Y==yr)
  for (i in 1:length(x))
  { Fkr[i]=sum((X0<=x[i])*ind_yr)/sum(ind_yr) }
  return(Fkr)
}

MV<-function(Xk,Y) {
  Fk0 <- Fk(Xk,Xk)
  Yr <- unique(Y)
  MVr <- c()
  for (r in 1:length(Yr)) {
    MVr[r] <- (sum(Y==Yr[r])/length(Y))*mean((Fkr(Xk,Y,Yr[r],Xk)-Fk0)^2)
  }
  MV <- sum(MVr)
  return(MV)
}

screen_mv <- function(i,input,y)
{return (abs(MV(Y=input[,i],Xk=y)))}






robustness_check <-function(no,effective_p ){
  print(no)
  set.seed(no)
  score <- sapply(c(1:as.numeric(dim(s)[2])),screen_add,s=s,l=l,r=r)
  res_ <- data.frame(screen=score)
  input_index <- c(1:677)[678 - rank(res_$screen)<1+effective_p]
  gr = model_fitting(fitting_prepare(res_,index=input_index),"gaussian")
  value_ground_tmp <- predict(gr,new_data=fitting_prepare(res_,index=input_index))
  value_ground <- as.numeric(exp(value_ground_tmp + 0.5 * rnorm(n=497,sd = gr$scale)))
  
  
  
  left = 100 * sample(c(rpois(n=248,lambda=20),rpois(n=249,lambda=10)))
  right = 100 * sample(c(rchisq(n=248,df=20) ,rchisq(n=249,df=40) ))
  observe_l =   1000 * floor((value_ground - left)/1000)
  observe_r = 500 * floor((value_ground + right)/500) + 500
  observe_mid = 0.5*observe_l + 0.5* observe_r
  observe_three = 2/3*observe_l + 1/3* observe_r
  observe_six = 1/3*observe_l + 2/3* observe_r
  
  
  #left = 400 * rpois(n=497,lambda=3)
  #right = 400 * rchisq(n=497,df=4)
  #observe_l =   100 * floor((value_ground - left)/100)
  #observe_r = 100 * floor((value_ground + right)/100)
  
  
  x_simu <- simulation_matrix(res_,input_index,s)
  
  
  
  score1 <- sapply(c(1:as.numeric(dim(x_simu)[2])),screen_add,s=x_simu,l=observe_l,r=observe_r)
  
  res1 <- data.frame(screen=score1)
  score2 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_ks,y=observe_mid,input=x_simu)
  res2 <- data.frame(screen=score2)
  score3 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_dc,y=observe_mid,input=x_simu)
  res3 <- data.frame(screen=score3)
  score4 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_mv,y=observe_mid,input=x_simu)
  res4 <- data.frame(screen=score4)
  
  
  
  score5 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_ks,y=observe_three,input=x_simu)
  res5 <- data.frame(screen=score5)
  score6 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_dc,y=observe_three,input=x_simu)
  res6 <- data.frame(screen=score6)
  score7 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_mv,y=observe_three,input=x_simu)
  res7 <- data.frame(screen=score7)
  score8 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_ks,y=observe_six,input=x_simu)
  res8 <- data.frame(screen=score8)
  score9 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_dc,y=observe_six,input=x_simu)
  res9 <- data.frame(screen=score9)
  score10 = sapply(c(1:as.numeric(dim(x_simu )[2])),screen_mv,y=observe_six,input=x_simu)
  res10 <- data.frame(screen=score10)
  
  res_all <- as.data.frame(cbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10))
  
  path = paste("../Intermediate/Dim",effective_p,"Trial",no,".csv",sep="")
  
  write.csv(res_all,path)
}



sapply(c(1:100),robustness_check,effective_p=30)

sapply(c(1:100),robustness_check,effective_p=60)

sapply(c(1:100),robustness_check,effective_p=120)
