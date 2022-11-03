
library(Icens)
library(interval)
library(MASS)
library(dplyr)
library(energy)
### basic functional definations for the discrete censored interval 
screen_npmle1 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))
  if (a$pf==1){
    
    data2 <-  transs(df2)
    tmp2 <-  merge(data2, basis,all.y = TRUE)
    sta <- cumsum(tmp2$X3)}
  else if  (b$pf==1){
    
    data1 <-  transs(df1)
    tmp1 <-  merge(data1, basis,all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  }
  else {data1 <-  transs(df1)
  data2 <-  transs(df2)
  tmp1 <- merge(data1, basis,all.y = TRUE)
  tmp2 <-  merge(data2, basis,all.y = TRUE)
  tmp <- merge(tmp1,tmp2,by='X1')
  tmp[is.na(tmp$X3.x),c('X3.x')]<-0
  tmp[is.na(tmp$X3.y),c('X3.y')]<-0
  sta <- cumsum(tmp$X3.x)-cumsum(tmp$X3.y)}
  
  return(sum(abs(sta)))
  
}

screen_npmle2 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))
  if (a$pf==1){
    
    data2 <-  transs(df2)
    tmp2 <-  merge(data2, basis,all.y = TRUE)
    sta <- cumsum(tmp2$X3)}
  else if  (b$pf==1){
    
    data1 <-  transs(df1)
    tmp1 <-  merge(data1, basis,all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  }
  else {
    
  data1 <-  transs(df1)
  data2 <-  transs(df2)
  tmp1 <- merge(data1, basis,all.y = TRUE)
  tmp2 <-  merge(data2, basis,all.y = TRUE)
  tmp <- merge(tmp1,tmp2,by='X1')
  tmp[is.na(tmp$X3.x),c('X3.x')]<-0
  tmp[is.na(tmp$X3.y),c('X3.y')]<-0
  sta <- cumsum(tmp$X3.x)-cumsum(tmp$X3.y)}
  
  return(sum(abs(sta)**2))
  
}
screen_npmle3 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))
  if (a$pf==1){
    
    data2 <-  transs(df2)
    tmp2 <-  merge(data2, basis,all.y = TRUE)
    sta <- cumsum(tmp2$X3)}
  else if  (b$pf==1){
    
    data1 <-  transs(df1)
    tmp1 <-  merge(data1, basis,all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  }
  else {data1 <-  transs(df1)
  data2 <-  transs(df2)
  tmp1 <- merge(data1, basis,all.y = TRUE)
  tmp2 <-  merge(data2, basis,all.y = TRUE)
  tmp <- merge(tmp1,tmp2,by='X1')
  tmp[is.na(tmp$X3.x),c('X3.x')]<-0
  tmp[is.na(tmp$X3.y),c('X3.y')]<-0
  sta <- cumsum(tmp$X3.x)-cumsum(tmp$X3.y)}
  
  return(max(abs(sta)))
}

int_split <- function(i,df){
  tmp1 <- as.numeric(df[i,][1])
  tmp2 <- as.numeric(df[i,][2])
  if (is.infinite(tmp1)) {start<- tmp2-1}
  else {start <- as.numeric(df[i,][1])}
  if (is.infinite(tmp2)) {end <- tmp1+1}
  else {end <- as.numeric(df[i,][2])}
  
  difs <- as.numeric(end-start-1)
  
  return (data.frame(cbind(start+c(0:difs),start+1+c(0:difs),rep(df[i,][3]/(difs+1),difs+1))))
}
transs <- function(df){
  dff <- df$X2-df$X1
  subs <- df[dff==1,c('X1','X2','X3')]  
  no <- which(dff!=1)
  if (length(no)==0) {res <- subs}
  else if (length(no)==1) {res <- data.frame(rbind(subs,int_split(no[1],df=df)))
  res$X1 <- as.numeric(res$X1)
  res$X2 <- as.numeric(res$X2)
  res$X3 <- as.numeric(res$X3)}
  else if (length(no)>=2) {ress <- data.frame(rbind(subs,int_split(no[1],df=df)))
  res <- data.frame(rbind(ress,int_split(no[2],df=df)))
  res$X1 <- as.numeric(res$X1)
  res$X2 <- as.numeric(res$X2)
  res$X3 <- as.numeric(res$X3)
  }
  return (res)
}

screen_ks<- function(i,y_mid,s)
{every = data.frame(cbind(y_mid,s[,i]))
part1 = every[every$V2==0,c('y_mid')]
part2 = every[every$V2==1,c('y_mid')]
return (as.numeric(ks.test(part1,part2)$statistic))
}




screen_dc <- function(i,y_mid,s)
{return (dcor(y_mid,s[,i]))}

screen_t <- function(i,y_mid,s)
{return (abs(cor(y_mid,s[,i])))}


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

screen_mv <- function(i,s,y_mid)
{return (abs(MV(Y=s[,i],Xk=y_mid)))}








Simu_Multi_Norm<-function(x_len, sd = sd, pho =pho){
  #初始化协方差矩阵
  V <- matrix(data = NA, nrow = x_len, ncol = x_len)
  
  #mean及sd分别为随机向量x的均值和方差
  
  #对协方差矩阵进行赋值pho(i,j) = pho^|i-j|
  for(i in 1:x_len){ ##遍历每一行
    for(j in 1:x_len){ ##遍历每一列
      V[i,j] <- pho^abs(i-j)
    }
  }
  
  V<-(sd^2) * V
  return(V)
}
sigma1 <- Simu_Multi_Norm(x_len = 1000,sd  = 1, pho = 0.2)
sigma2 <- Simu_Multi_Norm(x_len = 1000,sd  = 1, pho = 0.5)
sigma3 <- Simu_Multi_Norm(x_len = 1000,sd  = 1, pho = 0.8)
simu1 <- function(kk,pp,nn,rat){
  print (kk)
  set.seed(kk)
  x <-mvrnorm(nn, mu = rep(0,1000),pp)
  s <- x>0
  prob <- 0.2
  #beta =  c(c(0.5,0.5,0.5,0.5,0.5),rep(0,995))
  beta = c(c(1,0.8,0.6),rep(0,15),c(-0.6,-0.8),rep(0,980))
  lambda <- (s %*% beta)
  y <- rpois(n=nn,exp(lambda))
  l = y-rpois(n=nn,2)-1
  r = y+2*floor(rchisq(df=3,n=nn))+1
  tot <- prob*nn
  nos <- sample(c(1:nn),size=tot)
  tt <- rat*(prob*nn)
  r[nos[1:tt]] <- Inf
  l[nos[tt+1:tot]] <- -Inf
  transform <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=0.5*(l[i]+r[i])}
    return (y)
  }
  transform2 <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=2/3*l[i]+1/3*r[i]}
    return (y)
  }
  transform3 <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=y=1/3*l[i]+2/3*r[i]}
    return (y)
  }
  y_mid <- sapply(c(1:nn),transform,l=l,r=r)
  y_3 <- sapply(c(1:nn),transform2,l=l,r=r)
  y_6 <- sapply(c(1:nn),transform3,l=l,r=r)
  intervals <- initcomputeMLE(l,r)
  
  max_tmp <- max(intervals$intmap[2,])
  min_tmp <- min(intervals$intmap[1,])
  if (is.finite(max_tmp)) {end <- max_tmp} 
  if (is.infinite(max_tmp)) {end <-max(intervals$intmap[1,])+1 }
  if (is.finite(min_tmp)) {start <- min_tmp} 
  if (is.infinite(min_tmp)) {start <-min(intervals$intmap[2,])-1 }
  
  X1<-c(start:as.numeric(end-1))
  X2<-c(as.numeric(start+1):end)
  basis<- cbind(X1,X2)
  a1 = Sys.time()
  s_npmle1 <- sapply(c(1:1000),screen_npmle1,l=l,r=r,s=s,basis=basis)
  b1 = Sys.time()
  t1 = b1-a1
  
  a2 = Sys.time()
  s_npmle2 <- sapply(c(1:1000),screen_npmle2,l=l,r=r,s=s,basis=basis)
  b2 = Sys.time()
  t2 = b2-a2
  
  a3 = Sys.time()
  s_npmle3 <- sapply(c(1:1000),screen_npmle3,l=l,r=r,s=s,basis=basis)
  b3 = Sys.time()
  t3 = b3-a3
  
  a4 = Sys.time()
  s_ks <- sapply(c(1:1000),screen_ks,y_mid=y_mid,s=s)
  b4 = Sys.time()
  t4 = b4-a4
  
  a5 = Sys.time()
  s_dc <- sapply(c(1:1000),screen_dc,y_mid=y_mid,s=s)
  b5 = Sys.time()
  t5 = b5-a5
  
  a6 = Sys.time()
  s_mv <- sapply(c(1:1000),screen_mv,y_mid=y_mid,s=s)
  b6 = Sys.time()
  t6 = b6-a6
  
  a7 = Sys.time()
  s_ks3 <- sapply(c(1:1000),screen_ks,y_mid=y_3,s=s)
  b7 = Sys.time()
  t7 = b7-a7
  
  a8 = Sys.time()
  s_dc3 <- sapply(c(1:1000),screen_dc,y_mid=y_3,s=s)
  b8 = Sys.time()
  t8 = b8-a8
  
  a9 = Sys.time()
  s_mv3 <- sapply(c(1:1000),screen_mv,y_mid=y_3,s=s)
  b9 = Sys.time()
  t9 = b9-a9
  
  a10 = Sys.time()
  s_ks6 <- sapply(c(1:1000),screen_ks,y_mid=y_6,s=s)
  b10 = Sys.time()
  t10 = b10-a10
  
  a11 = Sys.time()
  s_dc6 <- sapply(c(1:1000),screen_dc,y_mid=y_6,s=s)
  b11 = Sys.time()
  t11 = b11-a11
  
  a12 = Sys.time()
  s_mv6 <- sapply(c(1:1000),screen_mv,y_mid=y_6,s=s)
  b12 = Sys.time()
  t12 = b12-a12
  
  a <- (1001-rank(s_npmle1))[c(1,2,3,19,20)]
  b <- (1001-rank(s_npmle2))[c(1,2,3,19,20)]
  c <- (1001-rank(s_npmle3))[c(1,2,3,19,20)]
  d<- (1001-rank(s_ks))[c(1,2,3,19,20)]
  e <- (1001-rank(s_dc))[c(1,2,3,19,20)]
  f <- (1001-rank(s_mv))[c(1,2,3,19,20)]
  d3<- (1001-rank(s_ks3))[c(1,2,3,19,20)]
  e3 <- (1001-rank(s_dc3))[c(1,2,3,19,20)]
  f3 <- (1001-rank(s_mv3))[c(1,2,3,19,20)]
  d6<- (1001-rank(s_ks6))[c(1,2,3,19,20)]
  e6 <- (1001-rank(s_dc6))[c(1,2,3,19,20)]
  f6 <- (1001-rank(s_mv6))[c(1,2,3,19,20)]
  times = (c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12))
  return (c(a,b,c,d,e,f,d3,e3,f3,d6,e6,f6,times))
  
}
simu2 <- function(kk,pp,nn,rat){
  print (kk)
  set.seed(kk)
  x <-mvrnorm(nn, mu = rep(0,1000),pp)
  s <- x>as.numeric(sapply(c(1:1000),function(i) quantile(x[,i],0.3)))
  prob <- 0.2
  #beta =  c(c(0.5,0.5,0.5,0.5,0.5),rep(0,995))
  beta = c(c(1,0.8,0.6),rep(0,15),c(-0.6,-0.8),rep(0,980))
  lambda <- (s %*% beta)
  y <- rpois(n=nn,exp(lambda))
  l = y-rpois(n=nn,2)
  r = y+2*floor(rchisq(df=3,n=nn))+1
  tot <- prob*nn
  nos <- sample(c(1:nn),size=tot)
  tt <- rat*(prob*nn)
  r[nos[1:tt]] <- Inf
  l[nos[tt+1:tot]] <- -Inf
  transform <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=0.5*(l[i]+r[i])}
    return (y)
  }
  transform2 <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=2/3*l[i]+1/3*r[i]}
    return (y)
  }
  transform3 <- function(i,l=l,r=r){
    if (is.infinite(l[i])) {y = r[i]}
    else if (is.infinite(r[i])) {y=l[i]}
    else {y=y=1/3*l[i]+2/3*r[i]}
    return (y)
  }
  y_mid <- sapply(c(1:nn),transform,l=l,r=r)
  y_3 <- sapply(c(1:nn),transform2,l=l,r=r)
  y_6 <- sapply(c(1:nn),transform3,l=l,r=r)
  intervals <- initcomputeMLE(l,r)
  
  max_tmp <- max(intervals$intmap[2,])
  min_tmp <- min(intervals$intmap[1,])
  if (is.finite(max_tmp)) {end <- max_tmp} 
  if (is.infinite(max_tmp)) {end <-max(intervals$intmap[1,])+1 }
  if (is.finite(min_tmp)) {start <- min_tmp} 
  if (is.infinite(min_tmp)) {start <-min(intervals$intmap[2,])-1 }
  
  X1<-c(start:as.numeric(end-1))
  X2<-c(as.numeric(start+1):end)
  basis<- cbind(X1,X2)
  a1 = Sys.time()
  s_npmle1 <- sapply(c(1:1000),screen_npmle1,l=l,r=r,s=s,basis=basis)
  b1 = Sys.time()
  t1 = b1-a1
  
  a2 = Sys.time()
  s_npmle2 <- sapply(c(1:1000),screen_npmle2,l=l,r=r,s=s,basis=basis)
  b2 = Sys.time()
  t2 = b2-a2
  
  a3 = Sys.time()
  s_npmle3 <- sapply(c(1:1000),screen_npmle3,l=l,r=r,s=s,basis=basis)
  b3 = Sys.time()
  t3 = b3-a3
  
  a4 = Sys.time()
  s_ks <- sapply(c(1:1000),screen_ks,y_mid=y_mid,s=s)
  b4 = Sys.time()
  t4 = b4-a4
  
  a5 = Sys.time()
  s_dc <- sapply(c(1:1000),screen_dc,y_mid=y_mid,s=s)
  b5 = Sys.time()
  t5 = b5-a5
  
  a6 = Sys.time()
  s_mv <- sapply(c(1:1000),screen_mv,y_mid=y_mid,s=s)
  b6 = Sys.time()
  t6 = b6-a6
  
  a7 = Sys.time()
  s_ks3 <- sapply(c(1:1000),screen_ks,y_mid=y_3,s=s)
  b7 = Sys.time()
  t7 = b7-a7
  
  a8 = Sys.time()
  s_dc3 <- sapply(c(1:1000),screen_dc,y_mid=y_3,s=s)
  b8 = Sys.time()
  t8 = b8-a8
  
  a9 = Sys.time()
  s_mv3 <- sapply(c(1:1000),screen_mv,y_mid=y_3,s=s)
  b9 = Sys.time()
  t9 = b9-a9
  
  a10 = Sys.time()
  s_ks6 <- sapply(c(1:1000),screen_ks,y_mid=y_6,s=s)
  b10 = Sys.time()
  t10 = b10-a10
  
  a11 = Sys.time()
  s_dc6 <- sapply(c(1:1000),screen_dc,y_mid=y_6,s=s)
  b11 = Sys.time()
  t11 = b11-a11
  
  a12 = Sys.time()
  s_mv6 <- sapply(c(1:1000),screen_mv,y_mid=y_6,s=s)
  b12 = Sys.time()
  t12 = b12-a12
  
  a <- (1001-rank(s_npmle1))[c(1,2,3,19,20)]
  b <- (1001-rank(s_npmle2))[c(1,2,3,19,20)]
  c <- (1001-rank(s_npmle3))[c(1,2,3,19,20)]
  d<- (1001-rank(s_ks))[c(1,2,3,19,20)]
  e <- (1001-rank(s_dc))[c(1,2,3,19,20)]
  f <- (1001-rank(s_mv))[c(1,2,3,19,20)]
  d3<- (1001-rank(s_ks3))[c(1,2,3,19,20)]
  e3 <- (1001-rank(s_dc3))[c(1,2,3,19,20)]
  f3 <- (1001-rank(s_mv3))[c(1,2,3,19,20)]
  d6<- (1001-rank(s_ks6))[c(1,2,3,19,20)]
  e6 <- (1001-rank(s_dc6))[c(1,2,3,19,20)]
  f6 <- (1001-rank(s_mv6))[c(1,2,3,19,20)]
  times = (c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12))
  return (c(a,b,c,d,e,f,d3,e3,f3,d6,e6,f6,times))
}


#res1_ <- sapply(c(1:2),simu1,pp=sigma1,n=200,rat=0.3)

res1 <- sapply(c(1:200),simu1,pp=sigma1,n=200,rat=0.3)
res2 <- sapply(c(1:200),simu1,pp=sigma2,n=200,rat=0.3)
#res3 <- sapply(c(1:200),simu1,pp=sigma3,n=200,rat=0.3)



res4 <- sapply(c(1:200),simu1,pp=sigma1,n=200,rat=0.7)
res5 <- sapply(c(1:200),simu1,pp=sigma2,n=200,rat=0.7)
#res6 <- sapply(c(1:200),simu1,pp=sigma3,n=200,rat=0.7)


res7 <- sapply(c(1:200),simu2,pp=sigma1,n=200,rat=0.3)
res8 <- sapply(c(1:200),simu2,pp=sigma2,n=200,rat=0.3)
#res9 <- sapply(c(1:200),simu2,pp=sigma3,n=200,rat=0.3)
res10 <- sapply(c(1:200),simu2,pp=sigma1,n=200,rat=0.7)
res11 <- sapply(c(1:200),simu2,pp=sigma2,n=200,rat=0.7)
#res12 <- sapply(c(1:200),simu2,pp=sigma3,n=200,rat=0.7)

write.csv(res1,"Full_res1_Time.csv")
write.csv(res2,"Full_res2_Time.csv")
#write.csv(res3,"Full_res3.csv")
write.csv(res4,"Full_res4_Time.csv")
write.csv(res5,"Full_res5_Time.csv")
#write.csv(res6,"Full_res6.csv")
write.csv(res7,"Full_res7_Time.csv")
write.csv(res8,"Full_res8_Time.csv")
#write.csv(res9,"Full_res9.csv")
write.csv(res10,"Full_res10_Time.csv")
write.csv(res11,"Full_res11_Time.csv")
#write.csv(res12,"Full_res12.csv")

cbind(median(sapply(c(1:200),function(i) max(res1[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res1[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res1[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res1[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res1[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res1[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res1[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res1[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res1[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res1[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res1[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res1[56:60,i]))))

cbind(median(sapply(c(1:200),function(i) max(res2[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res2[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res2[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res2[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res2[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res2[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res2[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res2[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res2[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res2[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res2[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res2[56:60,i]))))


cbind(median(sapply(c(1:200),function(i) max(res4[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res4[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res4[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res4[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res4[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res4[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res4[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res4[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res4[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res4[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res4[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res4[56:60,i]))))

cbind(median(sapply(c(1:200),function(i) max(res5[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res5[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res5[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res5[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res5[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res5[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res5[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res5[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res5[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res5[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res5[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res5[56:60,i]))))



cbind(median(sapply(c(1:200),function(i) max(res7[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res7[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res7[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res7[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res7[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res7[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res7[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res7[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res7[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res7[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res7[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res7[56:60,i]))))

cbind(median(sapply(c(1:200),function(i) max(res8[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res8[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res8[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res8[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res8[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res8[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res8[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res8[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res8[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res8[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res8[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res8[56:60,i]))))



cbind(median(sapply(c(1:200),function(i) max(res10[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res10[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res10[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res10[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res10[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res10[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res10[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res10[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res10[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res10[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res10[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res10[56:60,i]))))
  
cbind(median(sapply(c(1:200),function(i) max(res11[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res11[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res11[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res11[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res11[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res11[26:30,i]))),
      median(sapply(c(1:200),function(i) max(res11[31:35,i]))),
      median(sapply(c(1:200),function(i) max(res11[36:40,i]))),
      median(sapply(c(1:200),function(i) max(res11[41:45,i]))),
      median(sapply(c(1:200),function(i) max(res11[46:50,i]))),
      median(sapply(c(1:200),function(i) max(res11[51:55,i]))),
      median(sapply(c(1:200),function(i) max(res11[56:60,i]))))


cbind(median(sapply(c(1:200),function(i) max(res2[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res2[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res2[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res2[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res2[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res2[26:30,i])))      )

cbind(median(sapply(c(1:200),function(i) max(res3[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res3[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res3[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res3[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res3[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res3[26:30,i])))      )

cbind(median(sapply(c(1:200),function(i) max(res4[1:5,i]))),
      median(sapply(c(1:200),function(i) max(res4[6:10,i]))),
      median(sapply(c(1:200),function(i) max(res4[11:15,i]))),
      median(sapply(c(1:200),function(i) max(res4[16:20,i]))),
      median(sapply(c(1:200),function(i) max(res4[21:25,i]))),
      median(sapply(c(1:200),function(i) max(res4[26:30,i])))      )


c((as.numeric(quantile(sapply(c(1:200),function(i) max(res1[1:5,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[1:5,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res1[6:10,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[6:10,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res1[11:15,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[11:15,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res1[16:20,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[16:20,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res1[21:25,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[21:25,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res1[26:30,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res1[26:30,i])),0.25)))/1.34)                                                                                                    

c((as.numeric(quantile(sapply(c(1:200),function(i) max(res2[1:5,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[1:5,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res2[6:10,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[6:10,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res2[11:15,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[11:15,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res2[16:20,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[16:20,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res2[21:25,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[21:25,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res2[26:30,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res2[26:30,i])),0.25)))/1.34)                                                                                                    

c((as.numeric(quantile(sapply(c(1:200),function(i) max(res3[1:5,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[1:5,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res3[6:10,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[6:10,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res3[11:15,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[11:15,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res3[16:20,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[16:20,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res3[21:25,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[21:25,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res3[26:30,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res3[26:30,i])),0.25)))/1.34)                                                                                                    

c((as.numeric(quantile(sapply(c(1:200),function(i) max(res4[1:5,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[1:5,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res4[6:10,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[6:10,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res4[11:15,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[11:15,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res4[16:20,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[16:20,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res4[21:25,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[21:25,i])),0.25)))/1.34,
  (as.numeric(quantile(sapply(c(1:200),function(i) max(res4[26:30,i])),0.75))-as.numeric(quantile(sapply(c(1:200),function(i) max(res4[26:30,i])),0.25)))/1.34)                                                                                                    

cbind(c(sum(sapply(c(1:200),function(i) max(res1[1:5,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res1[6:10,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res1[11:15,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res1[16:20,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res1[21:25,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res1[26:30,i])<37))/200),
      c(sum(sapply(c(1:200),function(i) max(res1[1:5,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res1[6:10,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res1[11:15,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res1[16:20,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res1[21:25,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res1[26:30,i])<75))/200),
      c(sum(sapply(c(1:200),function(i) max(res1[1:5,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res1[6:10,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res1[11:15,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res1[16:20,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res1[21:25,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res1[26:30,i])<113))/200))



cbind(c(sum(sapply(c(1:200),function(i) max(res2[1:5,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res2[6:10,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res2[11:15,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res2[16:20,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res2[21:25,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res2[26:30,i])<37))/200),
      c(sum(sapply(c(1:200),function(i) max(res2[1:5,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res2[6:10,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res2[11:15,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res2[16:20,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res2[21:25,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res2[26:30,i])<75))/200),
      c(sum(sapply(c(1:200),function(i) max(res2[1:5,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res2[6:10,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res2[11:15,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res2[16:20,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res2[21:25,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res2[26:30,i])<113))/200))


cbind(c(sum(sapply(c(1:200),function(i) max(res3[1:5,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res3[6:10,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res3[11:15,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res3[16:20,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res3[21:25,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res3[26:30,i])<37))/200),
      c(sum(sapply(c(1:200),function(i) max(res3[1:5,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res3[6:10,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res3[11:15,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res3[16:20,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res3[21:25,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res3[26:30,i])<75))/200),
      c(sum(sapply(c(1:200),function(i) max(res3[1:5,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res3[6:10,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res3[11:15,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res3[16:20,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res3[21:25,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res3[26:30,i])<113))/200))



cbind(c(sum(sapply(c(1:200),function(i) max(res4[1:5,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res4[6:10,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res4[11:15,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res4[16:20,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res4[21:25,i])<37))/200,
        sum(sapply(c(1:200),function(i) max(res4[26:30,i])<37))/200),
      c(sum(sapply(c(1:200),function(i) max(res4[1:5,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res4[6:10,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res4[11:15,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res4[16:20,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res4[21:25,i])<75))/200,
        sum(sapply(c(1:200),function(i) max(res4[26:30,i])<75))/200),
      c(sum(sapply(c(1:200),function(i) max(res4[1:5,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res4[6:10,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res4[11:15,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res4[16:20,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res4[21:25,i])<113))/200,
        sum(sapply(c(1:200),function(i) max(res4[26:30,i])<113))/200))

rbind(sapply(c(1:5),function(i) sum(res1[i,]<37)/200),
      sapply(c(6:10),function(i) sum(res1[i,]<37)/200),
      sapply(c(11:15),function(i) sum(res1[i,]<37)/200),
      sapply(c(16:20),function(i) sum(res1[i,]<37)/200),
      sapply(c(21:25),function(i) sum(res1[i,]<37)/200),
      sapply(c(26:30),function(i) sum(res1[i,]<37)/200))


rbind(sapply(c(1:5),function(i) sum(res2[i,]<37)/200),
      sapply(c(6:10),function(i) sum(res2[i,]<37)/200),
      sapply(c(11:15),function(i) sum(res2[i,]<37)/200),
      sapply(c(16:20),function(i) sum(res2[i,]<37)/200),
      sapply(c(21:25),function(i) sum(res2[i,]<37)/200),
      sapply(c(26:30),function(i) sum(res2[i,]<37)/200))


rbind(sapply(c(1:5),function(i) sum(res3[i,]<37)/200),
      sapply(c(6:10),function(i) sum(res3[i,]<37)/200),
      sapply(c(11:15),function(i) sum(res3[i,]<37)/200),
      sapply(c(16:20),function(i) sum(res3[i,]<37)/200),
      sapply(c(21:25),function(i) sum(res3[i,]<37)/200),
      sapply(c(26:30),function(i) sum(res3[i,]<37)/200))

rbind(sapply(c(1:5),function(i) sum(res4[i,]<37)/200),
      sapply(c(6:10),function(i) sum(res4[i,]<37)/200),
      sapply(c(11:15),function(i) sum(res4[i,]<37)/200),
      sapply(c(16:20),function(i) sum(res4[i,]<37)/200),
      sapply(c(21:25),function(i) sum(res4[i,]<37)/200),
      sapply(c(26:30),function(i) sum(res4[i,]<37)/200))
################################################################
################################################################
################################################################


