
rm(list=ls())

########### two continuous X ###########

f.pall = function(seed, nsample, parm){
 
  library(SuperLearner)
  library(np)
  library(rgl)
  
  f.sim.gener = function(n, parm){
    #---- generate data -----
    
    {
      distr.y = parm$distr.y 
      qx1 = parm$qx1
      qx2 = parm$qx2
      
      x1.mean=parm$x1.mean
      x1.sd=parm$x1.sd
      x2.mean=parm$x2.mean
      x2.sd=parm$x2.sd
      sig.w1 = sig.w2 = parm$sig.w
      mu.w1 = mu.w2 = parm$mu.w
      bet = parm$bet  
      B = parm$B
      DELTA = parm$DELTA
      sig.y = parm$sig.y
      SL_library = parm$SL_library
      SL_library_wr = parm$SL_library_wr
      alp = parm$alp
      #----- some functions ----
      # linear model for delta
      lin.delta = function(X,a) as.vector(X%*%a) 
      lin.der = function(X,a) X
      
      # model based on a logit link for delta
      lgt.delta = function(X,a) {
        u = as.vector(exp(X%*%a))
        (u-1)/(u+1)
      }
      lgt.der = function(X,a) {
        p = as.vector(plogis(X%*%a))
        2*p*(1-p)*X	
      }
      
      #---- generate data -----
      x1 = rnorm(n, x1.mean, x1.sd)
      x2 = rnorm(n, x2.mean, x2.sd)
      w1 = rnorm(n, mu.w1, sig.w1)
      w2 = rnorm(n, mu.w2, sig.w2)
    }
    
    prob.a = 0.5
    a = rbinom(n, 1, prob.a)
    
    if (distr.y=="gaussian"){
      y = rnorm(n, cbind(1, x1, x2, w1, w2, x1*w1, x1*a, x2*a)%*%bet,sig.y)  
      FAMILY = gaussian() 
      modl = list(delta = lin.delta, der = lin.der)

    }
    if (distr.y=="binomial"){
      y = rbinom(n, 1, plogis(cbind(1, x1, x2, w1, w2, x1*w1, x1*a, x2*a)%*%bet))
      FAMILY = binomial() 
      modl = list(delta = lgt.delta, der = lgt.der)
    }
    dat = data.frame(y=y, x1=x1, x2=x2, w1 = w1, w2 = w2, a = a)
    
    #------- start -----
    
    f_B.genl = function(dat, emeth = "nonpar", family=gaussian(), modl, Vh=NULL, bet0=NULL, SL_library){
      #----h.hat (w)
      #family = gaussian() or binomial()
      dat0 = dat1 = dat
      dat0$a = 0
      dat1$a = 1
      p = mean(dat$a)
      if(emeth=="true"){#--true h(w): 0. h(w)=0
        pred  =cbind(model.matrix(y~x1+x2+w1+w2+x1:w1+x1:a+x2:a, data = dat0)%*%bet, 
                     model.matrix(y~x1+x2+w1+w2+x1:w1+x1:a+x2:a, data = dat1)%*%bet)
        if (distr.y=="binomial") pred = plogis(pred)
      }
      if(emeth=="non"){ #--estimate h(w): 1. h(w)=0
        pred = cbind(rep(0,n),rep(0,n))
      }
      if(emeth=="par"){#--estimate h(w): 2. parametric
        par_fit = glm(y~x1+x2+w1+w2+x1:w1+x1:a+x2:a, data = dat, family = family)
        pred  =cbind(predict(par_fit, newdata=dat0),predict(par_fit, newdata=dat1))
        if (distr.y=="binomial") pred = plogis(pred)
      }
      if(emeth== "nonpar"){#--estimate h(w): 3. nonparametric
        
        sl_fit = SuperLearner(Y = dat$y, X = subset(dat,select=c(-y)), family = family,
                              SL.library = SL_library)
        pred  =cbind(predict(sl_fit, newdata=dat0[,-1])$pred, predict(sl_fit, newdata=dat1[,-1])$pred)
        
      }
     
      h.hat = as.vector(pred%*%c(p,1-p))

      
      #---- Delta
      p = mean(dat$a)
      D.hat= (dat$a/p-(1-dat$a)/(1-p))*(dat$y-h.hat)
      
      X = cbind(1,dat$x1,dat$x2) 
      np = ncol(X)
      n = nrow(X)
      
      if (is.null(bet0)) bet0 = rep(0,np)
      if (is.null(Vh)) Vh = X
      
      g = modl$delta
      f.EE = function(bet) {
        delt = as.vector(g(X,bet))
        U = Vh*(D.hat-delt)
        sum(colSums(U)^2)
      }
      fit = optim(bet0, f.EE)
      bet = matrix(fit$par,ncol=1)
      
      der.delta = modl$der(X,bet)
      Der = -t(Vh)%*%der.delta/n
      
      EE = Vh*(D.hat-as.vector(g(X,bet)))
      Ver = t(EE)%*%EE/n
      
      inv.Der = solve(Der)
      Sig = inv.Der%*%Ver%*%t(inv.Der)/n
      
      Ca = sqrt(qchisq(1-alp, np))
      
      nqx = length(qx1)
      qx = cbind(1,rep(qx1, nqx),rep(qx2, each=nqx))
      ss = matrix(sqrt(diag(qx%*%Sig%*%t(qx))),ncol=1)

      DL = g(qx%*%bet-ss*Ca, 1)
      DU = g(qx%*%bet+ss*Ca, 1)
      
      BL = (DL>0)      
      BU = (DU>0)
      
      return(data.frame(BL, BU, DL, DU))
    }
    
    #-- run
    B0 = f_B.genl(dat, emeth= "true", family=FAMILY, modl, Vh=NULL, bet0=NULL, SL_library) 
    B1 = f_B.genl(dat, emeth= "non", family=FAMILY, modl, Vh=NULL, bet0=NULL, SL_library)
    B2 = f_B.genl(dat, emeth= "par", family=FAMILY, modl, Vh=NULL, bet0=NULL, SL_library)
    B3 = f_B.genl(dat, emeth= "nonpar", family=FAMILY, modl, Vh=NULL, bet0=NULL, SL_library) 

   
    cb = list(true = B0, zero = B1, par = B2, nonpar = B3)

    f.cp.continuous = function(x, B){
      BL = x[,1]; BU=x[,2]
      cp = ifelse(all(B[BL],BU[B]),1,0)
      scp.L = ifelse(all(B[BL]),1,0)
      scp.U = ifelse(all(BU[B]),1,0)
      size = mean(BU&(!BL))
      return(c(cp, scp.L, scp.U, size))
    }
    res = sapply(cb, f.cp.continuous, B=B)
    return(c(res))
  }
  
  #------- setting -------
  set.seed(seed)

  n=nsample
  
  res = f.sim.gener(n, parm)
  
  return(res)
}



#-------seting------
nsample= 200
{
  parm = list(  
                distr.y = "gaussian", #"binomial",
                x1.mean = 1, x1.sd = 2, x2.mean = -1, x2.sd = 1,
                sig.w = 0.5, mu.w=0, sig.y=0.5,
                bet = c(0, 0.2, 0.2, 0.2, -0.2, 0.2, 0.2, 0.4)*3,
                alp = 0.05,
                SL_library = c("SL.lm","SL.glm", "SL.glm.interaction",
                               "SL.randomForest")
  )
  distr.y = parm$distr.y 
  
  qx1 = qnorm((1:99)*0.01, parm$x1.mean, parm$x1.sd)
  parm$qx1 = qx1
  
  qx2 = qnorm((1:99)*0.01, parm$x2.mean, parm$x2.sd)
  parm$qx2 = qx2
  
  sig.w1 = sig.w2 = parm$sig.w
  mu.w1 = mu.w2 = parm$mu.w
  bet = parm$bet 
  
  if (distr.y=="gaussian"){
    f_delta.true = function(x1,x2){cbind(0, 0, 0, 0, 0, 0, x1, x2) %*% bet}
  }
  if (distr.y=="binomial"){
    f_delta.true = function(x1,x2){
      
      f2 = function(w1){
        f1 = function(w2,w1){(plogis(cbind(1, x1, x2, w1, w2, x1*w1, x1, x2) %*%bet)
                              -plogis(cbind(1, x1, x2, w1, w2, x1*w1, 0, 0) %*%bet))*dnorm(w1, mu.w1, sig.w1)*dnorm(w2, mu.w2, sig.w2)}
        t1 = sapply(w1, function(w1){integrate(f1,-4*parm$sig.w,4*parm$sig.w,w1=w1)$value})
      }
      delta.true = integrate(f2,-4*parm$sig.w,4*parm$sig.w)$value
    }
  }
  
  DELTA = matrix(NA, length(qx1),length(qx2))
  for(i in 1:length(qx1)){
    for(j in 1:length(qx2)){
      DELTA[i,j] = f_delta.true(qx1[i],qx2[j]) 
    }
  }
  B = (DELTA>0)
  
  parm$DELTA = DELTA
  parm$B=B
}

#---run

print(f.pall(3, nsample= nsample, parm = parm))

