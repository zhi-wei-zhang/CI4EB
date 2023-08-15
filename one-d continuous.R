
rm(list=ls())

########### one dimensional continuous X ###########

f.pall = function(seed, nsample, parm){
  
  library(SuperLearner)
  library(np)

  f.sim.contin = function(n, parm){
    #---- generate data -----
    distr.y = parm$distr.y 
    qx = parm$qx
    
    sig.w1 = sig.w2 = parm$sig.w
    mu.w1 = mu.w2 = parm$mu.w
    sig.y = parm$sig.y
    x.mean = parm$x.mean
    x.sd = parm$x.sd
    bet = parm$bet 
    B = parm$B
    DELTA = parm$DELTA
    
    SL_library = parm$SL_library
    
    x = rnorm(n, x.mean, x.sd)
    w1 = rnorm(n, mu.w1, sig.w1)
    w2 = rnorm(n, mu.w2, sig.w2)
    
    prob.a = 0.5
    a = rbinom(n, 1, prob.a)
    nqx = length(qx)
    
    if (distr.y=="gaussian"){
      y = rnorm(n, cbind(1, x, w1, w2, x*w1, x*a) %*%bet, sig.y)  
      FAMILY = gaussian() 
    }
    if (distr.y=="binomial"){
      y = rbinom(n, 1, plogis(cbind(1, x, w1, w2, x*w1, x*a) %*%bet))  
      FAMILY = binomial() 
    }
    
    dat = data.frame(y=y, x=x, w1 = w1, w2 = w2, a = a)
    
    alp = parm$alp
    #------- start -----
    
    f_h.hat = function(dat, emeth="nonpar", family=gaussian(), SL_library, distr.y){
      #----h.hat (w)
      #family = gaussian() or binomial()
      dat0 = dat1 = dat
      dat0$a = 0
      dat1$a = 1
      p = mean(dat$a)   
      
      if(emeth=="true"){#--true h(w): 0. h(w)=0
        pred  =cbind(model.matrix(y~x+w1+w2+x:w1+x:a, data = dat0)%*%bet, 
                     model.matrix(y~x+w1+w2+x:w1+x:a, data = dat1)%*%bet)
        if (distr.y=="binomial") pred = plogis(pred)
      }
      if(emeth=="non"){ #--estimate h(w): 1. h(w)=0
        pred = cbind(rep(0,n),rep(0,n))
      }
      if(emeth=="par"){#--estimate h(w): 2. parametric 
        par_fit = glm(y~x+w1+w2+x:w1+x:a, data = dat, family = family)
        pred  = cbind(predict(par_fit, newdata=dat0),predict(par_fit, newdata=dat1))
        if (distr.y=="binomial") pred = plogis(pred)
      }
      if(emeth== "nonpar"){#--estimate h(w): 3. nonparametric
        
        
        sl_fit = SuperLearner(Y = dat$y, X = subset(dat,select=c(-y)), family = family,
                              SL.library = SL_library) 
        
        pred  =cbind(predict(sl_fit, newdata=dat0[,-1])$pred, predict(sl_fit, newdata=dat1[,-1])$pred)
      }
      h.hat = as.vector(pred%*%c(p,1-p))
      return(h.hat)
    }
    
    f_B = function(dat, emeth= "nonpar", family=gaussian(), SL_library, distr.y){
      
      h.hat = f_h.hat(dat, emeth, family, SL_library, distr.y)
      p = mean(dat$a)
      D.hat= (dat$a/p-(1-dat$a)/(1-p))*(dat$y-h.hat)
      
      bws = 1
      bw.fit = npregbw(D.hat~x, data = data.frame(h.hat = h.hat,x=dat$x), regtype ="ll",
                       ckertype = "gaussian",  bws=bws, bandwidth.compute=F )
      delt.fit = npreg(bw.fit)
      # summary(bw.fit)
      
      bn = bw.fit$bw 
      
      
      delt = predict(delt.fit, newdata = data.frame(x=qx),se.fit=T)
      delta.hat = delt$fit
      sigma.hat = delt$se.fit
      
      xu = x.mean+3*x.sd 
      xl = x.mean-3*x.sd 
      qn = 2*log((xu-xl)/bn)-log(2)-2*log(2*pi)
      
      Ca = sqrt(qn-2*log(-1/2*log(1-alp)))
      
      DL = (delta.hat - Ca*sigma.hat)
      DU = (delta.hat + Ca*sigma.hat)
      
      BL = (DL>0)
      BU = (DU>0)
      
      return(data.frame(BL, BU, DL, DU))
    }
    
    #-- run
    B0 = f_B(dat, emeth= "true", family=FAMILY, SL_library, distr.y) 
    B1 = f_B(dat, emeth= "non", family=FAMILY, SL_library, distr.y) 
    B2 = f_B(dat, emeth= "par", family=FAMILY, SL_library, distr.y)
    B3 = f_B(dat, emeth= "nonpar", family=FAMILY, SL_library, distr.y)
    
    cb = list(true = B0, zero = B1, par = B2,  nonpar = B3)
    
    f.cp.continuous = function(x, B, DELTA){
      BL = (x[,1]==1); BU=(x[,2]==1)
      cp = ifelse(all(B[BL],BU[B]),1,0)
      scp.L = ifelse(all(B[BL]),1,0)
      scp.U = ifelse(all(BU[B]),1,0)
      size = mean(BU&(!BL))
      return(c(cp, scp.L, scp.U, size))
    }
    
    res = sapply(cb, f.cp.continuous, B=B, DELTA=DELTA)
    return(c(res))
  }
  
  #------- setting  ----
  set.seed(seed)
  
  n=nsample
  
  res = f.sim.contin(n, parm)
  
  return(res)
}


#-------seting------
nsample= 200

{
  parm = list(  
    distr.y = "gaussian", # "binomial",
    x.mean = 1, x.sd = 2, 
    mu.w=0, sig.w = 0.5, sig.y = 0.5,
    bet = c(0, 0.2, -0.2, 0.2, 0.1, 0.2)*3, # cbind(1, x, w1, w2, x*w1, x*a) %*%bet
    alp = 0.05,
    SL_library = c("SL.lm","SL.glm", "SL.glm.interaction", 
                   "SL.randomForest")
  )
  distr.y = parm$distr.y #"gaussian"
  
  qx = qnorm((1:99)*0.01, parm$x.mean, parm$x.sd)
  parm$qx = qx
  
  sig.w1 = sig.w2 = parm$sig.w
  mu.w1 = mu.w2 = parm$mu.w
  bet = parm$bet 
  
  if (distr.y=="gaussian"){
    f_delta.true = function(x){cbind(0, 0, 0, 0, 0, x) %*% bet}
  }
  if (distr.y=="binomial"){
    f_delta.true = function(x){
      
      f2 = function(w1){
        f1 = function(w2,w1){(plogis(cbind(1, x, w1, w2, x*w1, x) %*%bet)
                              -plogis(cbind(1, x, w1, w2, x*w1, 0) %*%bet))*dnorm(w1, mu.w1, sig.w1)*dnorm(w2, mu.w2, sig.w2)}
        t1 = sapply(w1, function(w1){integrate(f1,-Inf,Inf,w1=w1)$value})
      }
      delta.true = integrate(f2,-Inf,Inf)$value
    }
  }
  DELTA = sapply(qx, f_delta.true)
  B = (DELTA>0)
  
  parm$DELTA = DELTA
  parm$B=B
  
}

#-------run------
print( f.pall(123, nsample= nsample, parm = parm) )

