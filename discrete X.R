

########### discrete X ###########

rm(list=ls())

f.pall = function(seed, nsample, parm){
 
  library(SuperLearner)

  
  f.sim.discr = function(n, parm){
    #---- generate data -----
    
    distr.y = parm$distr.y 
    K = parm$K
    prob.x = parm$prob.x
    field.x = parm$field.x
    
    sig.w1 = sig.w2 = parm$sig.w
    mu.w1 = mu.w2 = parm$mu.w
    bet = parm$bet 
    B = parm$B
    DELTA = parm$DELTA
    
    SL_library = parm$SL_library
    SL_library_wr = parm$SL_library_wr
    x = sample( field.x, n, replace = T, prob = prob.x)#1:4*0.1  rep(0.25,4)
    w1 = rnorm(n, mu.w1, sig.w1)
    w2 = rnorm(n, mu.w2, sig.w2)

    prob.a = 0.5
    a = rbinom(n, 1, prob.a)
    
    sig.y = parm$sig.y
    
    if (distr.y=="gaussian"){
      y = rnorm(n, cbind(1, x, w1, w2, x*w1, x*a) %*%bet,sig.y)  
      FAMILY = gaussian() 
    }
    if (distr.y=="binomial"){
      y = rbinom(n, 1, plogis(cbind(1, x, w1, w2, x*w1, x*a) %*%bet))  
      FAMILY = binomial() 
    }
    
    dat = data.frame(y=y, x=factor(x),w1=w1, w2 = w2, a = a)
    
    alp = parm$alp
    Ca = qnorm((1+(1-alp)^(1/K))/2)
    
    #------- start -----
    
    f_h.hat = function(dat, emeth="nonpar", family=gaussian(), SL_library ){
      #----h.hat (w)
      #family = gaussian() or binomial()
      dat0 = dat1 = dat
      dat0$a = 0
      dat1$a = 1
      p = mean(dat$a)
      if(emeth=="true"){#--true h(w): 0. h(w)=0
        dat0$x = dat1$x =  as.numeric(as.character(dat$x))
        pred  =cbind(model.matrix(y~x+w1+w2+x:w1+x:a, data = dat0)%*%bet, 
                     model.matrix(y~x+w1+w2+x:w1+x:a, data = dat1)%*%bet)
        if (distr.y=="binomial") pred = plogis(pred)
      }
      if(emeth=="zero"){ #--estimate h(w): 1. h(w)=0
        pred = cbind(rep(0,n),rep(0,n))
      }
      if(emeth=="par"){#--estimate h(w): 2. parametric
        par_fit = glm(y~x+w1+w2+x:w1+x:a, data = dat, family = family)
        pred  =cbind(predict(par_fit, newdata=dat0),predict(par_fit, newdata=dat1))
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
    
    f_B = function(dat, emeth= "par", family=gaussian(), SL_library, Ca){
      
      h.hat = f_h.hat(dat, emeth, family,SL_library)
      p = mean(dat$a)
      D.hat= (dat$a/p-(1-dat$a)/(1-p))*(dat$y-h.hat)
      
      #--- delta(x) and sigma(x) 
      delta.hat = aggregate(D.hat, by = list(dat$x), FUN = mean)$x
      sigma.hat = aggregate(D.hat, by = list(dat$x), FUN = function(x){n=length(x);sd(x)*sqrt(n-1)/n})$x
      
      DL = (delta.hat - Ca*sigma.hat)
      DU = (delta.hat + Ca*sigma.hat)
      
      BL = (DL>0)
      BU = (DU>0)
      
      z = (delta.hat-DELTA)/sigma.hat
      max.D = (max(abs(z))<Ca)
      
      return(data.frame(BL, BU, DL, DU, max.D, z)) #
    }
    
    #-- run
    B1 = f_B(dat, emeth= "zero", family=FAMILY, SL_library, Ca) #1. h(w)=0
    B2 = f_B(dat, emeth= "par", family=FAMILY, SL_library, Ca) #2a. parametric and correct
    B3 = f_B(dat, emeth= "nonpar", family=FAMILY, SL_library, Ca) #3a. nonparametric without CV for sigma.hat
    B0 = f_B(dat, emeth= "true", family=FAMILY, SL_library, Ca) #1. h(w)=0

    cb = list(true = B0, zero = B1, par = B2, nonpar = B3)

    f.cp = function(x, B, DELTA, prob.x){
      BL = x[,1]; BU=x[,2]
      cp = ifelse(all(B[BL],BU[B]),1,0)
      scp.L = ifelse(all(B[BL]),1,0)
      scp.U = ifelse(all(BU[B]),1,0)
      size = sum(prob.x[BU&(!BL)])
  
      return(c(cp, scp.L, scp.U, size)) 
    }
    
    res = sapply(cb, f.cp, B=B, DELTA=DELTA, prob.x =prob.x)
    return(c(res))
  }

  set.seed(seed)

  n=nsample
  K = parm$K
  
  res = f.sim.discr(n, parm)
  return(res)
}


#-------seting------
nsample= 200

{
  parm = list(  
    distr.y = "gaussian", #"binomial",
    K = 4, x.mean = 1, x.sd = 2,
    mu.w=0, sig.w = 0.5, sig.y=0.5,
    bet = c(0, 0.2, 0.2, -0.2, 0.1, 0.2)*3, # cbind(1, x, w1, w2, x*w1, x*a) %*%bet
    alp = 0.1,
    SL_library = c("SL.lm","SL.glm", "SL.glm.interaction",  
                   "SL.randomForest")
  )
  
  distr.y = parm$distr.y
  K = parm$K
  prob.x = rep(1/K, K)
  x.mean = parm$x.mean
  x.sd = parm$x.sd
  field.x= qnorm(seq(0,1,length.out=K+2)[-c(1,K+2)],x.mean, x.sd)

  
  parm$prob.x = prob.x
  parm$field.x = field.x
  
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
  DELTA = sapply(field.x, f_delta.true)
  B = (DELTA>0)
  
  parm$DELTA = DELTA
  parm$B=B
  
}

#-------run------

print( f.pall(123, nsample= nsample, parm = parm) )


