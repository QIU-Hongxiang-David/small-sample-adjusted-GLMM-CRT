library(magrittr)
library(dplyr)

expit<-binomial()$linkinv

gen.cluster.size<-function(n.cluster,mean.size,coef.var){
    #n.cluster: number of clusters
    #mean.size: mean cluster size
    #coef.var: coefficient of variation of cluster size
    if(coef.var==0){
        rep(mean.size,n.cluster)
    }else{
        var.size<-(mean.size*coef.var)^2
        
        #generate [cluster size]-5 and then add 5 back to ensure all cluster sizes>=1
        mean.size<-mean.size-5
        
        if(mean.size<var.size){
            cluster.sizes<-rnbinom(n.cluster,mu=mean.size,size=mean.size^2/(var.size-mean.size))
        }else if(mean.size==var.size){
            cluster.sizes<-rpois(n.cluster,mean.size)
        }else{
            #for mean.size>var.size, can only approximate with binomial distribution
            prob<-1-var.size/mean.size
            cluster.sizes<-rbinom(n.cluster,size=round(mean.size/prob),prob=prob)
        }
        
        sort(cluster.sizes)+5
    }
}

gen.cluster.covariate<-function(n) rbinom(n,1,.5)
gen.individual.covariate<-function(n) rnorm(n)

#true coefficients
#true intercept is 0
beta.cluster<-.8
beta.individual<-.7

calc.binary.random.effect.var<-function(ICC){
    ICC/(1-ICC)*pi^2/3
}

calc.Poisson.random.effect.var<-function(ICC,true.model,max.sigma2=20,...){
    #mean rate is E[Y]=e^(beta.individual^2/2+sigma^2/2)*E[e^(beta.cluster*cluster.X)] where sigma^2 is random effect variance
    #ICC=sigma^2/(sigma^2+log(1+1/E[Y]))
    #NEED TO UPDATE THIS FUNCTION ACCORDING TO DISTRIBUTION OF COVARIATES!!!
    
    #use uniroot to find sigma^2
    #max.sigma2: max sigma^2 in uniroot
    #...: other arguments to uniroot
    uniroot(function(sigma2){
        denominator<-exp(sigma2/2)
        if(true.model %in% c("B","D")){
            denominator<-denominator*exp(beta.individual^2/2)
        }
        if(true.model %in% c("C","D")){
            denominator<-denominator*.5*(1+exp(beta.cluster))
        }
        sigma2/(sigma2+log1p(1/denominator))-ICC
    },interval=c(0,max.sigma2),...)$root
}

generate.data<-function(n.cluster,mean.size,coef.var,ICC,outcome.type=c("binary","count"),true.model,max.sigma2=20,...){
    outcome.type<-match.arg(outcome.type)
    
    random.effect.var<-ifelse(outcome.type=="binary",
                              calc.binary.random.effect.var(ICC),
                              calc.Poisson.random.effect.var(ICC,true.model,max.sigma2=max.sigma2,...))
    
    cluster.sizes<-gen.cluster.size(n.cluster,mean.size,coef.var)
    
    n1<-floor(n.cluster/2)
    treatments1<-sample(c(rep(0:1,each=floor(n1/2)),sample(0:1,n1-floor(n1/2)*2)),size=n1)
    
    n2<-n.cluster-n1
    if(n1%%2==1 && n2%%2==1){
        treatments2<-sample(c(rep(0:1,each=floor(n2/2)),ifelse(sum(treatments1)%%2==1,0,1)),size=n2)
    }else{
        treatments2<-sample(c(rep(0:1,each=floor(n2/2)),sample(0:1,n2-floor(n2/2)*2)),size=n2)
    }
    
    treatments<-c(treatments1,treatments2)
    
    lapply(1:n.cluster,function(i){
        ni<-cluster.sizes[i]
        
        output<-data.frame(cluster.id=i,individual.id=1:ni,
                           cluster.X=gen.cluster.covariate(1),
                           individual.X=gen.individual.covariate(ni),
                           useless.cluster.X=gen.cluster.covariate(1),
                           useless.individual.X=gen.individual.covariate(ni),
                           treatment=treatments[i])%>%
            mutate(linear.part=rnorm(1,sd=sqrt(random.effect.var)))
        if(true.model %in% c("B","D")){
            output<-output%>%mutate(linear.part=linear.part+beta.individual*individual.X)
        }
        if(true.model %in% c("C","D")){
            output<-output%>%mutate(linear.part=linear.part+beta.cluster*cluster.X)
        }
        output<-output%>%mutate(Y=if(outcome.type=="binary"){
                       rbinom(ni,1,expit(linear.part))
                   }else{
                       rpois(ni,exp(linear.part))
                   })%>%
            select(-linear.part)
        
        output
    })%>%do.call(what=rbind)
}
