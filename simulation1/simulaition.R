#copy all library() to clusterEvalQ
library(magrittr)
library(lme4)
library(performance)
library(dplyr)

#need to create a folder cache/ to save temporary files in the working directory before running the simulation

run.all.methods<-function(data,outcome.type=c("binary","count"),true.model,alpha=.05){
    #data has following variables:
    #cluster.id
    #cluster.X
    #individual.X
    #useless.cluster.X
    #useless.individual.X
    #treatment
    #Y
    
    outcome.type<-match.arg(outcome.type)
    family<-if(outcome.type=="binary"){
        binomial()
    }else{
        poisson()
    }
    
    
    formula<-if(true.model=="A"){
        Y~treatment+(1|cluster.id)
    }else if(true.model=="B"){
        Y~treatment+individual.X+(1|cluster.id)
    }else if(true.model=="C"){
        Y~treatment+cluster.X+(1|cluster.id)
    }else if(true.model=="D"){
        Y~treatment+individual.X+cluster.X+(1|cluster.id)
    }
    
    output<-NULL
    
    
    #In the df below, BW and IO correspond to BW1 and BW2 in the paper, respectively
    
    
    #Model 1: correctly specified model 
    output<-rbind(output,tryCatch({
        df.residual<-nrow(data)-length(all.vars(formula))+1
        df.containment<-nrow(data)-n_distinct(data$cluster.id)-ifelse(true.model %in% c("B","D"),1,0)
        df.BW<-n_distinct(data$cluster.id)-length(all.vars(formula))+1
        df.IO<-df.BW+ifelse(true.model %in% c("B","D"),1,0)
        
        model<-lmer(formula,data=data)
        LMM.ICC<-tryCatch(icc(model)$ICC_adjusted,error=function(e) 0)
        
        model<-glmer(formula,data=data,family=family)
        # print(summary(model))
        
        #Wald t-test
        t.value<-coef(summary(model))["treatment","z value"]
        Wald.residual<-pt(abs(t.value),df.residual,lower.tail=FALSE)*2<=alpha
        Wald.containment<-pt(abs(t.value),df.containment,lower.tail=FALSE)*2<=alpha
        Wald.BW<-pt(abs(t.value),df.BW,lower.tail=FALSE)*2<=alpha
        Wald.IO<-pt(abs(t.value),df.IO,lower.tail=FALSE)*2<=alpha
        
        
        #LRT
        model.null<-update(model,~.-treatment)
        F.value<-anova(model,model.null)$Chisq[2]
        LRT.residual<-pf(F.value,1,df.residual,lower.tail=FALSE)<=alpha
        LRT.containment<-pf(F.value,1,df.containment,lower.tail=FALSE)<=alpha
        LRT.BW<-pf(F.value,1,df.BW,lower.tail=FALSE)<=alpha
        LRT.IO<-pf(F.value,1,df.IO,lower.tail=FALSE)<=alpha
        
        
        data.frame(model="1",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=c(Wald.residual,Wald.containment,Wald.BW,Wald.IO,
                            LRT.residual,LRT.containment,LRT.BW,LRT.IO),
                   LMM.ICC=LMM.ICC)
    },error=function(e){
        data.frame(model="1",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=rep(NA,8),LMM.ICC=NA)
    }))
    
    
    
    
    
    
    #Model 2: adjust for a useless person-level covariate
    formula.model2<-update(formula,.~.+useless.individual.X)
    output<-rbind(output,tryCatch({
        df.residual<-nrow(data)-length(all.vars(formula.model2))+1
        df.containment<-nrow(data)-n_distinct(data$cluster.id)-ifelse(true.model %in% c("B","D"),2,0)
        df.BW<-n_distinct(data$cluster.id)-length(all.vars(formula.model2))+1
        df.IO<-df.BW+ifelse(true.model %in% c("B","D"),2,1)

        model<-lmer(formula.model2,data=data)
        LMM.ICC<-tryCatch(icc(model)$ICC_adjusted,error=function(e) 0)

        model<-glmer(formula.model2,data=data,family=family)
        # print(summary(model))

        #Wald t-test
        t.value<-coef(summary(model))["treatment","z value"]
        Wald.residual<-pt(abs(t.value),df.residual,lower.tail=FALSE)*2<=alpha
        Wald.containment<-pt(abs(t.value),df.containment,lower.tail=FALSE)*2<=alpha
        Wald.BW<-pt(abs(t.value),df.BW,lower.tail=FALSE)*2<=alpha
        Wald.IO<-pt(abs(t.value),df.IO,lower.tail=FALSE)*2<=alpha


        #LRT
        model.null<-update(model,~.-treatment)
        F.value<-anova(model,model.null)$Chisq[2]
        LRT.residual<-pf(F.value,1,df.residual,lower.tail=FALSE)<=alpha
        LRT.containment<-pf(F.value,1,df.containment,lower.tail=FALSE)<=alpha
        LRT.BW<-pf(F.value,1,df.BW,lower.tail=FALSE)<=alpha
        LRT.IO<-pf(F.value,1,df.IO,lower.tail=FALSE)<=alpha


        data.frame(model="2",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=c(Wald.residual,Wald.containment,Wald.BW,Wald.IO,
                            LRT.residual,LRT.containment,LRT.BW,LRT.IO),LMM.ICC=LMM.ICC)
    },error=function(e){
        data.frame(model="2",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=rep(NA,8),LMM.ICC=NA)
    }))






    #Model 3: adjust for a useless cluster-level covariate
    formula.model3<-update(formula,.~.+useless.cluster.X)
    output<-rbind(output,tryCatch({
        df.residual<-nrow(data)-length(all.vars(formula.model3))+1
        df.containment<-nrow(data)-n_distinct(data$cluster.id)-ifelse(true.model %in% c("B","D"),1,0)
        df.BW<-n_distinct(data$cluster.id)-length(all.vars(formula.model3))+1
        df.IO<-df.BW+ifelse(true.model %in% c("B","D"),1,0)

        model<-lmer(formula.model3,data=data)
        LMM.ICC<-tryCatch(icc(model)$ICC_adjusted,error=function(e) 0)

        model<-glmer(formula.model3,data=data,family=family)
        # print(summary(model))

        #Wald t-test
        t.value<-coef(summary(model))["treatment","z value"]
        Wald.residual<-pt(abs(t.value),df.residual,lower.tail=FALSE)*2<=alpha
        Wald.containment<-pt(abs(t.value),df.containment,lower.tail=FALSE)*2<=alpha
        Wald.BW<-pt(abs(t.value),df.BW,lower.tail=FALSE)*2<=alpha
        Wald.IO<-pt(abs(t.value),df.IO,lower.tail=FALSE)*2<=alpha


        #LRT
        model.null<-update(model,~.-treatment)
        F.value<-anova(model,model.null)$Chisq[2]
        LRT.residual<-pf(F.value,1,df.residual,lower.tail=FALSE)<=alpha
        LRT.containment<-pf(F.value,1,df.containment,lower.tail=FALSE)<=alpha
        LRT.BW<-pf(F.value,1,df.BW,lower.tail=FALSE)<=alpha
        LRT.IO<-pf(F.value,1,df.IO,lower.tail=FALSE)<=alpha


        data.frame(model="3",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=c(Wald.residual,Wald.containment,Wald.BW,Wald.IO,
                            LRT.residual,LRT.containment,LRT.BW,LRT.IO),LMM.ICC=LMM.ICC)
    },error=function(e){
        data.frame(model="3",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=rep(NA,8),LMM.ICC=NA)
    }))
    
    
    
    
    
    
    
    #Model 4: adjust for a useless person- and cluster-level covatiate
    formula.model4<-update(formula,.~.+useless.cluster.X+useless.individual.X)
    output<-rbind(output,tryCatch({
        df.residual<-nrow(data)-length(all.vars(formula.model4))+1
        df.containment<-nrow(data)-n_distinct(data$cluster.id)-ifelse(true.model %in% c("B","D"),2,0)
        df.BW<-n_distinct(data$cluster.id)-length(all.vars(formula.model4))+1
        df.IO<-df.BW+ifelse(true.model %in% c("B","D"),2,1)
        
        model<-lmer(formula.model4,data=data)
        LMM.ICC<-tryCatch(icc(model)$ICC_adjusted,error=function(e) 0)
        
        model<-glmer(formula.model4,data=data,family=family)
        # print(summary(model))
        
        #Wald t-test
        t.value<-coef(summary(model))["treatment","z value"]
        Wald.residual<-pt(abs(t.value),df.residual,lower.tail=FALSE)*2<=alpha
        Wald.containment<-pt(abs(t.value),df.containment,lower.tail=FALSE)*2<=alpha
        Wald.BW<-pt(abs(t.value),df.BW,lower.tail=FALSE)*2<=alpha
        Wald.IO<-pt(abs(t.value),df.IO,lower.tail=FALSE)*2<=alpha
        
        
        #LRT
        model.null<-update(model,~.-treatment)
        F.value<-anova(model,model.null)$Chisq[2]
        LRT.residual<-pf(F.value,1,df.residual,lower.tail=FALSE)<=alpha
        LRT.containment<-pf(F.value,1,df.containment,lower.tail=FALSE)<=alpha
        LRT.BW<-pf(F.value,1,df.BW,lower.tail=FALSE)<=alpha
        LRT.IO<-pf(F.value,1,df.IO,lower.tail=FALSE)<=alpha
        
        
        data.frame(model="4",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=c(Wald.residual,Wald.containment,Wald.BW,Wald.IO,
                            LRT.residual,LRT.containment,LRT.BW,LRT.IO),LMM.ICC=LMM.ICC)
    },error=function(e){
        data.frame(model="4",method=rep(c("Wald","LRT"),each=4),
                   ddf=rep(c("residual","containment","BW","IO"),times=2),
                   reject=rep(NA,8),LMM.ICC=NA)
    }))
    
    
    
    output
}






source("generate_data.R")
n_replicates<-5e3 #number of runs per scenario
n_replicates_per_job<-5e2 #number of runs per job
seed<-153719
simulation.settings<-expand.grid(n.cluster=c(10,20),mean.size=c(50,100),coef.var=c(0,.75,1.5),ICC=c(.001,.01,.05,.1,.2),outcome.type=c("binary","count"),true.model=c("A","B","C","D"),simid=1:(n_replicates/n_replicates_per_job),stringsAsFactors=FALSE)%>%
    arrange(n.cluster,mean.size,coef.var,outcome.type,ICC,true.model,simid)%>%
    mutate(seed=seed+1:nrow(.))

run.once<-function(i,...){
    if(file.exists(paste0("cache/sim",i,".csv"))){
        return(NULL)
    }
    set.seed(simulation.settings$seed[i])
    result<-lapply(1:n_replicates_per_job,function(dummy) {
        data<-generate.data(simulation.settings$n.cluster[i],simulation.settings$mean.size[i],simulation.settings$coef.var[i],simulation.settings$ICC[i],simulation.settings$outcome.type[i],simulation.settings$true.model[i],...)
        run.all.methods(data,simulation.settings$outcome.type[i],simulation.settings$true.model[i],...)
    })%>%do.call(what=rbind)
    result%<>%mutate(n.cluster=simulation.settings$n.cluster[i],mean.size=simulation.settings$mean.size[i],coef.var=simulation.settings$coef.var[i],ICC=simulation.settings$ICC[i],outcome.type=simulation.settings$outcome.type[i],true.model=simulation.settings$true.model[i])
    
    data.table::fwrite(result,paste0("cache/sim",i,".csv"))
}



library(parallel)
n_cores<-detectCores(logical=FALSE)-1 #=4-1 on the machine on which we ran the simulation
cl<-makeCluster(n_cores)
clusterExport(cl,varlist=ls())
clusterEvalQ(cl,{
    suppressMessages({
        library(magrittr)
        library(lme4)
        library(performance)
        library(dplyr)
    })
})

dummy<-parLapply(cl,1:nrow(simulation.settings),run.once)

stopCluster(cl)

simulation.results<-lapply(list.files("cache/"),function(file.name){
    read.csv(paste0("cache/",file.name))
})%>%do.call(what=bind_rows)

data.table::fwrite(simulation.results,"simulation_results.csv")

type.I.error<-simulation.results%>%
    group_by(model,method,ddf,n.cluster,mean.size,coef.var,ICC,outcome.type,true.model)%>%
    summarize(type.I.error=mean(reject,na.rm=TRUE))%>%ungroup

data.table::fwrite(type.I.error,"type.I.error.csv")


LMM.ICC<-simulation.results%>%
    group_by(model,n.cluster,mean.size,coef.var,ICC,outcome.type,true.model)%>%
    summarize(mean.LMM.ICC=mean(LMM.ICC,na.rm=TRUE))%>%ungroup
data.table::fwrite(LMM.ICC,"LMM.ICC.csv")
