library(tidyverse)
library(lubridate)
library(gridExtra)
library(pROC)
library(data.table)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(viridis)
library(purrr)
library(broom)
library(profvis)
library(furrr)
library(mice)
library(MASS)
library(cvAUC)
library(slider)
plan(multiprocess)

setwd('//your_path') # set working directory 


dat_joined <- readRDS("dat_afe.RDS")# load data file 

weather=readRDS("weather_weighted.rds") # weighted weather data 

weather=weather %>% arrange(site,date)


dat_joined=dat_joined %>% arrange(site,date)
dat_joined=dat_joined %>% mutate(any_breast_fed=factor(case_when((f4a_breastfed==1|f4a_breastfed==2)~1,TRUE~0))) # binary breastfeeding
dat_joined$dt=as.numeric(dat_joined$date)-13848 # set date starting to 0 
dat_joined=dat_joined %>% mutate(cont=case_when(site %in% c(1,2,3,4) ~ 1,
                                                TRUE ~ 2),
                                 no_split=1)

dat_joined$site=as.numeric(dat_joined$site)
dat_joined$index=1:dim(dat_joined)[1]


dat_joined=dat_joined %>% mutate(base_age = case_when(base_age <= 24 ~ base_age,
                                                  base_age>24 & base_age < 36 ~ as.integer(24),
                                                  base_age>=36 & base_age < 48 ~ as.integer(36),
                                                  base_age>=48 & base_age < 60 ~ as.integer(48))
) 

countries=c("The Gambia","Mali","Mozambique","Kenya","India","Bangladesh","Pakistan")



# kernel density function +.001 at the end avoids numerical issues from areas with very little density 

dens_rat=function(pred_train1,pred_train0,preds_test,adjust){
  d1=density(pred_train1,from=0,to=1,adjust=adjust)
  d0=density(pred_train0,from=0,to=1,adjust=adjust)
  ind1=sapply(preds_test,function(x) which.min(abs(d1$x-x)))
  ind0=sapply(preds_test,function(x) which.min(abs(d0$x-x)))
  return(log(d1$y[ind1]+.001)-log(d0$y[ind0]+.001))
}

# 2d kernel density 

dens_rat2d=function(clim_pred_train1,clim_pred_train0,clin_pred_train1,clin_pred_train0,clim_preds_test,clin_preds_test){
  d1=kde2d(clim_pred_train1,clin_pred_train1,n=512,lims=c(0,1,0,1))
  d0=kde2d(clim_pred_train0,clin_pred_train0,n=512,lims=c(0,1,0,1))
  ind1x=sapply(clim_preds_test,function(x) which.min(abs(d1$x-x)))
  ind0x=sapply(clim_preds_test,function(x) which.min(abs(d0$x-x)))
  ind1y=sapply(clin_preds_test,function(y) which.min(abs(d1$y-y)))
  ind0y=sapply(clin_preds_test,function(y) which.min(abs(d0$y-y)))
  return(log(d1$z[cbind(ind1x,ind1y)]+.001)-log(d0$z[cbind(ind0x,ind0y)]+.001))
}


# createst quantiles for testing and training and keeps track of indices 
Create_TrnTst=function(dat_joined,p,percent=.80,clin_band,clim_band,split_by){

est_cutoff=dat_joined %>% group_by(dat_joined[[split_by]]) %>% summarize(cutoff=quantile(dt,percent))
train=bind_rows(map2(dat_joined %>% split(.[[split_by]]), est_cutoff$cutoff,function(x,y) x %>% filter(dt<y)))
test=dat_joined %>% filter(!(index %in% train$index))

train_ind=split(train$index,train[split_by]) 
test_ind=split(test$index,test[split_by]) 
orig_ind=train %>% split(.[split_by]) %>% purrr::map(. %>% group_by(dt,site) %>% dplyr::select(dt,site,index,viral_only) %>% nest(index,viral_only)) 
return(list(train_ind,test_ind,orig_ind))
}

# this function fits all models 
Post_Test_Odds=function(dat_joined=dat_joined,split_by="site",train_ind,test_ind,orig_ind,clin_band,clim_band,adjust=1,choices,
                        sens=.5,spec=.5){

bin_cols=c('f4a_drh_vomit','f4a_drh_blood','any_breast_fed')
num_cols=c('base_age','f4b_muac','viral_only_ma')   
names=c('base_age','f4a_drh_vomit','f4b_muac','f4a_drh_blood','any_breast_fed');names=c(names,"f4b_eyes","f4b_resp","wealth_index","f4b_temp","f4a_ppl_house")

pt_list=names[c(1:5)] # we only are using the first five variables of interest 

### Weather aggreation 

weather=bind_rows(weather %>% arrange(site,dt) %>% split(.$site) %>% purrr::map(. %>% tidyr::complete(dt=seq(min(dt),max(dt)))))
weather$site=as.numeric(weather$site)
weather$temp.aggr=unlist(weather %>% split(.$site) %>% future_map(function(x) x %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% slide_index_dbl(.x=x$temp,.i=x$dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1)))
weather$rain.aggr=unlist(weather %>% split(.$site) %>% future_map(function(x) x %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% slide_index_dbl(.x=x$rain,.i=x$dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1)))

#weather = weather %>% replace(is.na(.), 0)

dat_joined=dat_joined %>% left_join(weather,by=c("site","dt"))#%>% left_join(train_clin,by=c("site","dt")) 

#handling some missing weather by using last available, only a few days 
dat_joined$temp.aggr[1502:1505]=dat_joined$temp.aggr[1501]
dat_joined$rain.aggr[1502:1505]=dat_joined$rain.aggr[1501]
dat_joined$temp.aggr[2344:2348]=dat_joined$temp.aggr[2343]
dat_joined$rain.aggr[2344:2348]=dat_joined$rain.aggr[2343]

# creating seasonal curves 
dat_joined=dat_joined %>% mutate(Seasonal_sine=sin(2*pi*.$dt/365.25),
Seasonal_cosine=cos(2*pi*.$dt/365.25))

# Models are fit by site so we make lists by site of each training and testing dataset 
train_joint=map2(dat_joined %>% split(.[[split_by]]),train_ind,~.x %>% filter(index %in% .y))
test_joint=map2(dat_joined %>% split(.[[split_by]]),test_ind,~.x %>% filter(index %in% .y))


###Get indices of viral==1 
vir_ind_train=lapply(train_joint,function(x) which(x$viral_only==1))
vir_ind_test=lapply(test_joint,function(x) which(x$viral_only==1))

train_dat=dat_joined %>% filter(index %in% unlist(train_ind)) %>% arrange_at(vars(split_by,'index'))
test_dat=dat_joined %>% filter(index %in% unlist(test_ind)) %>% arrange_at(vars(split_by,'index'))

##
## This let's us only fit the models we choose 
if(any(choices[,3])){

  
mods=train_joint %>% purrr::map(~glm(viral_only~temp.aggr+rain.aggr,data=.,family="binomial"))


vir_ind_train=lapply(train_joint,function(x) which(x$viral_only==1))
vir_ind_test=lapply(test_joint,function(x) which(x$viral_only==1))


train_preds_clim=map2(mods,train_joint,~predict(.x,newdata=.y,type="response"))
test_preds_clim= map2(mods,test_joint,~predict(.x,newdata=.y,type="response"))
OR_Clim=pmap(list(train_preds_clim,vir_ind_train,test_preds_clim), function(trn,ind,tst) dens_rat(trn[ind],trn[-ind],tst,adjust))

test_dat$OR_Clim=unlist(OR_Clim)

}else{test_dat$OR_Clim=0}


### We stopped using this model but left the details here 
#if(any(choices[,4])){
#mods2=train_joint %>% purrr::map(~glm(as.formula(paste('viral_only~',paste(c(bin_cols,num_cols),"_p",collapse='+',sep=""),sep="")),data=.,family="binomial"))
#train_preds_clin=map2(mods2,train_joint,~predict(.x,newdata=.y,type="response"))
#test_preds_clin= map2(mods2,test_joint,~predict(.x,newdata=.y,type="response"))
#OR_Clin=pmap(list(train_preds_clin,vir_ind_train,test_preds_clin), function(trn,ind,tst) dens_rat(trn[ind],trn[-ind],tst,adjust))
#test_dat$OR_Clin=unlist(OR_Clin)
#} else{OR_Clin=0}

if(any(choices[,5])){
  mods3=train_joint %>% purrr::map(~glm(viral_only~Seasonal_sine+Seasonal_cosine,data=.,family="binomial"))
  train_preds_Seas=map2(mods3,train_joint,~predict(.x,newdata=.y,type="response"))
  test_preds_Seas= map2(mods3,test_joint,~predict(.x,newdata=.y,type="response"))
  OR_Seas=pmap(list(train_preds_Seas,vir_ind_train,test_preds_Seas), function(trn,ind,tst) dens_rat(trn[ind],trn[-ind],tst,adjust))
  test_dat$OR_Seas=unlist(OR_Seas)
}else{test_dat$OR_Seas=0}#OR_Clim_Clin=1}
#mods2=train_joint %>% purrr::map(~glm(viral_only~Seasonal_sine+Seasonal_cosine+temp.aggr+rain.aggr,data=.,family="binomial"))

if(any(choices[,6])){
  OR_Clim_Clin=pmap(list(train_preds_clim,train_preds_Seas,vir_ind_train,test_preds_clim,test_preds_Seas), function(v,w,x,y,z) dens_rat2d(v[x],v[-x],w[x],w[-x],y,z))
  test_dat$OR_Clim_Clin=unlist(OR_Clim_Clin)
}else{test_dat$OR_Clim_Clin=0}


#

###Additional diagnostic code - this does a bootstrap so the analysis at the end has to be done slightly differently 
if(any(choices[,7])){
  
  addl_diag=
    1:500 %>% purrr::map(function(x) test_dat %>% dplyr::select(caseid,site,viral_only) %>% sample_n(dim(test_dat)[1],replace=T)) %>% 
    purrr::map(function (d) d %>% mutate(test_post=case_when(viral_only==1 ~ rbinom(n(),1,sens),
                                                             TRUE ~ as.integer(1)-rbinom(n(),1,spec)),
                                         OR_Vir=case_when(test_post==1 ~ log((sens/(1-spec))),TRUE~log(((1-sens)/spec)))))
  

}else{
addl_diag=test_dat %>% dplyr::select(caseid,site,viral_only)}#OR_Clim_Clin=1}


### Current patient model - we always fit it but it can still be excluded later 
out=glm(as.formula(paste('viral_only~',paste(pt_list,collapse='+'),sep="")),family="binomial",data=train_dat)


preds_train=predict(out,type="response")
preds_test=predict(out,test_dat,type="response")
test_dat$pred=preds_test

pred1=preds_train[which(train_dat$viral_only==1)]
pred0=preds_train[which(train_dat$viral_only==0)]

OR=dens_rat(pred1,pred0,preds_test,adjust)
test_dat$OR=OR


### this function is used for pre-test odds for new sites 
# weighted_avg_fun=function(x,y) {
#   #wts=rep(1,clin_band)
#   #wts=1/(1:3)
#   wts=Wendland(1:clin_band,theta=clin_band,k=1,dimension=1)
#   if(length(y)==1) {NA} else {
#     today=y[length(y)]
#     inds=today-y
#     x=x[-length(x)]
#     weighted.mean(x,w=wts[inds],na.rm=T)}
# }


# using historical "training" data 
tst_daymo=train_dat %>% split(.$site) %>% purrr::map(. %>% mutate(daymo=(as.numeric(ymd(.$date.x)) %% 365)) %>% arrange(daymo)) 
tst2=purrr::map(tst_daymo,function(x) 0:364 %>% purrr::map(function(z) x %>% filter(daymo>z-30,daymo<z+30) %>% summarize(mean(viral_only))))
priordf=data.frame(site=rep(1:7,each=365),seas_ma=unlist(tst2),daymo=rep(1:365,7)) %>% replace(is.na(.), 0) %>% mutate(preOdds=log(seas_ma/(1-seas_ma)))

test_dat$daymo=as.numeric(ymd(test_dat$date.x)) %% 365
test_dat=test_dat  %>% left_join(priordf,by=c("site","daymo"))




return(list(test_dat,addl_diag))
}

#Calculates OR predictions 
Calc_AUC=function(test_dat,Choice){
  pre=Choice[1]
  cur_clin=Choice[2]
  Clim=Choice[3]
  Clin=Choice[4]
  Seas=Choice[5]
  Clim_Clin=Choice[6]
  Vir=Choice[7]
  
  results=exp((if(pre==1) test_dat$preOdds else 0)+(if(cur_clin==1) test_dat$OR else 0)+(if(Clim==1) test_dat$OR_Clim else 0)+(if(Clin==1) test_dat$OR_Clin else 0)+
                (if(Seas==1) test_dat$OR_Seas else 0) + (if(Clim_Clin==1) test_dat$OR_Clim_Clin else 0) + (if(Vir==1) test_dat$OR_Vir else 0))#*(if(Hol==1) Hol else 1)
  test_dat$pred_OR=results
  return(test_dat)
}
#choices=matrix(c(c(1,1,0,0,0)),ncol=5,byrow=T)


formula=c("Pre-test","PresPtnt","Climate","PriorPtnt","Seasonal","Joint","Addl. Diagnostic")

#Pick which models to fit 
choices=matrix(c(c(1,1,0,0,0,0,0),
                 c(0,1,0,0,0,0,0),
                 c(1,1,1,0,0,0,0),
                 c(0,1,1,0,0,0,0),
                 c(1,1,0,0,0,1,0),
                 c(0,1,0,0,0,1,0),
                 c(1,1,0,0,1,0,0),
                 c(0,1,0,0,1,0,0),
                 c(1,1,1,0,1,0,0),
                 c(0,1,1,0,1,0,0)
),
ncol=7,byrow=T)

nms=purrr::map(lapply(1:nrow(choices), function(i) choices[i,]), ~paste(formula[which(.==1)],collapse=" * "))


# This is where you run the functions and get results 
cln_bnd=28
clm_bnd=14
set.seed(5)
result=data.frame(caseid=NA,abx=NA,site=NA,iter=NA,sens=NA,spec=NA,true=NA,choice=NA,clim_band=NA,clin_band=NA,pred_OR=NA)
#itera1=.8#
for (itera1 in .8){#seq(.7,.95,by=.05)){
  #for (each in 1:dim(sensspec)[1]){
  sens=.5#sensspec[each,1]
  spec=.5#sensspec[each,2]
  
  print(itera1)
  
  splby="site" 
  #print(clm_bnd)
  inds=Create_TrnTst(dat_joined,p,percent=itera1,cln_bnd,clm_bnd,splby)
  #for (cln_bnd in c(seq(35,70,by=7))){
  
  ORS=Post_Test_Odds(dat_joined=dat_joined,split_by=splby,inds[[1]],inds[[2]],inds[[3]],clin_band=cln_bnd,clim_band=clm_bnd,adjust=1,choices,sens,spec)  #mat[itera1,itera2]=as.numeric(roc(ORS[[1]]$viral_only,ORS[[7]])$auc)
  for (itera2 in 1:dim(choices)[1]){
    sv=Calc_AUC(ORS[[1]],Choice=choices[itera2,])
    result=rbind(result,data.frame(caseid=ORS[[1]]$caseid,abx=ORS[[1]]$antibiotic,site=ORS[[1]]$site,iter=itera1,sens=sens,spec=spec,true=ORS[[1]]$viral_only,choice=itera2,clim_band=clm_bnd,clin_band=cln_bnd,pred_OR=sv$pred_OR))
  }
}
#}
result=result[-1,]


#write.csv(result,"cvPTO_for_clinband-WTS.csv",row.names = F)
#result=read.csv("cvPTO_for_clinband.csv")

result
res=result  %>% group_by(iter,choice) %>% dplyr::summarize(AUC=as.numeric(roc(true,pred_OR,ci=T)$auc),
                                                           low=as.numeric(roc(true,pred_OR,ci=T)$ci[1]),
                                                           high=as.numeric(roc(true,pred_OR,ci=T)$ci[3]),
                                                           width=high-low)

res$form=unlist(nms[as.numeric(res$choice)])
order=(res %>% group_by(form) %>% summarize(AUC=mean(AUC)) %>% arrange((AUC)))$form

res$form=factor(res$form,levels=order)

ggplot(data=res,
       aes(x = AUC,y = form, xmin = low, xmax = high),xlim=c(.725,.925)) + 
  geom_point(aes(col=form), show.legend = FALSE) +
  geom_errorbarh(aes(col=form),cex=1, show.legend = FALSE)+  
  scale_color_viridis_d(begin=.15,end=.85) + facet_wrap(~iter,ncol=2) + theme_bw() +
  scale_x_continuous(breaks=seq(0.725,.925,by=.025),labels=seq(0.725,.925,by=.025)) + 
  theme(aspect.ratio = .7,text = element_text(size=14) ) + 
  ylab("Tests") + ggtitle("Forest Plots of Rolling-Origin-Recalibration Evaluation.")
