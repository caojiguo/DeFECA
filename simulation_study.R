required_packages = c('dplyr','fda','ggplot2','Rsolnp','devtools','caret')
(missed_packages = setdiff(required_packages,rownames(installed.packages())))

if(length(missed_packages)){
	sapply(missed_packages, install.packages)
} 

devtools::install_github("YunlongNie/DeFEC")



library(fda);library(dplyr);library(DeFEC);library(caret)


data(pcs_simulation)
pc1 = pcs_simulation[[1]]
pc2 = pcs_simulation[[2]]

# data generation 
seedset =1
sigma_sd  = 1
set.seed(seedset)
nsim = 200
nobs = sample(x=3:7,size=nsim, replace=TRUE)
s=(rgamma(n=nsim,1,0.03))
s = s-mean(s)
coefsy= matrix(rep(s, each =pc1$basis$nbasis), pc1$basis$nbasis)*(matrix(rep(coef(pc1), length(s)), length(s),byrow=TRUE)%>%t)

s2=(rgamma(n=nsim,1,0.1))
s2 = s2-mean(s2)
coefsy2= matrix(rep(s2, each =pc2$basis$nbasis), pc2$basis$nbasis)*(matrix(rep(coef(pc2), length(s)), length(s),byrow=TRUE)%>%t)
derivfds = fd(coefsy+coefsy2, pc2$basis)
initial_values = rnorm(nsim,sd=3)

(range_time = pc1$basis$rangeval)
timegrid =seq(range_time[1],range_time[2], by=0.05)
timepoints_index = lapply(1:nsim, function(x) {sample(nobs[x],x=1:length(timegrid),replace=F)%>%sort})
timepoints  = lapply(1:nsim, function(x){
	timegrid[timepoints_index[[x]]]
})

data(basis_intemat_simulation)
basis_intematlist = lapply(1:length(timepoints), function(x){
	basis_intemat[timepoints_index[[x]],]%>%t
})

observed= lapply(1:length(timepoints),function(i){
 colSums(rep((derivfds%>%coef)[,i],each=length(timepoints[[i]]))%>%matrix(nrow=pc2$basis$nbasis,byrow = T)*basis_intematlist[[i]])+initial_values[i]+rnorm(length(timepoints[[i]]), sd=sigma_sd)
})
 



# estimating the DeFEC

nbasis=12
data(basis_intemat12)
spline_basis =create.bspline.basis(rangeval=c(0.2,11.1),nbasis=nbasis)


gamma_pool = c(0,1e-3,1e-1,1,10)

# do the cv to select the gamma for the first PC
set.seed(1010)
kfolds=10
library(DeFEC)

# we use a 10 fold CV to select the gamma for each DeFEC
folds = caret::createFolds(1:length(observed),k=kfolds)
cv_res1 = select_gamma(evalc=5, gamma_pool, folds, obs=observed, timep=timepoints,basis_intelist=basis_intematlist,spline_basis=spline_basis,thres=1e-2)

(gamma_select = cv_res1$gamma_select)


pc1_multi= FPCfirst_multi_start(evaltimes=5,obs=observed, timep=timepoints,basis_intelist =basis_intematlist, spline_basis=spline_basis, gamma=gamma_select,threshold=1e-2,minit=20,maxeval=1e3)
pc1_fit = pc1_multi[[which.min(sapply(pc1_multi, function(x) x$options$value))]]
previous_beta=list()
previous_beta[[1]]=pc1_fit$beta

fit_pc1 = predict_DeFEC(betalist=previous_beta,ylist=observed, tlist=timepoints,basis_intelist=basis_intematlist,spline_basis=spline_basis, nminus=1)
residuals1 = fit_pc1$residuals
previous_beta=list()
previous_beta[[1]]=pc1_fit$beta


cv_res2 = select_gamma(evalc=5, gamma_pool, folds, obs=residuals1, timep=timepoints,basis_intelist=basis_intematlist,spline_basis=spline_basis,thres=1e-2)

cv_res2$cv_res%>%rowMeans
(gamma_select=cv_res2$gamma_select)

pc2_multi= FPCfirst_multi_start(evaltimes=5,obs=residuals1, timep=timepoints,basis_intelist =basis_intematlist, spline_basis=spline_basis, gamma=gamma_select,threshold=1e-2,minit=20,maxeval=1e3)
pc2_fit = pc2_multi[[which.min(sapply(pc2_multi, function(x) x$options$value))]]


previous_beta=list()
previous_beta[[1]]=pc1_fit$beta
previous_beta[[2]]=pc2_fit$beta

DeFEC_fd = fd(do.call(cbind,previous_beta), spline_basis)
fdaggD(DeFEC_fd)

DeFEC_pred = predict_DeFEC(betalist=previous_beta,ylist=observed, tlist=timepoints,basis_intelist=basis_intematlist,spline_basis, nminus=2)

# we can plot the estimated derivative function with the true one 
plot(DeFEC_pred$yfd_fit[2])
lines(derivfds[2],col=2)

# we can plot the fitted value with the true one 
plot1 = DeFECplot(DeFEC_pred,DeFEC_fd,1:10)
plot1$fitplot



