# install all the required packages 
required_packages = c('dplyr','fda','ggplot2','Rsolnp','devtools')
(missed_packages = setdiff(required_packages,rownames(installed.packages())))

if(length(missed_packages)){
	sapply(missed_packages, install.packages)
} 

devtools::install_github("YunlongNie/DeFEC")

# load all the packages
library(dplyr);library(fda);library(ggplot2);library(Rsolnp)
library(DeFEC)
data(fat)

# observed contains a list with each element being the measurements for each subject, timepoints contains the corresponding time points
observed = fat$observed
timepoints = fat$timepoints

nbasis=8
spline_basis=create.bspline.basis(rangeval=c(2.1,11.1),breaks=quantile(summary(timepoints%>%do.call(c,.)),seq(0,1, len=nbasis-2)))
basis_intematlist = create_basis_intematlist(fat$timepoints, spline_basis)

# estimate the first three DeFECs 
DeFEC_fit = DeFEC(observed, timepoints, K=3,spline_basis, basis_intematlist ,thresh=1e-3, gammas=c(0,0.1,0.01), maxit=5e2)
# we can plot each DeFEC
DeFEC_fd = fd(do.call(cbind,DeFEC_fit), spline_basis)
fdaggD(DeFEC_fd)

# we can estimate the derivative function for each objects
DeFEC_pred = predict_DeFEC(betalist=DeFEC_fit,ylist=observed, tlist=timepoints,basis_intelist=basis_intematlist,spline_basis, nminus=2)

# we can plot the deriavtive estimation and the trajectory estimation for the first subject 
plot1 = DeFECplot(DeFEC_pred,DeFEC_fd,1)

plot1$fitplot
plot1$derivplot

