#######################################################################################
### dca is a function that calculates the net benefit associated with 
### various treatment strategies.
###
### see vickers and elkin (2006). Medical decision making Nov-Dec;26(6):565-74.
###
### The function has 8 inputs:
###	yvar:		binary response, nx1 column vector
###	xmatrix:	predictors of yvar, nxp matrix
###	xstart:	starting values for x-axis (threshold probability)
###			between 0 and 1; default is 0.01
###	xstop:	stopping values for x-axis (threshold probability)
###			between 0 and 1; default is 0.99
###	xby:		increment for threshold probability
###			default is 0.01
###	ymin:		minimum value for net benefit when plotting default decision curve
###			default is -0.05
###	ymax:		maximum value for net benefit when plotting default decision curve
###			default is 1.0
###	prob:		indication of whether the predictors are given as probabilities
###			px1 column vector, N for no, Y for yes
### 
### The function outputs a matrix with the following columns:
### 	threshold: 	threshold probablity
###			default is 1 to 100 in steps of 1
###	none:		net benefit of treating no patients
###	all:		net benefit of treating all patients
###	modelp1:	net benefit of treating patients according to the 1st predictor
###	modelp2:	net benefit of treating patients according to the 2nd predictor
###	... and so on for all p predictors specified
### 
### Examples:
###	dca(yvar=cancer, xmatrix=tpsa, prob="N")
###	dca(yvar=cancer, xmatrix=tpsa, xstart=0.10, xstop=0.30, prob="N")
###	dca(yvar=cancer, xmatrix=cbind(tpsa, fpsa), prob=c("N", "N"))
###	dca(yvar=cancer, xmatrix=cbind(tpsa, model1), prob=c("N", "Y"))
###
###
#######################################################################################


dca <- function (yvar, xmatrix, xstart=0.01, xstop=0.99, xby=0.01, ymin=-0.05, ymax=1.0, prob) 
{
# Check if the length of "prob" and xmatrix can match each other or not
  if (length(prob)==1) { if (!is.vector(xmatrix)) {stop("Lengths of PROB and XMATRIX do not match")}}
  else {if (length(prob)!=ncol(xmatrix)) {stop("Lengths of PROB and XMATRIX do not match")}}


# check that the inputs are correctly specified

	# the failure variable is coded as 0/1 or a probability
	if (length(yvar[yvar==0]) + length(yvar[yvar==1]) != length(yvar) ) {
		if ( max(yvar)<=1 && min(yvar)>=0 ) {
			print("Decision curve is computed by the predicted probabilities and without reference to actual outcome")
			}
		else {
			stop("failure variable must have values ranging from 0 to 1")
			}
		}


	# xstart is between 0 and 1
	if (xstart<0 | xstart>1) {
		stop("xstart must lie between 0 and 1")
		}

	# xstop is between 0 and 1
	if (xstop<0 | xstop>1) {
		stop("xstop must lie between 0 and 1")
		}
		
	# xby is between 0 and 1
	if (xby<=0 | xby>=1) {
		stop("xby must lie between 0 and 1")
		}

	# xstart is before xstop
	if (xstart>=xstop) {
		stop("xstop must be larger than xstart")
		}

	# prob is coded as "Y" or "N"
	for (i in 1:length(prob)) {
		if (prob[i]!="N" & prob[i]!="Y") {
			stop("prob input misspecified")
			}
		}
	# if prob is coded as Y, then the corresponding xvar is a probability
	if (length(prob)>1) { 
		for (i in 1:length(prob)) {
			if (prob[i]=="Y" & ( min(cbind(xmatrix)[,i])<0 | max(cbind(xmatrix)[,i])>1 )  ) {
				stop("prob input misspecified")
				}
			}
		}


# get rid of the missing values among predictors and response
  xy.matrix<-data.frame(yvar,xmatrix)
  for (i in 1:ncol(xy.matrix)) {
      xy.matrix<-xy.matrix[!is.na(xy.matrix[,i]),]
  }


# assign failure and predictor variables according to inputs
fail <- xy.matrix[,1]
nobs = nrow(xy.matrix)
print(paste("After deleting the missing, the Nobs of this analysis is",nobs,sep=" "))


# initialize the result matrix
result.out<-NULL


# Use loops to calculate the net benefit for each model for all threshold probabilities
for (i in 2:ncol(xy.matrix)) {
	pred<-xy.matrix[,i]
	model.n<-i-1

	# if the xvar is not a probability, then fit model and get predicted probability
	if (prob[model.n]=="Y") {
		p.n <- pred
		}
	else if (prob[model.n]=="N") {
		fit <- glm(fail ~ pred, family = binomial)
		p.n <- fit$fitted.values
		}

	# initialize variables that contain net benefit
	if(model.n==1) {
		all<-none<-NULL
		}
	modelp <- NULL

	# calculate net benefit for each threshold
	threshold<-seq(xstart,xstop,xby)
	for (j in 1:length(threshold)) {
		thres.val<-threshold[j]
		temp.mean <-mean(fail[p.n>=thres.val])
		temp.obs <-length(fail[p.n>=thres.val])
		tp <- temp.mean * temp.obs
		fp <- (1 - temp.mean) * temp.obs

		# append the net benefit to the output vector; 
		tempnb <- (tp - fp * thres.val / (1 - thres.val) ) / nobs 
		# if there are no tp or fp, then the net benefit = 0
		if (tempnb=="NaN") {
			tempnb <- 0
			}
		# if <= ymin, then blank out
		if (tempnb <= ymin) {
			modelp <- c(modelp, NA)
			}
		else {
			modelp <- c(modelp, tempnb)
			}

		# Calculate the "all" and "none" for one time since they should be constant
		px <- mean(fail)
		tempnb <- px - (1-px) * thres.val / (1 - thres.val)
		if (model.n==1) {
			# if all net benefit <= ymin, then blank out
			if (tempnb <= ymin) {
				all <- c(all, NA)
				}
			else {
				all <- c(all, tempnb )   
				}
			none <- c(none, 0)
			 }
	        }
	# append the results to the matrix
	result.out<-cbind(result.out,modelp)

     }
  
# put the threshold as a %
threshold <- threshold * 100

# save the results matrix as a data frame
result.out<-data.frame(as.data.frame(result.out),all,none,threshold)

# assign column names to the results matrix
names(result.out)<-c(paste("modelp",seq(1,model.n),sep=""),"all","none","threshold")

# plot the default decision curve

	# start by plotting the net benefit for treating all
	plot(result.out$threshold, result.out$none, type="l", lwd=2, xlim=c(xstart*100, xstop*100), ylim=c(ymin, ymax), xlab="Threshold probability (%)", ylab="Net benefit")

	# then add the net benefits for treat all; the default is to disregard net benefits < ymin
	nb <- result.out$all
	which.ymin <- which(nb > ymin)
	which.threshold <- result.out$threshold[which.ymin]
	which.nb <- nb[which.ymin]
	lines(which.threshold, which.nb, type="l", col=8, lwd=2)

	# initialize the legend label, color, and width using the standard specs of the none and all lines
	legendlabel <- c("None", "All")
	legendcolor <- c(17, 8)
	legendwidth <- c(2, 2)
	legendpattern <- c(1, 1)
	
	# then add the net benefits for treating according to each model; the default is to disregard net benefits < ymin
	for (i in 1:(ncol(xy.matrix) - 1)) {
		nb <- result.out[, i]
		which.ymin <- which(nb > ymin)
		which.threshold <- result.out$threshold[which.ymin]
		which.nb <- nb[which.ymin]
		lines(which.threshold, which.nb, type="l", col=i, lty=2)
		# add each model to the legend
		legendlabel <- c(legendlabel, paste("Model", i,sep=" "))
		legendcolor <- c(legendcolor, i)
		legendwidth <- c(legendwidth, 1)
		legendpattern <- c(legendpattern, 2)
		}

	# then add the legend
	legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)

# return the results matrix so that the results can be saved and used for plots
return(result.out)


}
