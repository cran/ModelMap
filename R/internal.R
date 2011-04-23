#############################################################################################
#############################################################################################
###################################### Get Rasts ############################################
#############################################################################################
#############################################################################################

getRasts<-function(){
	## This function prompts user to browse to and select raster layers used for model predictors.
	## Returns vector of rasters.  

	## Adds to file filters to Cran R Filters table.
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))

	userPrompt <- "Yes"
      rastnmVector <- {}
      while(userPrompt == "Yes"){
		## Gets raster format from user.
      	rasterType = select.list(c("Imagine Image", "ArcInfo Grid"), title="Select type of raster.")
		if(rasterType==""){stop("Type of Raster must be selected")}
		if(rasterType == "Imagine Image"){	
			rasts <- choose.files(caption="Select image", filters = Filters["img",], multi = TRUE)
				if(is.null(rasts)){
					stop("")}
			}
		if(rasterType == "ArcInfo Grid"){
			rasts = choose.dir(default=getwd(),caption="Select grid")
			if(is.null(rasts)){
				stop("")}
		}
		if(is.null(rasterType)){
			stop("")
		}

		## Prompts user to select more.
		userPrompt = select.list(c("Yes", "No"), title="Another?")
		if(userPrompt==""){userPrompt<-"No"}
		## Compiles list of rasters.
		rastnmVector = c(rastnmVector, rasts)
	}
	return(rastnmVector)
}



#############################################################################################
###################################### Diagnostics ##########################################
#############################################################################################


diagnostics.function<-function(	model.obj=NULL,
						model.type,
						predList,
						PRED,
						MODELfn=NULL,
						MODELpredfn=NULL,
						response.name,
						response.type,
						folder=getwd(),
						device.type="default",
						jpeg.res=72,
						device.width=7,
						device.height=7,
						main=basename(MODELpredfn),
						cex=par$cex,
						req.sens,
						req.spec,
						FPC,
						FNC	){


### check predictions for -9999 ###

PRED.NA<-PRED[PRED$pred==-9999,]
Npred.NA<-nrow(PRED.NA)
PRED<-PRED[PRED$pred!=-9999,]
if(Npred.NA>0){
	warning(Npred.NA," datapoints not included in diagnostic plots because categorical factored predictor contained levels not found in training data")}


### Write out tables ###

if(response.type == "binary"){
	OPTTHRESHfn<-paste(MODELpredfn,"_optthresholds.csv",sep="")
	PREDPREVfn<-paste(MODELpredfn,"_prevalence.csv",sep="")

	opt.thresh<-error.threshold.plot(	PRED,opt.methods=optimal.thresholds(),plot.it=FALSE,
							req.sens=req.sens,req.spec=req.spec,FPC=FPC,FNC=FNC)

	pred.prev<-predicted.prevalence(PRED, threshold = opt.thresh$threshold)
	pred.prev<-cbind(opt.thresh$opt.methods, pred.prev)

	write.table(	opt.thresh,file=OPTTHRESHfn,sep=",",row.names=FALSE)
	write.table(	pred.prev,file=PREDPREVfn,sep=",",row.names=FALSE)
}


if(response.type == "continuous"){

	COR.pearson<-cor(PRED$pred,PRED$obs,method="pearson")
	COR.pearson<-round(COR.pearson,2)
	
	COR.spearman<-cor(PRED$pred,PRED$obs,method="spearman")
	COR.spearman<-round(COR.spearman,2)

	Resid<-PRED$obs-PRED$pred
	n<-length(Resid)
	MSE<-mean(Resid^2)
	RMSD<-(sum((Resid^2))/(n-1))^.5
	RMSD<-round(RMSD,2)

	lm.pred<-lm(obs~pred,data=PRED)
	b<-round(lm.pred$coefficients[1],2)
	m<-round(lm.pred$coefficients[2],2)
}

### make graphs ###

if(!"none"%in%device.type){
for(i in 1:length(device.type)){

	### Output filenames ###

	if(device.type[i] == "default"){
		SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot",sep="")
		IMPORTANCEfn<-paste(MODELpredfn,"_importance",sep="")
		THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots",sep="")
	}

	if(device.type[i] == "jpeg"){
		SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot.jpg",sep="")
		IMPORTANCEfn<-paste(MODELpredfn,"_importance.jpg",sep="")
		THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots.jpg",sep="")
	}

	if(device.type[i] == "pdf"){
		SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot.pdf",sep="")
		IMPORTANCEfn<-paste(MODELpredfn,"_importance.pdf",sep="")
		THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots.pdf",sep="")
	}

	if(device.type[i] == "postscript"){
		SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot.ps",sep="")
		IMPORTANCEfn<-paste(MODELpredfn,"_importance.ps",sep="")
		THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots.ps",sep="")
	}

	if(device.type[i] == "win.metafile"){
		SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot.emf",sep="")
		IMPORTANCEfn<-paste(MODELpredfn,"_importance.emf",sep="")
		THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots.emf",sep="")
	}


	###################   model.obj Based   ##############################


	if(model.type=="RF"){

		if(i==1){library(randomForest)}

		print("IMPORTANCEfn:")
		print(IMPORTANCEfn)

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)
							print("Importance Device")}
		if(device.type[i]=="jpeg"){jpeg(filename=IMPORTANCEfn,width = device.width, height = device.height, res=jpeg.res, unit="in")}
		if(device.type[i]=="postscript"){postscript(file=IMPORTANCEfn,width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=IMPORTANCEfn,width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=IMPORTANCEfn,width = device.width, height = device.height, 
								pointsize = 12,restoreConsole = TRUE)}
 
		opar<-par(cex=cex)
		varImpPlot(model.obj,main="Relative Influence",cex=cex)
		mtext(main,side=3,line=-4,cex=1.3*cex,outer=TRUE)
		par(opar)

		if(device.type[i]!="default"){dev.off()}
	}

	if(model.type=="SGB"){

		if(i==1){library(gbm)}
	
		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=IMPORTANCEfn,width = device.width, height = device.height, res=jpeg.res, unit="in")}
		if(device.type[i]=="postscript"){postscript(file=IMPORTANCEfn,width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=IMPORTANCEfn,width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=IMPORTANCEfn,width = device.width, height = device.height, 
								pointsize = 12,restoreConsole = TRUE)}

		opar<-par(las=1,mar=(c(5, 11, 4, 2) + 0.1),cex=cex)
		summary(model.obj)
		par(las=0)
		mtext("Relative Influence",side=3,line=.7,cex=1.5*cex)
		mtext(main,side=3,line=2.7,cex=1.5*cex)
		mtext("Predictors",side=2,line=10,cex=1*cex)
		par(opar)

		if(device.type[i]!="windows"){dev.off()}
	}


	###################   PRED Based   ##############################

	### binary ###

	if(response.type == "binary"){
	
		if(i==1){library("PresenceAbsence")}

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=THRESHOLDPLOTSfn,width = device.width, height = device.height, res=jpeg.res, unit="in")}
		if(device.type[i]=="postscript"){postscript(file=THRESHOLDPLOTSfn,width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=THRESHOLDPLOTSfn,width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=THRESHOLDPLOTSfn,width = device.width, height = device.height, 
							pointsize = 12,restoreConsole = TRUE)}

		opar<-par(cex=cex)
		presence.absence.summary(PRED,main=main,legend.cex=cex,opt.legend.cex=cex)
		par(opar)

		if(device.type[i]!="default"){dev.off()}	
	}

### Continuous ###

	if(response.type == "continuous"){

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)
							print("Scatterplot Device")}
		if(device.type[i]=="jpeg"){jpeg(filename=SCATTERPLOTfn,width = device.width, height = device.height, res=jpeg.res, unit="in")}
		if(device.type[i]=="postscript"){postscript(file=SCATTERPLOTfn,width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=SCATTERPLOTfn,width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=SCATTERPLOTfn,width = device.width, height = device.height, 
										pointsize = 12,restoreConsole = TRUE)}

		opar<-par(pty="s",cex=cex)
		lim<-range(PRED$obs,PRED$pred)
		plot(PRED$pred,PRED$obs,xlab="predicted",ylab="observed",xlim=lim,ylim=lim,main="")
		abline(a=0,b=1,lty=2)
		abline(lm.pred)

		mtext(main,side=3,line=1.5,cex=1.3*cex)

		mtext(paste("RMSD:",RMSD," ",sep=""),side=1,line=-4.5,adj=1,cex=.8*cex)	
		mtext(paste("pearson's cor: ",COR.pearson," ",sep=""),side=1,line=-3.5,adj=1,cex=.8*cex)
		mtext(paste("spearman's cor: ",COR.spearman," ",sep=""),side=1,line=-2.5,adj=1,cex=.8*cex)
		mtext(paste("obs = ",m,"(pred) + ",b," ",sep=""),side=1,line=-1.5,adj=1,cex=.8*cex)

		par(opar)

		if(device.type[i]!="default"){dev.off()}
	}
}
}
}

#############################################################################################
###################################### Get Rasts ############################################
#############################################################################################

getRasts<-function(){
	## This function prompts user to browse to and select raster layers used for model predictors.
	## Returns vector of rasters.  

	## Adds to file filters to Cran R Filters table.
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))

	userPrompt <- "Yes"
      rastnmVector <- {}
      while(userPrompt == "Yes"){
		## Gets raster format from user.
      	rasterType = select.list(c("Imagine Image", "ArcInfo Grid"), title="Select type of raster.")
		if(rasterType==""){stop("Type of Raster must be selected")}
		if(rasterType == "Imagine Image"){	
			rasts <- choose.files(caption="Select image", filters = Filters["img",], multi = TRUE)
				if(is.null(rasts)){
					stop("")}
			}
		if(rasterType == "ArcInfo Grid"){
			rasts = choose.dir(caption="Select grid")
			if(is.null(rasts)){
				stop("")}
		}
		if(is.null(rasterType)){
			stop("")
		}

		## Prompts user to select more.
		userPrompt = select.list(c("Yes", "No"), title="Another?")
		if(userPrompt==""){userPrompt<-"No"}
		## Compiles list of rasters.
		rastnmVector = c(rastnmVector, rasts)
	}
	return(rastnmVector)
}


#############################################################################################
################################ SGB - Model Creation #######################################
#############################################################################################

model.SGB<-function(	qdata,
				predList,
				response.name,
				response.type,
				seed=NULL,
				n.trees=NULL,                 # number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                 	 	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				train.fraction = 1.0,       	# fraction of data for training,
                 	 	n.minobsinnode = 10,         	# minimum total weight needed in each node
				keep.data=TRUE,
				var.monotone = NULL
){


## This function generates a model using gbm.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: SGB model

if(!is.null(seed)){
	set.seed(seed)}

if(response.type=="binary"){distribution="bernoulli"}
if(response.type=="continuous"){distribution="gaussian"}

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.y<-qdata[,response.name]
if(response.type=="binary"){qdata.y[qdata.y>0]<-1}

if(is.null(n.trees)){

	if(keep.data==FALSE){
		warning("keep.data reset to TRUE because required by gbm.more() function needed for OOB determination of optimal number of trees.")}


	SGB <- gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=100,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction,          	
				train.fraction = train.fraction,       	           		
				n.minobsinnode = n.minobsinnode,
				keep.data=TRUE,
				var.monotone=var.monotone)

	# check performance using an out-of-bag estimator
	best.iter <- suppressWarnings(gbm.perf(SGB,method="OOB",plot.it=FALSE))
      
	# iterate until a sufficient number of trees are fit

	while(SGB$n.trees - best.iter < 10){
     	 	# do 100 more iterations
      	SGB <- gbm.more(SGB,100)          
      	best.iter <- suppressWarnings(gbm.perf(SGB,method="OOB",plot.it=FALSE))
	}
	SGB$best.iter <- best.iter
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package. OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap.")

}else{
	SGB <- gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=n.trees,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction,          	
				train.fraction = train.fraction,       	           		
				n.minobsinnode = n.minobsinnode,
				keep.data=keep.data,
				var.monotone=var.monotone)
	if(train.fraction<1){
		SGB$best.iter <- suppressWarnings(gbm.perf(SGB,method="test",plot.it=FALSE))
		if(SGB$best.iter>0.9*n.trees){
			warning("Best number of trees is ", SGB$best.iter, " and total number trees tested was ", n.trees, ". You may want to explore increasing the 'n.trees' argument.")
		}
	}

}

SGB$response<-response.name

return(SGB)

}




#############################################################################################
################################### SGB - Predict ###########################################
#############################################################################################

prediction.SGB<-function(	prediction.type,
					qdata,
					train,
					response.name=deparse(substitute(SGB$response.name)),
					SGB,
					na.action="na.omit",
					n.trees
					){

## This function makes predictions for SGB model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

response.type<-switch(SGB$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")

if(response.type=="unknown"){
	stop("supplied model.obj has an unknown response type")
}

if(is.null(response.name)){
	stop("must provide response name")}

predList<-SGB$var.names

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.y<-qdata[,response.name]

if(response.type=="binary"){qdata.y[qdata.y>0]<-1}

if(prediction.type=="TRAIN"){

	if(!is.null(SGB$levels)){
		print("assigning levels")
		print(SGB$levels)
		for(p in names(SGB$levels)){
			qdata.x[,p]<-factor(qdata.x[,p],levels=SGB$levels[[p]])
		}
	}

	pred<-predict.gbm(	object=SGB,
					newdata=qdata.x,
					n.trees=n.trees,
					type="response",
					single.tree=FALSE)

	#print(cbind(qdata.x,pred))
}

if(prediction.type=="TEST"){


	## check -9999 

	not9<-!apply(qdata.x[,]==-9999,1,any)
	qdata.x.not9<-qdata.x[not9,]
	qdata.x.not9<-data.frame(qdata.x.not9)
	qdata.x.not9<-qdata.x.not9[,match(SGB$var.names,names(qdata.x.not9))]

	## checking for NA's
	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("predictor",names(qdata.x)[p],"contains NA's"))
		}
		if(na.action=="na.roughfix"){
			warning("Replacing NA predictors with median value or most common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with NA predictors")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
			}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}

	## Factors

	if(!is.null(SGB$levels)){
		missing.levels<-SGB$levels
		for(p in names(SGB$levels)){
			missing.levels[[p]]<-unique(qdata.x.not9[,p][!qdata.x.not9[,p]%in%SGB$levels[[p]]])
			qdata.x.not9[,p]<-factor(qdata.x.not9[,p],levels=SGB$levels[[p]])
		}
	}

	## checking for missing factor categories
	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("categorical factored predictor",names(qdata.x)[p],"contains  levels",paste(missing.levels[[names(qdata.x)[p]]],collapse=", "), "not found in training data"))
		}

		if(na.action=="na.roughfix"){
			warning("Replacing predictor categories not found in training data with most common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with predictor categories not found in the training data")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
			}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}		

	## make predictions
	pred<-rep(-9999,nrow(qdata.x))
	pred[not9]<-predict.gbm(	object=SGB,
						newdata=qdata.x.not9,
						n.trees=n.trees,
						type="response",
						single.tree=FALSE)
}
SGB.PRED<-data.frame(cbind(obs=qdata.y,pred=pred))

row.names(SGB.PRED)<-row.names(qdata)

return(SGB.PRED)
}



#############################################################################################
############### RF - Model Creation - Categorical (binary) response ######################
#############################################################################################

rF.binary<-function(	qdata,
				predList,
				response.name,
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
				seed=NULL){

## This function generates a presence/absence (binary categorical) model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model


if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]
qdata.y[qdata.y>0]<-1
qdata.y<-as.factor(qdata.y)

print("about to start tuning")

if(is.null(mtry)){

	A<-list(	x=qdata.x, y=qdata.y,
			doBest=FALSE,
			importance=TRUE,
			proximity=TRUE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

print("finished tuning")

A<-list(	x=qdata.x, y=qdata.y,
		importance=TRUE,
		proximity=TRUE,
		mtry=mtry,
		ntree=ntree,
		replace=replace,
		strata=strata,
		sampsize=sampsize)

A<-A[!sapply(A, is.null)]

RF<-do.call("randomForest", A)

#rast.factors<-names(is.fact[is.fact])

RF$response<-response.name

return(RF)
}



#############################################################################################
################## RF - Predict - Categorical (binary) response ##########################
#############################################################################################

prediction.rF.binary<-function(	prediction.type,
						qdata,
						train,
						response.name=RF$response,
						RF,
						na.action="na.omit"
						){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}

predList<-row.names(RF$importance)

qdata.x<-qdata[,match(predList,names(qdata))]

qdata.y<-qdata[,response.name]
qdata.y[qdata.y>0]<-1
qdata.y<-as.factor(qdata.y)

if(prediction.type=="OOB"){
	pred<-predict(RF, type="vote")[,"1"]
}
if(prediction.type=="TEST"){

	## check -9999 

	not9<-!apply(qdata.x[,]==-9999,1,any)
	qdata.x.not9<-qdata.x[not9,]
	qdata.x.not9<-data.frame(qdata.x.not9)
	qdata.x.not9<-qdata.x.not9[,match(rownames(RF$importance),names(qdata.x.not9))]

	## checking for NA's

	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("predictor",names(qdata.x)[p],"contains NA's"))
		}
		if(na.action=="na.roughfix"){
			warning("Replacing NA predictors with median value or most common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with NA predictors")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
			}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}

	## Factors

	if(!is.null(RF$levels)){
			missing.levels<-RF$levels
			for(p in names(RF$levels)){
				missing.levels[[p]]<-unique(qdata.x.not9[,p][!qdata.x.not9[,p]%in%RF$levels[[p]]])
				qdata.x.not9[,p]<-factor(qdata.x.not9[,p],levels=RF$levels[[p]])}
		}

	## checking for missing factor categories

	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("categorical factored predictor",names(qdata.x)[p],"contains  levels",paste(missing.levels[[names(qdata.x)[p]]],collapse=", "), "not found in training data"))
		}

		if(na.action=="na.roughfix"){
			warning("Replacing predictor categories not found in training data with most common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with predictor categories not found in the training data")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
				}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}

	pred<-rep(-9999,nrow(qdata.x))
	pred[not9]<-predict(RF, qdata.x.not9,type="vote")[,"1"]
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),
					pred=pred))

row.names(RF.PRED)<-row.names(qdata)

return(RF.PRED)
}



#############################################################################################
#################### RF - Model Creation - Continuous Response ##############################
#############################################################################################

rF.continuous<-function(	qdata,
					predList,
					response.name,
					ntree=500,
					mtry=NULL,
					replace=TRUE,
					strata=NULL,
					sampsize = NULL,
					seed=NULL){


## This function generates a continuous response model using Random Forests.
##	Inputs: Full dataset, training indices, predictor names, response name, and seed (optional)
##	Output: Random Forest model

if(!is.null(seed)){
	set.seed(seed)}

qdata.x<-qdata[,match(predList,names(qdata))]

is.fact<-sapply(qdata.x,is.factor)
if(any(is.fact)){
	qdata.x[,is.fact]<-lapply(qdata.x[,is.fact,drop=FALSE],factor)}

qdata.y<-qdata[,response.name]



if(is.null(mtry)){

	A<-list(	x=qdata.x, y=qdata.y,
			doBest=FALSE,
			importance=TRUE,
			proximity=TRUE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

A<-list(	x=qdata.x, y=qdata.y,
		importance=TRUE,
		proximity=TRUE,
		mtry=mtry,
		ntree=ntree,
		replace=replace,
		strata=strata,
		sampsize=sampsize)

A<-A[!sapply(A, is.null)]

RF<-do.call("randomForest", A)

#is.fact<-sapply(qdata.x,is.factor)
#rast.factors<-names(is.fact[is.fact])

RF$response<-response.name

return(RF)

}


#############################################################################################
########################## RF - Predict - Continuous Response ###############################
#############################################################################################

prediction.rF.continuous<-function(	prediction.type,
						qdata,
						response.name=RF$response,
						RF,
						na.action="na.omit"
						){

## This function makes predictions to test data for Random Forest continuous model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

print("MAKING PREDICTIONS SUBFUNCTION")

if(is.null(response.name)){
	stop("must provide response name")}

predList<-row.names(RF$importance)

#print("predList:")
#print(predList)

#print("names(qdata):")
#print(names(qdata))

qdata.x<-qdata[,match(predList,names(qdata))]
qdata.y<-qdata[,response.name]

if(prediction.type=="OOB"){
	pred<-RF$predicted
}
if(prediction.type=="TEST"){

	## checking -9999 

	not9<-!apply(qdata.x[,]==-9999,1,any)

	#print("not9:")
	#print(not9)

	qdata.x.not9<-qdata.x[not9,]
	qdata.x.not9<-data.frame(qdata.x.not9)
	qdata.x.not9<-qdata.x.not9[,match(rownames(RF$importance),names(qdata.x.not9))]

	## checking for NA's
	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("predictor",names(qdata.x)[p],"contains NA's"))
		}
		if(na.action=="na.roughfix"){
			warning("Replacing NA predictors with median value or most common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with NA predictors")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
			}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}

	## Factors

	print("factor stuff")

	if(!is.null(RF$levels)){
		missing.levels<-RF$levels

		for(p in names(RF$levels)){
			missing.levels[[p]]<-unique(qdata.x.not9[,p][!qdata.x.not9[,p]%in%RF$levels[[p]]])
			qdata.x.not9[,p]<-factor(qdata.x.not9[,p],levels=RF$levels[[p]])}

		}

	## Checking for missing factor categories

	if(any(is.na(qdata.x.not9))){
		qdata.x.na<-apply(is.na(qdata.x.not9),2,sum)
		for(p in (1:length(qdata.x.na))[qdata.x.na>0]){
			warning(paste("categorical factored predictor",names(qdata.x)[p],"contains  levels",paste(missing.levels[[names(qdata.x)[p]]],collapse=", "), "not found in training data"))
		}
		print("na.action:")
		print(na.action)
		if(na.action=="na.roughfix"){
			warning("Replacing predictor categories not found in training data withmost common category")
			qdata.x.not9<-na.roughfix(qdata.x.not9)
		}else{
			if(na.action=="na.omit"){
				warning("Returning -9999 for data points with predictor categories not found in the training data")
				notna<-!apply(is.na(qdata.x.not9),1,any)
				qdata.x.not9<-qdata.x.not9[notna,]
				not9[not9]<-notna
				}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
	}

	pred<-rep(-9999,nrow(qdata.x))
	pred[not9]<-predict(RF, qdata.x.not9)
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),
					pred=pred))

row.names(RF.PRED)<-row.names(qdata)


return(RF.PRED)
}



#############################################################################################
######################## Model Creation - Wrapper Function ##################################
#############################################################################################


create.model<-function(	qdata,
				model.type=NULL,		# "RF", "GAM", "SGB"
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				predList,
				response.name=NULL,
				response.type,			# "binary", "continuous",
				seed=NULL,
				keep.data=TRUE,

			# RF arguments:
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
			
			# SGB arguments:
				n.trees=NULL,                 # number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				train.fraction = 1.0,       	# fraction of data for training,
                  	n.minobsinnode = 10,         	# minimum total weight needed in each node
				var.monotone = NULL

){


### Set Seed ###
if(!is.null(seed)){
	set.seed(seed)}

if(model.type=="RF"){
	if(response.type=="binary"){
		print("calling rF.binary")
		model.obj<-rF.binary(	qdata=qdata,
						predList=predList,
						response.name=response.name,
						ntree=ntree,
						mtry=mtry,
						replace=replace,
						strata=strata,
						sampsize=sampsize,
						seed=NULL)}
	if(response.type=="continuous"){
		print("calling rF.continuous")
		model.obj<-rF.continuous(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							strata=strata,
							sampsize=sampsize,
							seed=NULL)}
}

if(model.type=="SGB"){

	print("calling model.SGB")
	model.obj<-model.SGB(	qdata=qdata,
					predList=predList,
					response.name=response.name,
					seed=NULL,
					response.type=response.type, 				
					n.trees=n.trees,                 	# number of trees
					shrinkage=shrinkage,   	      	# shrinkage or learning rate,
                  		interaction.depth=interaction.depth,# 1: additive model, 2: two-way interactions, etc.
					bag.fraction=bag.fraction,          # subsampling fraction, 0.5 is probably best
					train.fraction=train.fraction,      # fraction of data for training,
                  		n.minobsinnode=n.minobsinnode,      # minimum total weight needed in each node
					keep.data=keep.data,
					var.monotone = var.monotone)

}

return(model.obj)
}

#############################################################################################
######################## Model Prediction - Wrapper Function ##################################
#############################################################################################


prediction.model<-function(	model.obj,
					model.type,
					qdata,
					folder=NULL,		# No ending slash, to output to working dir = getwd()
					response.name,
					response.type,
					unique.rowname="row_index",	# Row identifier

				# Model Evaluation Arguments
					prediction.type=NULL,
					MODELpredfn=NULL,
					na.action="na.omit",	# also used for mapping
					v.fold=FALSE,

				# SGB arguments
					n.trees
){


### make predictions ###
if(prediction.type!="CV"){
	if(model.type=="RF"){
			
		if(response.type=="binary"){
			PRED<-prediction.rF.binary(	prediction.type=prediction.type,
								qdata=qdata,
								response.name=response.name,
								RF=model.obj,
								na.action=na.action)}
		if(response.type=="continuous"){
			PRED<-prediction.rF.continuous(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj,
									na.action=na.action)}
	}

	if(model.type=="SGB"){
		print("calling prediction.SGB")
		PRED<-prediction.SGB(	prediction.type=prediction.type,
						qdata=qdata,
						response.name=response.name,
						SGB=model.obj,
						na.action=na.action,
						n.trees=n.trees)
	}
}else{
	print(paste("Begining ",v.fold,"-fold cross validation:",sep=""))
	n.data=nrow(qdata)
	n.per.fold<-floor(n.data/v.fold)
	cv.index<-sample(rep(1:v.fold,(n.per.fold+1))[1:n.data])

	###debugging###

	#cv.info<-cbind(qdata,cv.index)
	#cv.info$Vfold<-cv.index
	#write.table(cv.info,file="Vfold.csv",sep=",",row.names=FALSE)

	###end debugging###
	
	PRED<-data.frame(matrix(0,0,2))
	names(PRED)<-c("obs","pred")

	if(model.type=="RF"){
		ntree<-model.obj$ntree
		mtry<-model.obj$mtry
		replace<-model.obj$call$replace
		predList<-row.names(model.obj$importance)}
	if(model.type=="SGB"){
		shrinkage<-model.obj$shrinkage
		interaction.depth<-model.obj$interaction.depth
		bag.fraction<-model.obj$bag.fraction
		train.fraction<-model.obj$train.fraction
		n.minobsinnode<-model.obj$n.minobsinnode
		predList<-model.obj$var.names}


	for(i in 1:v.fold){
		train.cv<-(1:nrow(qdata))[cv.index!=i]
		if(model.type=="RF"){
			if(response.type=="binary"){
				RF.cv<-rF.binary(	qdata=qdata[train.cv,],
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							seed=NULL)
				print(paste("        calling prediction.rF.binary for fold",i))
				PRED.cv<-prediction.rF.binary(	prediction.type="TEST",
										qdata=qdata[-train.cv,],
										response.name=response.name,
										RF=RF.cv,
										na.action=na.action)
			}
			if(response.type=="continuous"){
				#print("generating new model")
				RF.cv<-rF.continuous(	qdata=qdata[train.cv,],
								predList=predList,
								response.name=response.name,
								ntree=ntree,
								mtry=mtry,
								replace=replace,
								seed=NULL)
				print(paste("        calling prediction.rF.continuous for fold",i))
				PRED.cv<-prediction.rF.continuous(	prediction.type="TEST",
										qdata=qdata[-train.cv,],
										response.name=response.name,
										RF=RF.cv,
										na.action=na.action)
			}
		}
		if(model.type=="SGB"){
			#print(paste("      making model for fold",i))

			
			SGB.cv<-model.SGB(	qdata=qdata[train.cv,],
							predList=predList,
							response.name=response.name,
							seed=NULL,
							response.type=response.type, 				
							n.trees=n.trees,	               	# number of trees
							shrinkage=shrinkage,   	      	# shrinkage or learning rate,
                  				interaction.depth=interaction.depth,# 1: additive model, 2: two-way interactions, etc.
							bag.fraction=bag.fraction,          # subsampling fraction, 0.5 is probably best
							train.fraction=train.fraction,      # fraction of data for training,
                  				n.minobsinnode=n.minobsinnode       # minimum total weight needed in each node
							)
			#print(paste("      making SGB predictions for fold",i))
			PRED.cv<-prediction.SGB(	prediction.type="TEST",
								qdata=qdata[-train.cv,],
								response.name=response.name,
								SGB=SGB.cv,
								na.action=na.action,
								n.trees=n.trees)
		}
		PRED.cv$VFold<-i
		PRED<-rbind(PRED,PRED.cv)
		#print(paste("     ending fold",i))
	}
PRED<-PRED[match(row.names(qdata),row.names(PRED)),]
}

PRED<-cbind(rownames(PRED),PRED)
colnames(PRED)[1]<-unique.rowname

PREDICTIONfn<-paste(MODELpredfn,".csv",sep="")
write.table(PRED,file=PREDICTIONfn,sep=",",row.names=FALSE)
return(PRED)
}

#############################################################################################
#############################################################################################
###################### Production Prediction - Random Forests ###############################
#############################################################################################
#############################################################################################


production.prediction<-function(	model.obj,
						model.type,
						rastLUT,
						na.action="na.omit",
						response.type,
						numrows=500,	
						map.sd=FALSE,
						asciifn,
						asciifn.mean,
						asciifn.stdev,
						asciifn.coefv,
						n.trees){



#####################################################################################
########################## Extract predictor names ##################################
#####################################################################################


## Make sure raster predictor names match names in training/test data.

if(model.type=="RF"){
	predList<-row.names(model.obj$importance)}
if(model.type=="SGB"){
	predList<-model.obj$var.names}

predLUT<-rastLUT[match(predList, rastLUT[,2]),]
if(any(predList!=predLUT[,2])){
	stop("predictor names from model do not match short names in rastLUT")}	

## Gets the raster layers and/or stacks necessary to run model.
rastnm.all<-unique(predLUT[,1])

########################################################################################
#################################### Set up ############################################
########################################################################################

## Initialize variables
rowcnt <- 0		# The total count of rows each time through loop
final <- FALSE	# Loop testing variable
offset <- 0		# Offset variable for importing rows

## Open raster
sp.rast <- open.SpatialGDAL(rastnm.all[1])

## Count rasters
rastcnt <- rep(1,length(rastnm.all))

## Get name of raster
rastnm = predLUT[1,2]

## Set variables for header of ASCII file
rastnm1 <- rastnm
ncols <- sp.rast@grid@cells.dim[1]
nrows <- sp.rast@grid@cells.dim[2]
xllcorner <- sp.rast@bbox[1]
yllcorner <- sp.rast@bbox[2]
cellsize <- sp.rast@grid@cellsize[1]
NODATA_value <- -9999

## Get dimensions of raster
sp.rast.dim<-dim(sp.rast@grod)
if(length(sp.rast.dim)==3){rastcnt[1]<-sp.rast.dim[3]}

## Close raster
close(sp.rast)

## Writes out header to an ASCII file
write(paste("ncols", ncols, sep = "         "), file=asciifn)
write(paste("nrows", nrows, sep = "         "), file=asciifn, append=TRUE)
write(paste("xllcorner", xllcorner, sep = "     "), file=asciifn, append=TRUE)
write(paste("yllcorner", yllcorner, sep = "     "), file=asciifn, append=TRUE)
write(paste("cellsize", cellsize, sep = "      "), file=asciifn, append=TRUE)
write(paste("NODATA_value", NODATA_value, sep = "  "), file=asciifn, append=TRUE)


if(map.sd && model.type=="RF" && response.type=="continuous"){
	## Writes out header to an ASCII file
	write(paste("ncols", ncols, sep = "         "), file=asciifn.mean)
	write(paste("nrows", nrows, sep = "         "), file=asciifn.mean, append=TRUE)
	write(paste("xllcorner", xllcorner, sep = "     "), file=asciifn.mean, append=TRUE)
	write(paste("yllcorner", yllcorner, sep = "     "), file=asciifn.mean, append=TRUE)
	write(paste("cellsize", cellsize, sep = "      "), file=asciifn.mean, append=TRUE)
	write(paste("NODATA_value", NODATA_value, sep = "  "), file=asciifn.mean, append=TRUE)

	## Writes out header to an ASCII file
	write(paste("ncols", ncols, sep = "         "), file=asciifn.stdev)
	write(paste("nrows", nrows, sep = "         "), file=asciifn.stdev, append=TRUE)
	write(paste("xllcorner", xllcorner, sep = "     "), file=asciifn.stdev, append=TRUE)
	write(paste("yllcorner", yllcorner, sep = "     "), file=asciifn.stdev, append=TRUE)
	write(paste("cellsize", cellsize, sep = "      "), file=asciifn.stdev, append=TRUE)
	write(paste("NODATA_value", NODATA_value, sep = "  "), file=asciifn.stdev, append=TRUE)

	## Writes out header to an ASCII file
	write(paste("ncols", ncols, sep = "         "), file=asciifn.coefv)
	write(paste("nrows", nrows, sep = "         "), file=asciifn.coefv, append=TRUE)
	write(paste("xllcorner", xllcorner, sep = "     "), file=asciifn.coefv, append=TRUE)
	write(paste("yllcorner", yllcorner, sep = "     "), file=asciifn.coefv, append=TRUE)
	write(paste("cellsize", cellsize, sep = "      "), file=asciifn.coefv, append=TRUE)
	write(paste("NODATA_value", NODATA_value, sep = "  "), file=asciifn.coefv, append=TRUE)
}


## Check all rasts for consistency

#print("checking rasts")

if(length(rastnm.all)>1){
for (r in 2:length(rastnm.all)) {

	rast<-rastnm.all[r]
	## Open raster
	sp.rast <- open.SpatialGDAL(rast)
	sp.rast.dim<-dim(sp.rast@grod)
	if(length(sp.rast.dim)==3){rastcnt[r]<-sp.rast.dim[3]}


	## Get name of raster
	rastnm <- basename(rast)
	rastnm <- strsplit(rastnm,".img")

# Check if all rasters have the same cellsize and extent.
	if (ncols != sp.rast@grid@cells.dim[1]){
		stop("Number of columns of", rastnm,"= ",sp.rast@grid@cells.dim[1]," Number of columns of",rastnm1,"= ",ncols)}
	if (nrows != sp.rast@grid@cells.dim[2]){
		stop("Number of rows of", rastnm,"= ",sp.rast@grid@cells.dim[2]," Number of rows of",rastnm1,"= ",nrows)}
	if (xllcorner != sp.rast@bbox[1]){
		warning("The xllcorner of", rastnm, "= ",sp.rast@bbox[1], " The xllcorner of",rastnm1,"= ",xllcorner,immediate.=TRUE)}
	if (abs(xllcorner - sp.rast@bbox[1]) > cellsize){
		warning("These images are misregistered by more than one cell",immediate.=TRUE)}
	if (yllcorner != sp.rast@bbox[2]){
		warning("The yllcorner of", rastnm, "= ",sp.rast@bbox[2], " The yllcorner of",rastnm1,"= ",yllcorner,immediate.=TRUE)}
	if (abs(yllcorner - sp.rast@bbox[2]) > cellsize){
		warning("These images are mis-registered by more than one cell",immediate.=TRUE)}
	if (cellsize != sp.rast@grid@cellsize[1]){
		warning("The cellsize of", rastnm,"= ",sp.rast@grid@cellsize[1]," The cellsize of",rastnm1,"= ",cellsize,immediate.=TRUE)}
		
	close(sp.rast)
}}

#print("done checking rasts")


#################################################################################

while (!final){

	print(paste("numrows =",numrows))
	print(paste("rowcnt =",rowcnt))

	if (rowcnt+numrows >= nrows) {
		numrows <- numrows - ((rowcnt+numrows) - nrows)
		offset <- rowcnt
		final <- TRUE
	}

	rast <- rastnm.all[1]
	rastnm <- predLUT[predLUT[,1]==rast,2][1]

	predcnt <- table(predLUT[,1])[match(rastnm.all,names(table(predLUT[,1])))]
	preds <- data.frame(matrix(-9999,numrows*ncols,sum(predcnt)))

	##define dematrix() function to turn matrix into vector
	dematrix <- function(m){m[1:length(m)]}

	p=0

	for(r in 1:length(rastnm.all)){

		#print(paste("rastnm =",rastnm.all[r]))

		predLUT.r<-predLUT[predLUT[,1]==rastnm.all[r],]

		openrast <- GDAL.open(rastnm.all[r])

		bands <- predLUT.r[,3]
		pred <- getRasterData(openrast, offset = c(offset, 0), band=bands, region.dim=c(numrows, ncols),as.is=TRUE)

		GDAL.close(openrast)
 
		if(predcnt[r]==1){
			preds[,p+1]<-pred[1:length(pred)]
			names(preds)[p+1]<-predLUT.r[,2]
		}else{
			if(numrows==1){
				preds[,(p+1):(p+predcnt[r])]<-pred
			}else{
				preds[,(p+1):(p+predcnt[r])]<-apply(pred,3,dematrix)}
			names(preds)[(p+1):(p+predcnt[r])]<-predLUT.r[,2]
			}

		p<-p+predcnt[r]

	}


	#preds<-preds[,match(rastLUT[,2],prednm.uncut)]
	preds<-preds[,match(predLUT[,2],names(preds))]


	
	test_pred<-rep(-9999,nrow(preds))

	## -9999 

	if(any(is.na(preds))){

		if(na.action=="na.roughfix"){
			warning("NA values in predictors repaced with median or most common category")
			preds<-na.roughfix(preds)
		}else{
			if(na.action=="na.omit"){
				warning("NA values in predictors repaced with  NODATA value of -9999 and no predictions will be made on these pixels")
				preds[is.na(preds)] <- -9999
				}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
		}
		#warning("NA values in predictors are being repaced with designated No Data value of -9999") 
		#preds[is.na(preds)]<- -9999
	}

	not9<-!apply(preds==-9999,1,any)
	all9<-!any(not9==TRUE)

	if(!all9){
		preds.not9<-preds[not9,]
		preds.not9<-data.frame(preds.not9)
		preds.not9<-preds.not9[,match(predList,names(preds.not9))]

		## Factors

		if(!is.null(model.obj$levels)){
			missing.levels<-model.obj$levels
			for(p in names(model.obj$levels)){
				missing.levels[[p]]<-unique(preds.not9[,p][!preds.not9[,p]%in%model.obj$levels[[p]]])
				preds.not9[,p]<-factor(preds.not9[,p],levels=model.obj$levels[[p]])}
		}

		if(any(is.na(preds.not9))){
			preds.na<-apply(is.na(preds.not9),2,sum)
			for(p in (1:length(preds.na))[preds.na>0]){
				warning(paste("categorical factored predictor",names(preds)[p],"contains levels",
					paste(missing.levels[[names(preds)[p]]],collapse=", "),  "not found in training data"))
		}

			if(na.action=="na.roughfix"){
				warning("Replacing categorical factored predictor levels not found in training data, with most common category that is found in training")
				preds.not9<-na.roughfix(preds.not9)
			}else{
				if(na.action=="na.omit"){
					warning("Returning -9999 for data points with predictor categories not found in the training data")
					notna<-!apply(is.na(preds.not9),1,any)
					preds.not9<-preds.not9[notna,]
					not9[not9]<-notna
					}else{stop("na.action must be either 'na.roughfix' or 'na.omit'")}
			}
		}
	

		## Model predictions

		print("making predictions")

		if(model.type=="RF"){
			if(response.type=="binary"){
				test_pred[not9]<-signif(predict(model.obj, preds.not9,type="vote")[,"1"],2)}
			if(response.type=="continuous"){
				test_pred[not9]<-predict(model.obj, preds.not9)}
		}
		if(model.type=="SGB"){
			test_pred[not9]<-predict.gbm(	object=model.obj,
								newdata=preds.not9,
								n.trees=n.trees,
								type="response",
								single.tree=FALSE)
		}
	}

	write(test_pred, file = asciifn, ncol=ncols, append=TRUE, sep=" ")
	
	#print("starting standard deviation stuff")

	if(map.sd && model.type=="RF" && response.type=="continuous"){

		print(paste("Are all pixels -9999?",all9))

		test_stdev<-rep(-9999,nrow(preds))
		test_mean <-rep(-9999,nrow(preds))		
		test_coefv<-rep(-9999,nrow(preds))

		if(!all9){

			print(paste("number of not -9999 pixels =",nrow(preds.not9)))			

			test_mat<-predict(model.obj, preds.not9, predict.all=TRUE)$individual
			#print("I've made the predictions! Now for the stdev...")
			test_stdev[not9]<-apply(test_mat,1,sd)
			test_mean[not9]<-apply(test_mat,1,mean)
			test_coefv[not9]<-test_stdev[not9]/test_mean[not9]
			test_coefv[test_stdev==0]<-0
			test_coefv[test_mean==0]<-0
			rm(test_mat)
		}

		write(test_mean,  file = asciifn.mean, ncol=ncols, append=TRUE, sep=" ")
		write(test_stdev, file = asciifn.stdev, ncol=ncols, append=TRUE, sep=" ")
		write(test_coefv, file = asciifn.coefv, ncol=ncols, append=TRUE, sep=" ")}

	rowcnt<-rowcnt + numrows
	offset<-rowcnt	
	}
}

