
#############################################################################################
#############################################################################################
################################## Importance Plot ##########################################
#############################################################################################
#############################################################################################

########## check model type #################

check.model.type<-function(model.obj){
	model.type.long<-attr(model.obj,"class")
	if(is.null(model.type.long)){model.type.long <- "unknown"}

	if("randomForest"%in%model.type.long){
		model.type<-"RF"
	}else{
		if("gbm"%in%model.type.long){
			model.type<-"SGB"
		}else{
			model.type<-"unknown"
		}
	}

	if(model.type=="unknown"){
		stop("model.obj is of unknown type")}

	return(model.type)
}


########## extract importance SGB ############

imp.extract.sgb<-function(model.obj,imp.type=NULL){
	imp.type.gbm<-switch(imp.type,"1"=permutation.test.gbm,"2"=relative.influence)
	IMP.SGB<-summary(model.obj, method=imp.type.gbm, plotit=FALSE)
	names(IMP.SGB)<-c("pred","imp")
	row.names(IMP.SGB)<-IMP.SGB$pred
	IMP.SGB<-IMP.SGB[order(IMP.SGB$imp,decreasing=FALSE),]
	return(IMP.SGB)
}

######### extract importance RF ############

imp.extract.rf<-function(model.obj,imp.type=NULL,class=NULL){
	IMP.RF<-importance(model.obj,type=imp.type,class=class)
	IMP.RF<-data.frame(pred=rownames(IMP.RF),imp=IMP.RF[,1])
	IMP.RF<-IMP.RF[order(IMP.RF[,2],decreasing=FALSE),]
	return(IMP.RF)
}

######### scale importance ###############

imp.scale<-function(IMP,scale.by="max"){

	if(scale.by=="max"){
		IMP$imp[IMP$imp<0]<-0
		IMP$imp<-(IMP$imp/max(IMP$imp))*1
	}

	if(scale.by=="sum"){
		IMP$imp[IMP$imp<0]<-0
		IMP$imp<-(IMP$imp/sum(IMP$imp))*1
	}

	return(IMP)
}



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
#############################################################################################
###################################### Diagnostics ##########################################
#############################################################################################
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

if(response.type == "categorical"){

	library("PresenceAbsence")
	CMXfn<-paste(MODELpredfn,"_cmx.csv",sep="")

	LEVELS<-unique(c(levels(PRED$pred),levels(PRED$obs)))

	CMX<-table(	predicted = factor(PRED$pred,levels=LEVELS),
				observed = factor(PRED$obs,levels=LEVELS))


	CMX.out<-matrix("cmx",nrow(CMX)+4,ncol(CMX)+4)
	CMX.out[1,(1:ncol(CMX))+2]<-"observed"
	CMX.out[2,(1:ncol(CMX))+2]<-colnames(CMX)

	CMX.out[(1:nrow(CMX))+2,1]<-"predicted"
	CMX.out[(1:nrow(CMX))+2,2]<-rownames(CMX)

	CMX.out[(1:nrow(CMX))+2,(1:ncol(CMX))+2]<-CMX

	###Kappa###
	KAPPA<-signif(Kappa(CMX),6)
	CMX.out[1,1:2]<-names(KAPPA)
	CMX.out[2,1]<-KAPPA[1,1]
	CMX.out[2,2]<-KAPPA[1,2]

	###Totals###
	CMX.out[1:2,ncol(CMX.out)-1]<-"total"
	CMX.out[nrow(CMX.out)-1,1:2]<-"total"
	
	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)-1]<-apply(CMX,1,sum)
	CMX.out[nrow(CMX.out)-1,(1:ncol(CMX))+2]<-apply(CMX,2,sum)
	CMX.out[nrow(CMX.out)-1,ncol(CMX.out)-1]<-sum(CMX)

	###marginals###
	CMX.diag<-diag(CMX)

	CMX.out[1:2,ncol(CMX.out)]<-"User"
	CMX.out[nrow(CMX.out),1:2]<-"Producer"

	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)]<-CMX.diag/apply(CMX,1,sum)
	CMX.out[nrow(CMX.out),(1:ncol(CMX))+2]<-CMX.diag/apply(CMX,2,sum)

	###pcc###
	CMX.out[ncol(CMX.out)-1,ncol(CMX.out)]<-"PCC"
	CMX.out[ncol(CMX.out),ncol(CMX.out)-1]<-"PCC"

	CMX.out[ncol(CMX.out),ncol(CMX.out)]<-sum(CMX.diag)/sum(CMX)


	write.table(CMX.out,file=CMXfn,sep=",",row.names=FALSE,col.names = FALSE)
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

#
if(!"none"%in%device.type){
for(i in 1:length(device.type)){

	### Output filenames ###

	SCATTERPLOTfn<-paste(MODELpredfn,"_scatterplot",sep="")
	IMPORTANCEfn<-paste(MODELpredfn,"_importance",sep="")
	ERRORfn<-paste(MODELpredfn,"_error",sep="")
	THRESHOLDPLOTSfn<-paste(MODELpredfn,"_thresholdplots",sep="")

	###################   model.obj Based   ##############################


	if(model.type=="RF"){

		if(i==1){library(randomForest)}

		###Importance Plot###
		#print(paste("IMPORTANCEfn =",IMPORTANCEfn)

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(IMPORTANCEfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(IMPORTANCEfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(IMPORTANCEfn,".pdf",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=paste(IMPORTANCEfn,".wmf",sep=""),width = device.width, height = device.height, 
								pointsize = 12,restoreConsole = TRUE)}
 
		opar<-par(cex=cex)
		varImpPlot(model.obj,main="Relative Influence",cex=cex)
		mtext(main,side=3,line=-4,cex=1.3*cex,outer=TRUE)
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}

		###Importance Plot - category graphs###
		if(response.type%in%c("binary","categorical")){
			Nplots<-table(model.obj$y)

			CATnames<-names(Nplots)
			CATfn<-make.names(CATnames)
			Nfig<-length(CATnames)

			for(j in 1:Nfig){

				IMPFIGfn<-paste(IMPORTANCEfn,CATfn[j],sep="_")

				if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
				if(device.type[i]=="jpeg"){jpeg(filename=paste(IMPFIGfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
				if(device.type[i]=="postscript"){postscript(file=paste(IMPFIGfn,".ps",sep=""),width = device.width, height = device.height)}
				if(device.type[i]=="pdf"){pdf(file=paste(IMPFIGfn,".pdf",sep=""),width = device.width, height = device.height)}
				if(device.type[i]=="win.metafile"){win.metafile(filename=paste(IMPFIGfn,".wmf",sep=""),width = device.width, height = device.height, 
										pointsize = 12,restoreConsole = TRUE)}
				opar<-par(cex=cex)			

				varImpPlot(model.obj,cex=cex,type=1,class=CATnames[j],main="")
				mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
				mtext(paste("Relative Influence -",CATnames[j],"-",Nplots[j],"plots"),side=3,line=-3.5,cex=1.3*cex,outer=TRUE)

				par(opar)

				if(!device.type[i]%in%c("default","none")){dev.off()}
			}
		}


		###Error Plot###
		if(response.type == "continuous"){
			if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
			if(device.type[i]=="jpeg"){jpeg(filename=paste(ERRORfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
			if(device.type[i]=="postscript"){postscript(file=paste(ERRORfn,".ps",sep=""),width = device.width, height = device.height)}
			if(device.type[i]=="pdf"){pdf(file=paste(ERRORfn,".pdf",sep=""),width = device.width, height = device.height)}
			if(device.type[i]=="win.metafile"){win.metafile(filename=paste(ERRORfn,".wmf",sep=""),width = device.width, height = device.height, 
									pointsize = 12,restoreConsole = TRUE)}
			opar<-par(cex=cex)
		
			Nmax<-length(model.obj$mse)
			Nplots<-length(model.obj$predicted)

			plot(	1:Nmax,
				model.obj$mse,
				type="l", 
				xlab="ntree",ylab="MSE")
			mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
			mtext(paste("OOB -",Nplots,"plots"),side=3,line=-3.5,cex=1.1*cex,outer=TRUE)

			par(opar)
			if(!device.type[i]%in%c("default","none")){dev.off()}
		}

		###Error Plot - category graphs###
		if(response.type%in%c("binary","categorical")){
			Nplots<-table(model.obj$y)
			Nplots<-c(sum(Nplots),Nplots)

			CATnames<-colnames(model.obj$err.rate)
			CATfn<-make.names(CATnames)
			Nfig<-length(CATnames)
			Nmax<-nrow(model.obj$err.rate)

			for(j in 1:Nfig){

				ERRORFIGfn<-paste(ERRORfn,CATfn[j],sep="_")

				if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
				if(device.type[i]=="jpeg"){jpeg(filename=paste(ERRORFIGfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
				if(device.type[i]=="postscript"){postscript(file=paste(ERRORFIGfn,".ps",sep=""),width = device.width, height = device.height)}
				if(device.type[i]=="pdf"){pdf(file=paste(ERRORFIGfn,".pdf",sep=""),width = device.width, height = device.height)}
				if(device.type[i]=="win.metafile"){win.metafile(filename=paste(ERRORFIGfn,".wmf",sep=""),width = device.width, height = device.height, 
										pointsize = 12,restoreConsole = TRUE)}
				opar<-par(cex=cex)

				plot(	1:Nmax,
					model.obj$err.rate[,j],
					type="l", 
					xlab="ntree",ylab="OOB err.rate")
				mtext(main,side=3,line=-2,cex=1.3*cex,outer=TRUE)
				mtext(paste("OOB -",CATnames[j],"-",Nplots[j],"plots"),side=3,line=-3.5,cex=1.1*cex,outer=TRUE)

				par(opar)
				if(!device.type[i]%in%c("default","none")){dev.off()}
			}
		}
	}

	if(model.type=="SGB"){

		if(i==1){library(gbm)}
	
		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(IMPORTANCEfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(IMPORTANCEfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(IMPORTANCEfn,".pdf",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=paste(IMPORTANCEfn,".wmf",sep=""),width = device.width, height = device.height, 
								pointsize = 12,restoreConsole = TRUE)}

		opar<-par(las=1,mar=(c(5, 11, 4, 2) + 0.1),cex=cex)
		summary(model.obj)
		par(las=0)
		mtext("Relative Influence",side=3,line=.7,cex=1.5*cex)
		mtext(main,side=3,line=2.7,cex=1.5*cex)
		mtext("Predictors",side=2,line=10,cex=1*cex)
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}


	###################   PRED Based   ##############################

	### binary ###

	if(response.type == "binary"){
	
		if(i==1){library("PresenceAbsence")}

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(THRESHOLDPLOTSfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(THRESHOLDPLOTSfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(THRESHOLDPLOTSfn,".pdf",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=paste(THRESHOLDPLOTSfn,".wmf",sep=""),width = device.width, height = device.height, 
							pointsize = 12,restoreConsole = TRUE)}

		opar<-par(cex=cex)
		presence.absence.summary(PRED,main=main,legend.cex=cex,opt.legend.cex=cex)
		par(opar)

		if(!device.type[i]%in%c("default","none")){dev.off()}	
	}

### Continuous ###

	if(response.type == "continuous"){

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(SCATTERPLOTfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(SCATTERPLOTfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(SCATTERPLOTfn,".pdf",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="win.metafile"){win.metafile(filename=paste(SCATTERPLOTfn,".wmf",sep=""),width = device.width, height = device.height, 
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

		if(!device.type[i]%in%c("default","none")){dev.off()}
	}
}
#
}
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
					response.name=deparse(substitute(SGB$response.name)),
					SGB,
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
		for(p in names(SGB$levels)){
			qdata.x[,p]<-factor(qdata.x[,p],levels=SGB$levels[[p]])
		}
	}

	pred<-predict.gbm(	object=SGB,
					newdata=qdata.x,
					n.trees=n.trees,
					type="response",
					single.tree=FALSE)
}

if(prediction.type=="TEST"){
	
	pred<-predict.gbm(	object=SGB,
					newdata=qdata.x,
					n.trees=n.trees,
					type="response",
					single.tree=FALSE)
}
SGB.PRED<-data.frame(cbind(obs=qdata.y,pred=pred))

return(SGB.PRED)
}



#############################################################################################
########################## RF - Model Creation - Binary response ############################
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
if(!is.numeric(qdata.y)){
	stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
qdata.y[qdata.y>0]<-1
qdata.y[qdata.y<0]<-0
qdata.y<-as.factor(qdata.y)

#print("about to start tuning")

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

#print("finished tuning")

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
############################ RF - Predict - Binary response #################################
#############################################################################################

prediction.rF.binary<-function(	prediction.type,
						qdata,
						response.name=RF$response,
						RF
						){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}

if(prediction.type=="OOB"){
	pred<-predict(RF, type="vote")[,"1"]
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){
	predList<-row.names(RF$importance)

	qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.numeric(qdata.y)){
		stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
	qdata.y[qdata.y>0]<-1
	qdata.y[qdata.y<0]<-0
	qdata.y<-as.factor(qdata.y)

	pred<-predict(RF, qdata.x,type="vote")[,"1"]
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),
					pred=pred))

return(RF.PRED)
}

#############################################################################################
######################## RF - Model Creation - Categorical response #########################
#############################################################################################

rF.categorical<-function(	qdata,
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
if(!is.factor(qdata.y)){
		qdata.y<-as.factor(qdata.y)}

#print("about to start tuning")

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

#print("finished tuning")

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
####################### RF - Predict - Categorical response #################################
#############################################################################################

prediction.rF.categorical<-function(	prediction.type,
							qdata,
							response.name=RF$response,
							RF
							){

## This function makes predictions to test data for Random Forest presence/absence (binary 
## categorical) model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.

if(is.null(response.name)){
	stop("must provide response name")}


if(prediction.type=="OOB"){
	pred<-predict(RF)
	vote<-predict(RF, type="vote")
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){

	predList<-row.names(RF$importance)

	qdata.x<-qdata[,match(predList,names(qdata))]

	qdata.y<-qdata[,response.name]
	if(!is.factor(qdata.y)){
			qdata.y<-as.factor(qdata.y)}

	pred<-predict(RF, qdata.x)
	vote<-predict(RF, qdata.x,type="vote")

}

RF.PRED<-data.frame(	obs=qdata.y,
				pred=pred,
				vote)
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
						RF
						){

## This function makes predictions to test data for Random Forest continuous model.
##	Inputs: Training data, training data indices, predictor names, response variable name,
##			the Random Forest model and what to do if NAs are in predictors (default).
##	Output: Observed and predicted values.


if(is.null(response.name)){
	stop("must provide response name")}

if(prediction.type=="OOB"){
	pred<-RF$predicted
	qdata.y<-RF$y
}
if(prediction.type=="TEST"){
	predList<-row.names(RF$importance)
	qdata.x<-qdata[,match(predList,names(qdata))]
	pred<-predict(RF, qdata.x)
	qdata.y<-qdata[,response.name]
}

RF.PRED<-data.frame(	cbind(obs=as.numeric(as.character(qdata.y)),pred=pred))

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
		#print("calling rF.binary")
		model.obj<-rF.binary(	qdata=qdata,
						predList=predList,
						response.name=response.name,
						ntree=ntree,
						mtry=mtry,
						replace=replace,
						strata=strata,
						sampsize=sampsize,
						seed=NULL)}
	if(response.type=="categorical"){
		#print("calling rF.categorical")
		model.obj<-rF.categorical(	qdata=qdata,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							strata=strata,
							sampsize=sampsize,
							seed=NULL)}
	if(response.type=="continuous"){
		#print("calling rF.continuous")
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

	#print("calling model.SGB")
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
					v.fold=FALSE,

					na.action=na.action,			
					NA.ACTION=NA.ACTION,
					model.na.action=model.na.action,			
					model.NA.ACTION=model.NA.ACTION,

				# SGB arguments
					n.trees
){

#print("prediction wrapper function")
### make predictions ###
if(prediction.type!="CV"){
	if(model.type=="RF"){
			
		if(response.type=="binary"){
			#print("starting binary predictions")
			PRED<-prediction.rF.binary(	prediction.type=prediction.type,
								qdata=qdata,
								response.name=response.name,
								RF=model.obj)}
		if(response.type=="categorical"){
			#print("starting cat categorical predictions")
			PRED<-prediction.rF.categorical(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj)}
		if(response.type=="continuous"){
			#print("starting continuous predictions")
			PRED<-prediction.rF.continuous(	prediction.type=prediction.type,
									qdata=qdata,
									response.name=response.name,
									RF=model.obj)}
	}

	if(model.type=="SGB"){
		#print("calling prediction.SGB")
		PRED<-prediction.SGB(	prediction.type=prediction.type,
						qdata=qdata,
						response.name=response.name,
						SGB=model.obj,
						n.trees=n.trees)
	}

	#print("got predictions, checking rownames")

	#If prediction type is OOB, model was roughfix and diagnostics is omit, must remove na rows from predictions
	if(!is.null(NA.ACTION) && prediction.type=="OOB" && model.NA.ACTION=="roughfix" && NA.ACTION=="omit"){

		print(paste("     nrow(qdata) =",length(rownames(qdata))))
		print(paste("          remove these rows =",paste(model.obj$na.action,collapse=",")))

		print(paste("     nrow(PRED) =", nrow(PRED)))
		PRED<-PRED[-model.obj$na.action,] 
		print(paste("     nrow(PRED) =", nrow(PRED)))
	}

	rownames(PRED)<-rownames(qdata)


}else{
	print(paste("Begining ",v.fold,"-fold cross validation:",sep=""))
	n.data=nrow(qdata)
	n.per.fold<-floor(n.data/v.fold)
	cv.index<-sample(rep(1:v.fold,(n.per.fold+1))[1:n.data])
	
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

	###start folds###
	for(i in 1:v.fold){
		train.cv<-(1:nrow(qdata))[cv.index!=i]
		qdata.train.cv<-qdata[train.cv,]
		qdata.test.cv<-qdata[-train.cv,]

		###check for factor levels not found in other folds###
		if(!is.null(model.obj$levels)){
			predFactor<-names(model.obj$levels)
			invalid.levels<-vector("list", length=length(predFactor))
			names(invalid.levels)<-predFactor

			for(p in predFactor){
				#valid.levels<-c(model.obj$levels[[p]],NA)
				train.levels<-unique(as.character(qdata.train.cv[,p]))
				test.levels<- unique(as.character(qdata.test.cv[,p]))

				invalid<-test.levels[!(test.levels%in%c(train.levels,NA))] 

				#print(test.levels%in%c(train.levels,NA))
				#print(!test.levels%in%c(train.levels,NA))
				#print(test.levels[!(test.levels%in%c(train.levels,NA))])
				
				invalid.levels[[p]]<-invalid

				qdata.test.cv[,p]<-factor(qdata.test.cv[,p],levels=train.levels)


				print(paste("     p                    =",p))
				print(paste("     train levels[p]      =",paste(train.levels,collapse=",")))
				print(paste("     test levels[p]       =",paste(test.levels,collapse=",")))
				print(paste("     invalid              =",paste(invalid,collapse=",")))
				print(paste("     invalid levels[[p]]  =",paste(invalid.levels[[p]],collapse=",")))

				print(paste("          invalid length =",length(invalid.levels[[p]])))
				
				print(paste("     levels qdata.test.cv =",paste(levels(qdata.test.cv[,p]),collapse=",")))

				if(length(invalid.levels[[p]])>0){
					print("     NA action triggered")
					N.invalid<-sum(is.na(qdata.test.cv[,p]))
					warning(paste(	"Factored predictor",p,"in fold",i,
								"contains", N.invalid, "data point(s) of levels",paste(invalid.levels[[p]],collapse=", "),
								"not found in other folds, these levels treated as NA"))
				}
			}
			
			print(qdata.test.cv[,predList])
			NA.pred<-apply(qdata.test.cv[,predList],1,function(x){any(is.na(x))})
			print(paste(     "NA.pred =",any(NA.pred)))
			print("TEST TEST TEST")
			if(any(NA.pred)){
				print("about to treat na's")
				print("NA.ACTION:")
				print(NA.ACTION)
				print("na.action:")
				print(na.action)
				qdata.test.cv<-na.action(qdata.test.cv)
				print("finished treating NA's")}
			print(qdata.test.cv[,predList])
		}

		###build models and make predictions###
		if(model.type=="RF"){
			if(response.type=="binary"){
				RF.cv<-rF.binary(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							seed=NULL)
				
				print(paste("        calling prediction.rF.binary for fold",i))
				PRED.cv<-prediction.rF.binary(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}
			if(response.type=="categorical"){
				RF.cv<-rF.categorical(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							ntree=ntree,
							mtry=mtry,
							replace=replace,
							seed=NULL)
				print(paste("        calling prediction.rF.categorical for fold",i))
				PRED.cv<-prediction.rF.categorical(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}

			if(response.type=="continuous"){
				##("generating new model")
				RF.cv<-rF.continuous(	qdata=qdata.train.cv,
								predList=predList,
								response.name=response.name,
								ntree=ntree,
								mtry=mtry,
								replace=replace,
								seed=NULL)
				print(paste("        calling prediction.rF.continuous for fold",i))
				PRED.cv<-prediction.rF.continuous(	prediction.type="TEST",
										qdata=qdata.test.cv,
										response.name=response.name,
										RF=RF.cv)
			}
		}
		if(model.type=="SGB"){
			#print(paste("      making model for fold",i))

			
			SGB.cv<-model.SGB(	qdata=qdata.train.cv,
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
								qdata=qdata.test.cv,
								response.name=response.name,
								SGB=SGB.cv,
								n.trees=n.trees)
		}
		PRED.cv$VFold<-i
		PRED<-rbind(PRED,PRED.cv)
		#print(paste("     ending fold",i))
	}
PRED<-PRED[match(row.names(qdata),row.names(PRED),nomatch=0),] #cv only
}

PRED<-cbind(rownames(PRED),PRED)
colnames(PRED)[1]<-unique.rowname

PREDICTIONfn<-paste(MODELpredfn,".csv",sep="")
write.table(PRED,file=PREDICTIONfn,sep=",",row.names=FALSE)
return(PRED)
}

#############################################################################################
#############################################################################################
############################# Production Prediction #########################################
#############################################################################################
#############################################################################################


production.prediction<-function(	model.obj,
						model.type,
						rastLUT,
						na.action=NULL,
						NA.ACTION=NULL,
						response.type,
						numrows=500,	
						map.sd=FALSE,
						asciifn,
						asciifn.mean,
						asciifn.stdev,
						asciifn.coefv,
						make.img=TRUE,
						n.trees){



#####################################################################################
########################## Extract predictor names ##################################
#####################################################################################


## Make sure raster predictor names match names in training/test data.

if(model.type=="RF"){
	predList<-row.names(model.obj$importance)}
if(model.type=="SGB"){
	predList<-model.obj$var.names}

predLUT<-rastLUT[match(predList, rastLUT[,2]),] #only the predictors form the model, in same order as in model
if(any(predList!=predLUT[,2])){
	stop("predictor names from model do not match short names in rastLUT")}	

## Gets the filenames of raster layers and/or stacks necessary to run model.
rastnm.all<-unique(predLUT[,1])

########################################################################################
#################################### Set up ############################################
########################################################################################

## Initialize variables
rowcnt <- 0		# The total count of rows each time through loop
final <- FALSE	# Loop testing variable
offset <- 0		# Offset variable for importing rows

## Open first raster
sp.rast <- open.SpatialGDAL(rastnm.all[1])

## Count the number of raster files
#rastcnt <- rep(1,length(rastnm.all))

## Get the basename of first raster
rast.basenm1 <- basename(predLUT[1,1])
rast.basenm1 <- strsplit(rast.basenm1,".img")

## Set variables for header of ASCII file
ncols <- sp.rast@grid@cells.dim[1]
nrows <- sp.rast@grid@cells.dim[2]
xllcorner <- sp.rast@bbox[1]
yllcorner <- sp.rast@bbox[2]
cellsize <- sp.rast@grid@cellsize[1]
NODATA_value <- -9999

## Get dimensions of raster ### we never ended up using this for anything
#sp.rast.dim<-dim(sp.rast@grod)
#if(length(sp.rast.dim)==3){rastcnt[1]<-sp.rast.dim[3]}

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

	## Open raster
	sp.rast <- open.SpatialGDAL(rastnm.all[r])
	#sp.rast.dim<-dim(sp.rast@grod)
	#if(length(sp.rast.dim)==3){rastcnt[r]<-sp.rast.dim[3]}


	## Get name of raster
	rast.basenm <- basename(rastnm.all[r])
	rast.basenm <- strsplit(rast.basenm,".img")

# Check if all rasters have the same cellsize and extent.
	if (ncols != sp.rast@grid@cells.dim[1]){
		stop("Number of columns of", rast.basenm,"= ",sp.rast@grid@cells.dim[1]," Number of columns of",rast.basenm1,"= ",ncols)}
	if (nrows != sp.rast@grid@cells.dim[2]){
		stop("Number of rows of", rast.basenm,"= ",sp.rast@grid@cells.dim[2]," Number of rows of",rast.basenm1,"= ",nrows)}
	if (xllcorner != sp.rast@bbox[1]){
		warning("The xllcorner of", rast.basenm, "= ",sp.rast@bbox[1], " The xllcorner of",rast.basenm1,"= ",xllcorner,immediate.=TRUE)}
	if (abs(xllcorner - sp.rast@bbox[1]) > cellsize){
		warning("These images are misregistered by more than one cell",immediate.=TRUE)}
	if (yllcorner != sp.rast@bbox[2]){
		warning("The yllcorner of", rast.basenm, "= ",sp.rast@bbox[2], " The yllcorner of",rast.basenm1,"= ",yllcorner,immediate.=TRUE)}
	if (abs(yllcorner - sp.rast@bbox[2]) > cellsize){
		warning("These images are mis-registered by more than one cell",immediate.=TRUE)}
	if (cellsize != sp.rast@grid@cellsize[1]){
		warning("The cellsize of", rast.basenm,"= ",sp.rast@grid@cellsize[1]," The cellsize of",rast.basenm1,"= ",cellsize,immediate.=TRUE)}
		
	close(sp.rast)
}}

#print("done checking rasts")


#predcnt <- table(predLUT[,1])[match(rastnm.all,names(table(predLUT[,1])))]	#not redundant-rastnm.all is defined as unique(predLUT[,1])  
													# table() sorts, while unique() preserves order of first occurance
predcnt <- table(predLUT[,1])									
predcnt <- predcnt[match(unique(predLUT[,1]),names(predcnt))]

##define dematrix() function to turn matrix into vector
dematrix <- function(m){m[1:length(m)]}

#################################################################################

while (!final){

	print(paste("numrows =",numrows))
	print(paste("rowcnt =",rowcnt))

	if (rowcnt+numrows >= nrows) {
		numrows <- numrows - ((rowcnt+numrows) - nrows)
		offset <- rowcnt
		final <- TRUE
	}

	preds <- data.frame(matrix(-9999,numrows*ncols,sum(predcnt)))

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

	# change column order from order in rastLUT to order in rastLUT??? redundant!!!
	preds<-preds[,match(predLUT[,2],names(preds))] 
	#preds<-preds[,match(predList,names(preds))] #maybe better? - order same as in predList. never mind, predLUT has already been sorted into same order as predList

	predictions_out<-rep(-9999,nrow(preds))
 
print(paste("NA.ACTION =",NA.ACTION))

	###Give warning of true/interior NA's
	if(any(is.na(preds))){
		N.na<-sum(is.na(preds))
		if(NA.ACTION=="roughfix"){
			warning(paste(N.na,"NA values in predictors repaced with median or most common category"))
		}else{
			if(NA.ACTION=="omit"){
				warning(paste(N.na,"NA values in predictors repaced with  NODATA value of -9999 and no predictions will be made on these pixels"))
			}else{
				stop("1: 'na.action' must be either \"na.roughfix\" or \"na.omit\"")
			}
		}
	} 

	### Turn missing values of factors into NA and give warnings of missing factors
	if(!is.null(model.obj$levels)){
		missing.levels<-model.obj$levels
		for(m in names(model.obj$levels)){
			missing.levels[[m]]<-unique(preds[,m][!preds[,m]%in%model.obj$levels[[m]]])

			if(!length(missing.levels[[m]])==0){
				if(NA.ACTION=="roughfix"){
					warning(paste("categorical factored predictor",m,"contains levels",
						paste(missing.levels[[m]],collapse=", "),  
						"not found in training data, pixels with these levels will be replaced with most common category"))
				}else{
					if(NA.ACTION=="omit"){
						warning(paste("categorical factored predictor",m,"contains levels",
							paste(missing.levels[[m]],collapse=", "),  
							"not found in training data, pixels with these levels will be returned as -9999"))
					}else{
						stop("2:'na.action' must be either \"na.roughfix\" or \"na.omit\"")
					}
				}
			}
			if(!"-9999"%in%model.obj$levels[[m]]){
				preds[,m]<-factor(preds[,m],levels=c(model.obj$levels[[m]],-9999))
			}else{
				preds[,m]<-factor(preds[,m],levels=c(model.obj$levels[[m]]))
			}
		}
	}


	###Actually do NA action
	if(any(is.na(preds))){
		if(NA.ACTION=="roughfix"){
			preds<-na.roughfix(preds)
		}else{
			if(NA.ACTION=="omit"){
				preds[is.na(preds)] <- -9999
			}else{
				stop("3:'na.action' must be either \"na.roughfix\" or \"na.omit\"")
			}
		}
	}


	###make predictions on every row of preds that does not contain any -9999 values
	not9<-!apply(preds==-9999,1,any)
	all9<-!any(not9==TRUE)

	if(!all9){
		preds.not9<-preds[not9,]
		preds.not9<-data.frame(preds.not9)

		for(m in names(model.obj$levels)){
			preds.not9[,m]<-factor(preds.not9[,m],levels=c(model.obj$levels[[m]]))}
		
		##sort col into same order as predList from model.obj. (Started as order in lookup table) 
		preds.not9<-preds.not9[,match(predList,names(preds.not9))] 

		## Model predictions

		print("making predictions")

		if(model.type=="RF"){
			if(response.type=="binary"){
				predictions_out[not9]<-signif(predict(model.obj, preds.not9,type="vote")[,"1"],2)}

			if(response.type=="categorical"){
				Ylev<-levels(model.obj$y)
				PREDCHUNK<-predict(model.obj, preds.not9)
				PREDCHUNK<-factor(PREDCHUNK,levels=Ylev) #this may be redundant-rf appears to do this, but not documented
				if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
					predictions_out[not9]<-unclass(PREDCHUNK)
				}else{
					predictions_out[not9]<-as.numeric(levels(PREDCHUNK))[PREDCHUNK]
				}
			}

			if(response.type=="continuous"){
				predictions_out[not9]<-predict(model.obj, preds.not9)}
		}
		if(model.type=="SGB"){
			predictions_out[not9]<-predict.gbm(	object=model.obj,
								newdata=preds.not9,
								n.trees=n.trees,
								type="response",
								single.tree=FALSE)
		}
	}

	write(predictions_out, file = asciifn, ncolumns=ncols, append=TRUE, sep=" ")
	
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

		write(test_mean,  file = asciifn.mean, ncolumns=ncols, append=TRUE, sep=" ")
		write(test_stdev, file = asciifn.stdev, ncolumns=ncols, append=TRUE, sep=" ")
		write(test_coefv, file = asciifn.coefv, ncolumns=ncols, append=TRUE, sep=" ")}

	rowcnt<-rowcnt + numrows
	offset<-rowcnt	
	}

	###############convert to Imagine Image files#################

	if(make.img){

		###use first predictor layer to define projection###
		RS.pred<-raster(rastnm.all[1])
		CRS<-projection(RS.pred,asText=FALSE)
		if(is.na(CRS@projargs)){
			warning("No projection associated with predictor 1 therefore no projection assigned to output image file")
		}
		###define image file name###
		imgfn<-strsplit(asciifn, split=".txt")[[1]]
		imgfn<-paste(imgfn,"img",sep=".")
		###read ascii grid and write imagine image###
		RS.out<-raster(asciifn,crs=CRS)
		setMinMax(RS.out)
		writeRaster(RS.out,imgfn,overwrite=TRUE,NAflag=-9999)

		if(map.sd && model.type=="RF" && response.type=="continuous"){

			imgfn.mean<-strsplit(asciifn.mean, split=".txt")[[1]]
			imgfn.mean<-paste(imgfn.mean,"img",sep=".")
			RS.out<-raster(asciifn.mean,crs=CRS)
			setMinMax(RS.out)
			writeRaster(RS.out,imgfn.mean,overwrite=TRUE,NAflag=-9999)

			imgfn.stdev<-strsplit(asciifn.stdev, split=".txt")[[1]]
			imgfn.stdev<-paste(imgfn.stdev,"img",sep=".")
			RS.out<-raster(asciifn.stdev,crs=CRS)
			setMinMax(RS.out)
			writeRaster(RS.out,imgfn.stdev,overwrite=TRUE,NAflag=-9999)

			imgfn.coefv<-strsplit(asciifn.coefv, split=".txt")[[1]]
			imgfn.coefv<-paste(imgfn.coefv,"img",sep=".")
			RS.out<-raster(asciifn.coefv,crs=CRS)
			setMinMax(RS.out)
			writeRaster(RS.out,imgfn.coefv,overwrite=TRUE,NAflag=-9999)
		}
	}
}

