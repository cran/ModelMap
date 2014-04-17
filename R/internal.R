
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
						prediction.type,
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

	CMXfn<-paste(MODELpredfn,"_cmx.csv",sep="")

	LEVELS.pred <- levels(PRED$pred)
	LEVELS.obs  <- levels(PRED$obs)

	print(paste("LEVELS.pred: ",paste(LEVELS.pred, collapse=" ")))
	print(paste("LEVELS.obs:  ",paste(LEVELS.obs,  collapse=" ")))


	LEVELS<-unique(c(LEVELS.pred,LEVELS.obs))

	CMX<-table(	predicted = factor(PRED$pred,levels=LEVELS),
			observed = factor(PRED$obs,levels=LEVELS))


	CMX.out<-matrix("cmx",nrow(CMX)+5,ncol(CMX)+4)
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
	CMX.out[nrow(CMX.out)-2,1:2]<-"total"
	
	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)-1]<-apply(CMX,1,sum)
	CMX.out[nrow(CMX.out)-2,(1:ncol(CMX))+2]<-apply(CMX,2,sum)
	CMX.out[nrow(CMX.out)-2,ncol(CMX.out)-1]<-sum(CMX)

	###marginals###
	CMX.diag<-diag(CMX)

	CMX.out[1:2,ncol(CMX.out)]<-"Commission"
	CMX.out[nrow(CMX.out)-1,1:2]<-"Omission"

	CMX.out[(1:nrow(CMX))+2,ncol(CMX.out)]<-1-(CMX.diag/apply(CMX,1,sum))
	CMX.out[nrow(CMX.out)-1,(1:ncol(CMX))+2]<-1-(CMX.diag/apply(CMX,2,sum))

	###pcc###
	CMX.out[nrow(CMX.out)-2,ncol(CMX.out)]<-"PCC"
	CMX.out[nrow(CMX.out)-1,ncol(CMX.out)-1]<-"PCC"

	CMX.out[nrow(CMX.out)-1,ncol(CMX.out)]<-sum(CMX.diag)/sum(CMX)

	###MAUC###

	#if(prediction.type=="CV"){
      #      PRED.mauc = PRED[4:(ncol(PRED)-1)]
	#}else{
	#	PRED.mauc = PRED[,4:ncol(PRED)]
	#}

	PRED.mauc <- PRED[,LEVELS.pred]

	if(any(!LEVELS.obs%in%LEVELS.pred)){

		LEVELS.new <- LEVELS.obs[!LEVELS.obs%in%LEVELS.pred]

		LEVELS.new.paste <-paste(LEVELS.new,collapse=" ")

		warning("Response categories: ", LEVELS.new.paste, " observed in test data were not included in training data")
		
		PRED.new  <- matrix(NA,nrow=nrow(PRED.mauc),ncol=length(LEVELS.new))
		colnames(PRED.new)<-LEVELS.new

		PRED.mauc<-cbind(PRED.mauc,PRED.new)
	}

	VOTE <- multcap(  response = PRED$obs,
                		predicted= as.matrix(PRED.mauc) )

	MAUC  <- HandTill2001::auc(VOTE)

	CMX.out[nrow(CMX.out),1]<-"MAUC"
	CMX.out[nrow(CMX.out),2]<-MAUC

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

		
		###Importance Plot###
		#print(paste("IMPORTANCEfn =",IMPORTANCEfn)

		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(IMPORTANCEfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(IMPORTANCEfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(IMPORTANCEfn,".pdf",sep=""),width = device.width, height = device.height)}
 
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

		
	
		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(IMPORTANCEfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(IMPORTANCEfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(IMPORTANCEfn,".pdf",sep=""),width = device.width, height = device.height)}


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
	
	
		if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
		if(device.type[i]=="jpeg"){jpeg(filename=paste(THRESHOLDPLOTSfn,".jpg",sep=""),width = device.width, height = device.height, res=jpeg.res, units="in")}
		if(device.type[i]=="postscript"){postscript(file=paste(THRESHOLDPLOTSfn,".ps",sep=""),width = device.width, height = device.height)}
		if(device.type[i]=="pdf"){pdf(file=paste(THRESHOLDPLOTSfn,".pdf",sep=""),width = device.width, height = device.height)}


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
				nTrain = NULL,
				#train.fraction = NULL,       	# fraction of data for training,
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
		warning("keep.data reset to TRUE because data needed for gbm.more() function needed for OOB determination of optimal number of trees")}

	SGB <- gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=100,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction,   
				nTrain = nTrain,       	
				#train.fraction = train.fraction,       	           		
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
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package however OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive however using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap")

}else{

	SGB <- gbm.fit(	x=qdata.x,
				y=qdata.y,        
				distribution=distribution,
				n.trees=n.trees,                	
				shrinkage=shrinkage, 
				interaction.depth=interaction.depth,		
				bag.fraction = bag.fraction, 
				nTrain = nTrain,         	
				#train.fraction = train.fraction,       	           		
				n.minobsinnode = n.minobsinnode,
				keep.data=keep.data,
				var.monotone=var.monotone)

	if(!is.null(nTrain)){
		if(nTrain<nrow(qdata.x)){
			SGB$best.iter <- suppressWarnings(gbm.perf(SGB,method="test",plot.it=FALSE))
			if(SGB$best.iter>0.9*n.trees){
				warning("best number of trees is ", SGB$best.iter, " and total number trees tested was ", n.trees, " therefore you may want to explore increasing the 'n.trees' argument")
			}
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
				proximity = NULL,
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

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

#print("finished tuning")

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
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
					proximity = NULL,
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

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

#print("finished tuning")

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
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
					proximity = NULL,
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

	A<-list(	x=quote(qdata.x), y=quote(qdata.y),
			doBest=FALSE,
			importance=TRUE,
			proximity=FALSE,
			plot=FALSE,
			replace=replace,
			strata=strata,
			sampsize=sampsize)

	A<-A[!sapply(A, is.null)]

	RT<-do.call("tuneRF", A)

	mtry<-RT[which.min(RT[,2]),1]
}

A<-list(	x=quote(qdata.x), y=quote(qdata.y),
		importance=TRUE,
		proximity=proximity,
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
				proximity=proximity,
			
			# SGB arguments:
				n.trees=NULL,                 # number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				nTrain = NULL,
				#train.fraction = NULL,       	# fraction of data for training,
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
						proximity=proximity,
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
							proximity=proximity,
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
							proximity=proximity,
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
					nTrain=nTrain,
					#train.fraction=train.fraction,      # fraction of data for training,
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


	#create qdata.y that has been transformed as response gets transformed

	qdata.y<-qdata[,response.name]
	if(response.type=="binary"){
		if(!is.numeric(qdata.y)){
			stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
		qdata.y[qdata.y>0]<-1
		qdata.y[qdata.y<0]<-0
	}
	if(response.type=="categorical"){
		if(!is.factor(qdata.y)){qdata.y<-as.factor(qdata.y)}
	}	



	#If prediction is OOB check that training data responses are the same as those used to build the model
	if(prediction.type=="OOB"){
		if(!isTRUE(all.equal(PRED$obs,qdata.y))){
			stop("Prediction type is 'OOB' but responses in 'qdata.trainfn' do not match data used to build the model")
		}
	}

	#Note - If this warning appears there is a bug
	if(!isTRUE(all.equal(PRED$obs,qdata.y))){
		stop("Something has gone wrong and observed responses are scrambled")
	}

	rownames(PRED)<-rownames(qdata)

}else{
	print(paste("Begining ",v.fold,"-fold cross validation:",sep=""))

	if(model.type=="RF"){
		ntree<-model.obj$ntree
		mtry<-model.obj$mtry
		replace<-model.obj$call$replace
		predList<-row.names(model.obj$importance)}
	if(model.type=="SGB"){
		shrinkage<-model.obj$shrinkage
		interaction.depth<-model.obj$interaction.depth
		bag.fraction<-model.obj$bag.fraction
		nTrain<-model.obj$nTrain
		n.minobsinnode<-model.obj$n.minobsinnode
		predList<-model.obj$var.names

		#deal with nTrain<nrow(qdata)
		print(paste("nTrain =",nTrain))
		print(paste("nrow(qdata) =",nrow(qdata)))
		if(nTrain<nrow(qdata)){
			qdata<-qdata[1:nTrain,]}	
	}

	n.data=nrow(qdata)
	n.per.fold<-floor(n.data/v.fold)
	cv.index<-sample(rep(1:v.fold,(n.per.fold+1))[1:n.data])
	
	PRED<-data.frame(matrix(0,0,2))
	names(PRED)<-c("obs","pred")

	###start folds###
	for(i in 1:v.fold){
		print(paste("starting fold", i))
		train.cv<-(1:nrow(qdata))[cv.index!=i]
		qdata.train.cv<-qdata[train.cv,]
		qdata.test.cv<-qdata[-train.cv,]

		###check for factor levels not found in other folds###
		print("starting factor checks")
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


				print(paste("     predictor            =",p))
				print(paste("     train levels[p]      =",paste(train.levels,collapse=",")))
				print(paste("     test levels[p]       =",paste(test.levels,collapse=",")))
				#print(paste("     invalid              =",paste(invalid,collapse=",")))
				#print(paste("     invalid levels[[p]]  =",paste(invalid.levels[[p]],collapse=",")))

				#print(paste("          invalid length =",length(invalid.levels[[p]])))
				
				#print(paste("     levels qdata.test.cv =",paste(levels(qdata.test.cv[,p]),collapse=",")))

				if(length(invalid.levels[[p]])>0){
					print("     NA action triggered")
					N.invalid<-sum(is.na(qdata.test.cv[,p]))
					warn.lev<-paste(invalid.levels[[p]],collapse=", ")
					warning(	"Factored predictor ",p," in fold ",i,
							" contains ", N.invalid, " data points of levels ",warn.lev,
							" not found in other folds and these levels treated as NA")
				}
			}
			
			#print(qdata.test.cv[,predList])
			NA.pred<-apply(qdata.test.cv[,predList],1,function(x){any(is.na(x))})
			#print(paste(     "NA.pred =",any(NA.pred)))
			#print("TEST TEST TEST")
			if(any(NA.pred)){
				print("about to treat na's")
				print("NA.ACTION:")
				print(NA.ACTION)
				#print("na.action:")
				#print(na.action)
				qdata.test.cv<-na.action(qdata.test.cv)
				print("finished treating NA's")}
			#print(qdata.test.cv[,predList])
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
			print(paste("      making SGB model for fold",i))

			
			SGB.cv<-model.SGB(	qdata=qdata.train.cv,
							predList=predList,
							response.name=response.name,
							seed=NULL,
							response.type=response.type, 				
							n.trees=n.trees,	               	# number of trees
							shrinkage=shrinkage,   	      	# shrinkage or learning rate,
                  				interaction.depth=interaction.depth,# 1: additive model, 2: two-way interactions, etc.
							bag.fraction=bag.fraction,          # subsampling fraction, 0.5 is probably best
							nTrain=nrow(qdata.train.cv),
							#train.fraction=train.fraction,      # fraction of data for training,
                  				n.minobsinnode=n.minobsinnode       # minimum total weight needed in each node
							)
			print(paste("      making SGB predictions for fold",i))
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
###################### Production Prediction - Sub Functions ################################
#############################################################################################
#############################################################################################



#############################################################################################
################################## Check filename ###########################################
#############################################################################################


FNcheck<-function(OUTPUTfn,folder,ERROR.NAME){
### checks filename for valid extensions, stops if there is more than 1 '.' in filename or invalid extension after the dot.
### If no, extension, warns and adds default extension of .img
### If no path, adds path from 'folder'.

	validExtensions<-c(".tif",
				".tiff",
				".grd",
				".asc",
				".nc",
				".cdf",
				".ncdf",
				".kml",
				".kmz",
				".big",
				".sgrd",
				".sdat",
				".bil",
				".bsq",
				".bip",
				".bmp",
				".gen",
				".bt",
				".envi",
				".ers",
				".img",
				".rst",
				".mpr",
				".rsw") 


	OUTPUTsplit<-strsplit(basename(OUTPUTfn),split="\\.")[[1]]
	OUTPUText<-paste(".",tail(OUTPUTsplit, 1),sep="")

	#if basename of output filename has more than 1 '.'
	if(length(OUTPUTsplit) > 2){
		stop(ERROR.NAME," ",OUTPUTfn," has more than one '.' The ",ERROR.NAME," should only use the character '.' to indicate the file extension.")}

	#if basename has exactly one dot, check if extension is valid
	if(length(OUTPUTsplit) == 2){
		if(!tolower(OUTPUText)%in%validExtensions){
			#OUTPUTfn<-paste(dirname(OUTPUTfn),"/",OUTPUTsplit[1],".img",sep="")
			stop(ERROR.NAME," extension ",OUTPUText," is invalid. See 'writeFormats()' for a list of valid raster file types.")}}

	#if basename has no extension, add default extension of '.img'
	if(length(OUTPUTsplit) == 1){
		OUTPUTfn<-paste(OUTPUTfn,".img",sep="")
		warning(ERROR.NAME," ",OUTPUTsplit," does not include an extension, using default extension of '.img'. New ",ERROR.NAME," is ",OUTPUTfn)
	}

	# if OUTPUTfn has no directory path, add folder to name
	if(identical(basename(OUTPUTfn),OUTPUTfn))
		{OUTPUTfn<-file.path(folder,OUTPUTfn)}

	return(OUTPUTfn)
}


#############################################################################################
################################# Fix projections ###########################################
#############################################################################################


projfix<-function(RAST,OUTPUTfn.noext){
#RAST	list of rasters suitable to be given to stack() function

	#find all projections
	PROJ<-lapply(RAST,projection)

	#create output file of all projections
	OUTPUTfn.proj<-paste(OUTPUTfn.noext,"_projections.txt",sep="")
	PROJout<-data.frame(predictor=names(PROJ),projection=sapply(PROJ,I))
	write.table(PROJout,file=OUTPUTfn.proj,row.names=FALSE,sep="\t")

	#find source filenames
	FN<-sapply(RAST,function(x){basename(filename(x))})

	#check if all rasters have projections
	HAVEPROJ<-sapply(RAST,function(x){is.na(projection(x)) || is.null(projection(x))})
	if(any(HAVEPROJ)){
		if(sum(HAVEPROJ)==1){
			warning.message<-paste("raster for predictor ",names(HAVEPROJ[HAVEPROJ])," is missing a projection, see '", OUTPUTfn.proj, "' for details", sep="")
		}else{
			warning.message<-paste("rasters for predictors ",paste(names(HAVEPROJ[HAVEPROJ]),collapse=" ")," are missing projections, see '", OUTPUTfn.proj, "' for details", sep="")
		}
		stop(warning.message)
	}


#	#check if projections are similar enough to reconcile
#	SAME<-sapply(PROJ,function(x,proj){projcompare(x,proj)},proj=PROJ[[1]])
#	SAME.R<-sapply(RAST,function(x,control.rast){compareRaster(x,control.rast,stopiffalse=FALSE)},control.rast=RAST[[1]])
#
#	if(all(SAME)){
#		if(all(SAME.R)){
#			#they all match, no need to do anything
#			RAST.out<-RAST
#		}else{
#			#set all projections to match layer 1
#			RAST.out<-RAST
#			RAST.out<-lapply(RAST,function(x,rast){projection(x)<-projection(rast); x},rast=RAST[[1]])
#			warning("one or more predictor layers had projections that needed to be reconciled, see '", OUTPUTfn.proj, "' for details")
#		}
#	}else{
#		stop("projections of predictor layers too different to reconcile, see '", OUTPUTfn.proj, "' for details")
#	}

	#check if projections are the same
	#SAME<-sapply(PROJ,function(x,proj){projcompare(x,proj)},proj=PROJ[[1]])
	SAME.R<-sapply(RAST,function(x,control.rast){compareRaster(x,control.rast,stopiffalse=FALSE)},control.rast=RAST[[1]])

	#if(all(SAME)){
	if(all(SAME.R)){
			#they all match, no need to do anything
			print("all predictor layer rasters match")
			RAST.out<-RAST
		#}else{
		#	#set all projections to match layer 1
		#	RAST.out<-RAST
		#	RAST.out<-lapply(RAST,function(x,rast){projection(x)<-projection(rast); x},rast=RAST[[1]])
		#	warning("one or more predictor layers had projections that needed to be reconciled, see '", OUTPUTfn.proj, "' for details")
		#}
	}else{
		stop("predictor layer rasters too different to reconcile, see '", OUTPUTfn.proj, "' for details")
	}

	return(RAST.out)
}

#############################################################################################
###################################### grd to gri ###########################################
#############################################################################################

grd2gri<-function(x){
	paste(strsplit(x,".grd")[[1]],".gri",sep="")}

#############################################################################################
###################################### Compare projections ##################################
#############################################################################################
# Note this function is not currently in use, it is currently replaced by raster package function compareRaster()

	projcompare <- function(layer1prj, layer2prj){
		########################################################################################
		## DESCRIPTION: Internal function to compare projections. If not the same, return FALSE
		##
		## ARGUMENTS:
		##   layer1prj - the projection name for layer1
		##   layer2prj - the projection name for layer1	
		##
		## NOTE:
		## Use function proj4string(layer) or projection(layer) to get the projection name 
		########################################################################################


		ellpsna <- FALSE
		datumna <- FALSE

		if(is.null(layer1prj) | is.null(layer2prj)){
			stop("check projection layers")
		}
	
		projnm1 <- strsplit(strsplit(layer1prj, "+proj=")[[1]][2], " +")[[1]][1]
		ellpsnm1 <- strsplit(strsplit(layer1prj, "+ellps=")[[1]][2], " +")[[1]][1]
		datumnm1 <- strsplit(strsplit(layer1prj, "+datum=")[[1]][2], " +")[[1]][1]
	
		projnm2 <- strsplit(strsplit(layer2prj, "+proj=")[[1]][2], " +")[[1]][1]
		ellpsnm2 <- strsplit(strsplit(layer2prj, "+ellps=")[[1]][2], " +")[[1]][1]
		datumnm2 <- strsplit(strsplit(layer2prj, "+datum=")[[1]][2], " +")[[1]][1]

		if(is.na(ellpsnm1) | is.na(ellpsnm2)){
			ellpsna <- TRUE
		}
		if(is.na(datumnm1) | is.na(datumnm2)){
			datumna <- TRUE
		}

		if(is.null(ellpsna) & is.null(datumna)){
			stop("CHECK PROJECTIONS. NEED ELLPSNM OR DATUM DEFINED IN BOTH PROJECTIONS")
		}else{
	
			## ASSUME WGS84 and NAD83 ARE EQUAL
			if(!is.na(datumnm1)){ if(datumnm1 == "WGS84"){ datumnm1 = "NAD83" } }
			if(!is.na(datumnm2)){ if(datumnm2 == "WGS84"){ datumnm2 = "NAD83" } }				

			if(projnm1 == projnm2){
				if(ellpsnm1 == ellpsnm2 | ellpsna){
					if(datumnm1 == datumnm2 | datumna){
						projmatch <- TRUE
					}else{
						projmatch <- FALSE
					}
				}else{
					projmatch <- FALSE
				}
			}else{
				projmatch <- FALSE
			}
		}
		return(projmatch)
	}

#############################################################################################
#############################################################################################
############################## color translation ############################################
#############################################################################################
#############################################################################################

col2trans<-function(col.names,alpha=0.5){
	col.out <- rgb(t(col2rgb(col.names)/255),alpha=alpha)
	return(col.out)
}

#############################################################################################
#############################################################################################
#################### Production Prediction - Actual Function ################################
#############################################################################################
#############################################################################################

production.prediction<-function(	model.obj,
						model.type,
						rastLUT,
						#na.action=NULL,
						NA.ACTION=NULL,
						response.type,
						keep.predictor.brick,							
						map.sd=FALSE,
						OUTPUTfn,
						OUTPUTfn.noext,
						#OUTPUTpath,
						OUTPUTname,
						OUTPUText,
						n.trees){

#############################################################################################
################################## Create File names ########################################
#############################################################################################

### Creat filename for native raster format map output

TMPfn.map <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_map_",sep=""))

### Creat filename for predictor raster brick

if(keep.predictor.brick){
	OUTPUTfn.brick <- paste(OUTPUTfn.noext,"_brick",sep="")
}else{
	OUTPUTfn.brick <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_brick_",sep=""))
}

### map sd filenames

if(map.sd && model.type=="RF" && response.type=="continuous"){

	OUTPUTfn.mean  <- paste(OUTPUTfn.noext,"_mean",OUTPUText,sep="")
	OUTPUTfn.stdev <- paste(OUTPUTfn.noext,"_stdev",OUTPUText,sep="")
	OUTPUTfn.coefv <- paste(OUTPUTfn.noext,"_coefv",OUTPUText,sep="")

	TMPfn.mean     <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_mean_",sep=""))
	TMPfn.stdev    <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_stdev_",sep=""))
	TMPfn.coefv    <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_coefv_",sep=""))	

}else{
	map.sd<-FALSE
}

#####################################################################################
########################## Extract predictor names ##################################
#####################################################################################


## Make sure raster predictor names match names in training/test data.

if(model.type=="RF"){
	predList<-row.names(model.obj$importance)}
if(model.type=="SGB"){
	predList<-model.obj$var.names}

predLUT<-rastLUT[match(predList, rastLUT[,2]),] #only the predictors from the model, in same order as in model
if(any(predList!=predLUT[,2])){
	stop("predictor names from model do not match short names in rastLUT")}	

anyFactor<-FALSE
if(!is.null(model.obj$levels)){
	anyFactor<-TRUE
	predFactor<-names(model.obj$levels)
	extraLevels<-vector("list",length(predFactor))
	names(extraLevels)<-predFactor
}

#####################################################################################
########################## Extract response classes #################################
#####################################################################################

if(response.type=="categorical"){
	Rclasses<-model.obj$classes
	Nrc<-length(Rclasses)
}


########################################################################################
################################ Set Data Type #########################################
########################################################################################

if(response.type=="binary"){data.type <- "FLT4S"}
if(response.type=="categorical"){
	data.type <- "INT2S"
	Ylev<-levels(model.obj$y)}
if(response.type=="continuous"){data.type <- "FLT4S"}

########################################################################################
########################### Build raster stack #########################################
########################################################################################

#require(raster)

RAST<-vector("list", 0)

for(p in predList){
	rastfn<-rastLUT[rastLUT[,2]==p,1]
	band<-  rastLUT[rastLUT[,2]==p,3]

	RAST[[p]]<-raster(rastfn,band=band)
}

RAST<-projfix(RAST,OUTPUTfn.noext=OUTPUTfn.noext)

#compareRaster(RAST,stopiffalse=TRUE, showwarning=TRUE) #, crs=FALSE)

RS<-stack(RAST)
RB<-brick(RS,values=TRUE,filename=OUTPUTfn.brick,overwrite=TRUE)

print("brick done")
########################################################################################
############################# Loop through rows ########################################
########################################################################################

out <- raster(RAST[[1]])
dataType(out) <- data.type
NAvalue(out) <- -9999

print("write start")
print(OUTPUTfn)
print(paste("datatype =",data.type))

out <- writeStart(out, filename=TMPfn.map, overwrite=TRUE, datatype=data.type)
#out <- writeStart(out, overwrite=TRUE,  dataType=data.type) #doesn't work


if(map.sd){
	out.mean <- raster(RAST[[1]])
	out.mean <- writeStart(out.mean, filename=extension(TMPfn.mean,""), overwrite=TRUE, datatype=data.type)

	out.stdev <- raster(RAST[[1]])
	out.stdev <- writeStart(out.stdev, filename=extension(TMPfn.stdev,""), overwrite=TRUE, datatype=data.type)

	out.coefv <- raster(RAST[[1]])
	out.coefv <- writeStart(out.coefv, filename=extension(TMPfn.coefv,""), overwrite=TRUE, datatype=data.type)
}

###############################################################################################################
print("starting loops")
for(r in 1:(dim(RB)[1])){
	print(paste("rows =",r))

	v <- data.frame(getValues(RB, r))

	### deal with factored predictors ###
	if(anyFactor){
		for(f in predFactor){
			f.lev<-model.obj$levels[[f]]
			f.extra<-v[,f][!v[,f]%in%f.lev]
			extraLevels[[f]]<-unique(c(extraLevels[[f]],f.extra))
			
			#v[,f][!v[,f]%in%f.lev] <- ??? #leaving this line just in case we decide to add rough fix back in

			v[,f]<-factor(v[,f],levels=f.lev)
		}
	}

	### deal with -9999 ###
	nonPredict <- apply(((v == -9999)|is.na(v)), 1, any)

	v.pred<-rep(NA,length=nrow(v))


	if(any(!nonPredict)){
		if(model.type=="RF"){
			if(response.type=="binary"){
				v.pred[!nonPredict] <- signif(predict(model.obj, v[!nonPredict,],type="vote")[,"1"],2)}

			if(response.type=="categorical"){
				PRED <- predict(model.obj, v[!nonPredict,])
				
				if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
					v.pred[!nonPredict] <- as.integer(PRED)
				}else{
					v.pred[!nonPredict] <- as.integer(as.character(PRED))
				}
	
			}

			if(response.type=="continuous"){
				v.pred[!nonPredict] <- predict(model.obj, v[!nonPredict,])}
		}
		if(model.type=="SGB"){
			v.pred[!nonPredict] <- predict.gbm(	object=model.obj,
									newdata=v[!nonPredict,],
									n.trees=n.trees,
									type="response",
									single.tree=FALSE)
		}
	}

	writeValues(out, v.pred, r)

	if(map.sd){
		v.mean  <- rep(NA,length=nrow(v))
		v.stdev <- rep(NA,length=nrow(v))
		v.coefv <- rep(NA,length=nrow(v))

		if(any(!nonPredict)){
			v.everytree <- predict(model.obj, v[!nonPredict,], predict.all=TRUE)$individual #only has rows for !nonPredict
			v.mean[!nonPredict]  <- apply(v.everytree,1,mean)
			v.stdev[!nonPredict] <- apply(v.everytree,1,sd)
			v.coefv[!nonPredict] <- v.stdev[!nonPredict]/v.mean[!nonPredict]
			v.coefv[!nonPredict][v.stdev[!nonPredict]==0]<-0
			v.coefv[!nonPredict][v.mean[!nonPredict]==0]<-0
			rm(v.everytree)
		}

		writeValues(out.mean, v.mean, r)
		writeValues(out.stdev, v.stdev, r)
		writeValues(out.coefv, v.coefv, r)
	}

}

################################################################################################

out <- writeStop(out)

out<-setMinMax(out)
writeRaster(out,OUTPUTfn,overwrite=TRUE, datatype=data.type)

if(map.sd){
	out.mean  <- writeStop(out.mean)
	out.stdev <- writeStop(out.stdev)
	out.coefv <- writeStop(out.coefv)
	
	out.mean  <- setMinMax(out.mean)
	writeRaster(out.mean, OUTPUTfn.mean,overwrite=TRUE,datatype=data.type)

	out.stdev <- setMinMax(out.stdev)
	writeRaster(out.stdev,OUTPUTfn.stdev,overwrite=TRUE,datatype=data.type)

	out.coefv <- setMinMax(out.coefv)
	writeRaster(out.coefv,OUTPUTfn.coefv,overwrite=TRUE,datatype=data.type)
}

###clean up tmp files

#print(OUTPUTfn.brick)
#print(TMPfn.map)

file.remove(TMPfn.map)
file.remove(grd2gri(TMPfn.map))

if(!keep.predictor.brick){
	file.remove(OUTPUTfn.brick)
	file.remove(grd2gri(OUTPUTfn.brick))
}
if(map.sd){
	file.remove(TMPfn.mean)
	file.remove(TMPfn.stdev)
	file.remove(TMPfn.coefv)
	file.remove(grd2gri(TMPfn.mean))
	file.remove(grd2gri(TMPfn.stdev))
	file.remove(grd2gri(TMPfn.coefv))
}

}

