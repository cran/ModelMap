


#############################################################################################
#############################################################################################
######################## Model Run Function == "THE BUTTON" #################################
#############################################################################################
#############################################################################################



model.map<-function(	model.obj=NULL,
				model.type=NULL,	# "RF", "GAM", "SGB"
				qdata.trainfn=NULL,
				qdata.testfn=NULL,
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				MODELfn=NULL,
				rastLUT=NULL,
				rastLUTfn=NULL,
				rastnmVector=NULL,
				predList=NULL,
				predFactor=FALSE,
				response.name=NULL,
				response.type=NULL,		# "binary", "continuous",
				unique.rowname=NULL,	# Row identifier
				seed=NULL,
			# Model Evaluation Arguments
				predict=NULL,		# Predicts over test set if available, otherwise training set
				MODELpredfn=NULL,
				na.action="na.omit",	# also used for mapping
				v.fold=NULL,
				diagnostics=predict,	# Diagnostic graphs of predictions (above). Predict must be TRUE
				device.type=NULL,	# options: "windows", "jpeg", "postscript", "win.metafile"
				DIAGNOSTICfn=NULL,
				jpeg.res=72,
				device.width=7,
				device.height=7,
				cex=par()$cex,
				req.sens,req.spec,FPC,FNC,
			# RF arguments:
				ntree=500,
				mtry=NULL,
			# SGB arguments:
				n.trees=NULL,                 	# number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  		interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				train.fraction = 1.0,       	# fraction of data for training,
                  		n.minobsinnode = 10,         	# minimum total weight needed in each node
			# Mapping arguments
				map=NULL,
				numrows = 500,		# Only used if mapping. Number of rows predicted at a time.				
				map.sd=FALSE,		# Only used if mapping. Generates 3 additional maps (mean,sd,coefvar)
				asciifn=NULL,
				asciifn.mean=NULL,
				asciifn.stdev=NULL,
				asciifn.coefv=NULL
){

## Note: Must have R version 2.5.1 or greater
## Note: Must have all rasters in same size, same resolution, and with -9999 for nodata values.
## Note: The names of the predictors must be the same as the predictor names in the training/test data.
##		If the predictors are part of a layer stack, the names must include the number of the
##		layer within the stack (Example: stacknameb1, stacknameb2, stacknameb3, ...)
## Note: If subsampled test data is desired, run training selection.R first.

## Note: Must have the following packages installed: rgdal, gbm, randomForest

# Note:	'rgdal' is only used to make map
#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.




#############################################################################################
################################### Check Platform ##########################################
#############################################################################################


Rplatform<-.Platform$OS.type


#############################################################################################
################################### Add to Filters Table ####################################
#############################################################################################


## Adds to file filters to Cran R Filters table.
if(.Platform$OS.type=="windows"){
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))}


#############################################################################################
########### Should validation predictions and diagnostics be made ###########################
#############################################################################################

print("About to check if should make validation predictions")

if(is.null(predict)){
	predict <- select.list(c("YES","NO"), title="Validation Predictions?")
	predict <- switch(predict,"YES"=TRUE,"NO"=FALSE,FALSE)

	if(predict==TRUE){
		diagnostics <- select.list(c("YES","NO"), title="Validation Diagnostics?")
		diagnostics <- switch(diagnostics,"YES"=TRUE,"NO"=FALSE,FALSE)
		if(diagnostics==TRUE){
			if(is.null(device.type)){
				device.type <- select.list(c("default","jpeg","pdf","postscript","win.metafile"), title="Diagnostic Output?", multiple = TRUE)
				device.type <- c(device.type,"default")
			}
			if(length(device.type)==0 || is.null(device.type)){
				device.type <- "default"
			}
		}
	}else{diagnostics<-FALSE}
}else{
	if(predict==TRUE && diagnostics==TRUE){
		if(is.null(device.type)){device.type<-"default"}
	}
}

if(!is.null(device.type)){
	device.type[device.type=="windows"]<-"default"
	if(any(!device.type%in%c("default","jpeg","pdf","postscript","win.metafile"))){
		stop("Illegal 'device.type'. Device types must be one or more of 'default', 'jpeg', 'pdf', 'postscript', or 'win.metafile'")
	}
	device.type<-sort(device.type)
	if("default"%in%device.type){
		device.type<-c(device.type[device.type!="default"],"default")
	}
}
#############################################################################################
################################ Select Output Folder #######################################
#############################################################################################

if(is.null(folder)){
	if(.Platform$OS.type=="windows"){
		folder<-choose.dir(default=getwd(), caption="Select directory")
	}else{
		folder<-getwd()}
}
	
#############################################################################################
###################################### Load Data ############################################
#############################################################################################

if(is.null(model.obj) || predict==TRUE){
	## If training data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.trainfn)){
		if(.Platform$OS.type=="windows"){
			qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
			if(is.null(qdata.trainfn)){stop("")}
		}else{stop("To create a model or make validation predictions, you must provide qdata.trainfn")}
	}

	## Check if file name is full path or basename
	if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
		if(identical(basename(qdata.trainfn),qdata.trainfn)){qdata.trainfn<-paste(folder,"/",qdata.trainfn,sep="")}
	}

	## Read in training data
	if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
		qdata.train<-qdata.trainfn
	}else{
		qdata.train<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)}

	## If test data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.testfn)){
		testfn.ask <- select.list(c("NO","YES"), title="Do you have a test file?")
		if (testfn.ask == "YES"){
			if(.Platform$OS.type=="windows"){
				qdata.testfn <- choose.files(caption="Select validation file", filters = Filters["csv",], multi = FALSE)
				if(is.null(qdata.testfn)){stop("")}
			}else{stop("If you have a separate test data set you must provide qdata.testfn")}
		}else{
			qdata.testfn <- FALSE}
	}		

	## Check qdata.testfn for format:
	##	if it is already a dataframe or matrix
	##	if it is FALSE (i.e. no test set)
	##	if file name is full path or basename
	##	if column headers do not match training data

	if(is.matrix(qdata.testfn)==TRUE || is.data.frame(qdata.testfn)==TRUE){
		qdata.test<-qdata.testfn
		if(!identical(names(qdata.test),names(qdata.train))){
			stop("Column names of training and test data must be the same")}
		qdata<-rbind(qdata.train,qdata.test)
		qdata<-data.frame(qdata)
		train<-1:nrow(qdata.train)
	}else{
		if (qdata.testfn == FALSE){
			qdata<-qdata.train
			qdata<-data.frame(qdata)
			train<-1:nrow(qdata)
		}else{
			if(identical(basename(qdata.testfn),qdata.testfn)){qdata.testfn<-paste(folder,"/",qdata.testfn,sep="")}

			qdata.test<-read.table(file=qdata.testfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)

			if(!identical(names(qdata.test),names(qdata.train))){
				stop("Column names of training and test data must be the same")}
			qdata<-rbind(qdata.train,qdata.test)
			qdata<-data.frame(qdata)
			train<-1:nrow(qdata.train)
		}
	}
}

#############################################################################################
############### Should cross - validation predictions be made ###############################
#############################################################################################

if(predict==TRUE){
	if(is.matrix(qdata.testfn)!=TRUE && is.data.frame(qdata.testfn)!=TRUE){
		if(qdata.testfn == FALSE){
			if(is.null(v.fold)){
				v.fold.ask <- select.list(c("NO","YES"), title="10-fold cross validation?")
				v.fold<-switch(v.fold.ask,"YES"=10,"NO"=FALSE,FALSE)
			}
		}else{
			v.fold<-FALSE
		}
	}
}

#############################################################################################
############################## Should map be made ###########################################
#############################################################################################

if(is.null(map)){
	map <- select.list(c("NO","YES"), title="Make Map?")
	map <- switch(map,"YES"=TRUE,"NO"=FALSE,FALSE)}


#####################################################################################
##################### Extract Model Type from model.obj #############################
#####################################################################################

if(is.null(model.obj)){
	if(is.null(model.type)){
		model.type <- select.list(c("RF","SGB"), title="Select model type.")}
	if(model.type=="" || is.null(model.type)){
		stop("model.type is required")}
}

if(!is.null(model.obj)){
	model.type.supplied<-model.type
	model.type.long<-attr(model.obj,"class")
	model.type<-switch(model.type.long,	"randomForest"="RF",
							"gbm"="SGB",
							"unknown")
	if(model.type=="unknown"){
		stop("model.obj is of unknown type")}
	if(!is.null(model.type.supplied)){
		if(model.type!=model.type.supplied){
			warning("model.obj is a model of type",model.type,"not the supplied type of",model.type.supplied)
		}
	}
}	
print(paste("model.type =",model.type))

if (model.type == "SGB") {
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package. OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap.")
}



#############################################################################################
#################################### Pick Response ##########################################
#############################################################################################

if(!is.null(model.obj)){
	if(!is.null(model.obj$response)){
		model.response.name<-model.obj$response
		if(!is.null(response.name)){
			if(response.name!=model.response.name){
				print(paste("supplied model.obj has response",model.response.name,"not",response.name))}}
		response.name<-model.response.name
	}
}

if(is.null(model.obj) || predict==TRUE){

## If the response variable is NULL, then the user selects variable from pop-up list.
if (is.null(response.name)){
	response.name <- select.list(names(qdata), title="Select response name.")
	if(response.name=="" || is.null(response.name)){
		stop("response.name is required")}	
}

print(paste("response.name =",response.name))

## If model.obj null, ask for response.type

if(is.null(model.obj)){
	if(is.null(response.type)){
		response.type <- select.list(c("continuous","binary"), title="Select response type.")}
	if(response.type=="" || is.null(response.type)){
		stop("response.type is required")}
	if(response.type=="categorical"){response.type<-"binary"}	
}


## If model.obj not NULL, then extract response.type from model.obj

if(!is.null(model.obj)){
	if(model.type=="RF"){
		model.response.type<-switch(model.obj$type,"regression"="continuous","classification"="binary","unknown")}
	if(model.type=="SGB"){
		model.response.type<-switch(model.obj$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")}
	if(model.response.type=="unknown"){stop("supplied model.obj has an unknown response type")}
	if(!is.null(response.type)){
		if(response.type!=model.response.type){
			print(paste("supplied model.obj response is of type",model.response.type,"not",response.type))}}
	response.type<-model.response.type
}
}
#############################################################################################
############################ unique row identifier ##########################################
#############################################################################################

if(predict==TRUE){
	if (is.null(unique.rowname)){
		unique.rowname <- select.list(c(names(qdata),"row_index"), title="Select unique row identifier.")	
		if(unique.rowname!="row_index" && unique.rowname!=""){
			rownames(qdata)<-qdata[,unique.rowname]
		}	
	}else{
		if(!(unique.rowname%in%names(qdata))){
			warning("unique.rowname",unique.rowname,"not found in qdata, row index numbers will be used instead")
			unique.rowname<-FALSE
		}
		if(unique.rowname!=FALSE){
			rownames(qdata)<-qdata[,unique.rowname]
		}
	}
}
#############################################################################################
################################## Ask for rastLUT ##########################################
#############################################################################################

if(is.null(rastLUT) && (is.null(predList) || map==TRUE )){

	if(!is.null(rastLUTfn)){
		if(identical(basename(rastLUTfn),rastLUTfn)){rastLUTfn<-paste(folder,"/",rastLUTfn,sep="")}
		rastLUT<-read.table(file=rastLUTfn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
	}else{
		if( .Platform$OS.type=="windows"){
			LUT.available<-select.list(c("YES","NO"), title="rastLUT available?")
			if(LUT.available=="YES"){
				rastLUTfn<-choose.files(caption="Select raster look up table", filters = Filters["csv",], multi = FALSE)
				if(!is.null(rastLUTfn)){	
						rastLUT<-read.table(file=rastLUTfn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)}
			}
		}
	}
}

if(!is.null(rastLUT)){
	if(is.factor(rastLUT[,1])){rastLUT[,1]<-as.character(rastLUT[,1])}
	if(is.factor(rastLUT[,2])){rastLUT[,2]<-as.character(rastLUT[,2])}
	if(is.factor(rastLUT[,3])){rastLUT[,3]<-as.numeric(as.character(rastLUT[,1]))}
}

if(!is.null(predList) && !is.null(rastLUT) && map==TRUE){
	pred.not.in.LUT<-!(predList%in%rastLUT[,2])
	if(any(pred.not.in.LUT)){
		stop("Predictors ",paste(predList[pred.not.in.LUT]," ",sep=""),"from predList are not found in rastLUT")}
}


#############################################################################################
################################ Load Libraries #############################################
#############################################################################################

## Loads necessary libraries.

if(map==TRUE || (is.null(rastLUT) && is.null(model.obj) && is.null(predList))){library(rgdal)}
if(model.type=="RF" || na.action=="na.roughfix"){library(randomForest)}
if(model.type=="SGB"){library(gbm)}

# Note:	'rgdal' is only used to make map
#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.



#############################################################################################
############################## Select Rasters/Predictors ####################################
#############################################################################################

## This section gets the list of predictors from the user (either explicitly or through
##	pop-up window user selection. If the predictor list is NULL, the getRasts function 
##	allows user to select the predictor rasters from a pop-up browser and select predictors
## 	from a list of raster bands. A lookup table is also generated of raster names, used
##	throughout the program.

if(is.null(model.obj) && is.null(predList)){
	if(is.null(rastLUT)){

		## Gets raster list from user selection.
		if(is.null(rastnmVector)){
			if(.Platform$OS.type=="windows"){
				rastnmVector = getRasts()
			}else{
				stop("To create a model, you must provide one of the following: predList or rastLUT or rastLUTfn or rastnmVector")
			}
		}
		if (is.null(rastnmVector)){
			stop("The raster list is empty")
		}


		## Loops through list of rasters and generates a lookup table with 3 columns: 
		## 	(1) fullname of raster 
		## 	(2) shortname of raster, including band number (if more than one band)
		## 	(3) band number (this value will be 1 if only one layer or band) 
		rastLUT <- {}
		for (i in 1:length(rastnmVector)) {
			rast <- rastnmVector[i]
			if (!is.na(rast)){
				## Gets spatial information for raster
				sp.rast <- open.SpatialGDAL(rast)
				if(is.null(sp.rast)){
					stop("Invalid raster")}
				## Gets number of bands for raster. If null, then only 1 band exists.
				numbands <- dim(sp.rast@grod)[3]
				if(is.na(numbands)){
					numbands <- 1}
				## Gets base name of raster.
				rastnm <- basename(rast)
				rastnm <- strsplit(rastnm,".img")
	
				## Compiles list of rasters and generates a lookup table of fullname and shortname.
				for (j in 1:numbands) {
					## To generate lookup table of fullname and shortname of rasters.
					if (numbands > 1){
						rastLUT <- rbind(rastLUT, c(rast, paste(rastnm,"b",j, sep=""), j))
					}else{
						rastLUT <- rbind(rastLUT, c(rast, paste(rastnm), 1))
					}

				}
			}
		}
		rastLUT<-data.frame(rastLUT,stringsAsFactors = FALSE)
		rastLUT[,3]<-as.numeric(rastLUT[,3])
	}

	## Presents list of possible predictors from raster lookup table to user for selection.
	predList = select.list(rastLUT[,2], "Select predictors", multiple = TRUE)
	if(length(predList)==0 || is.null(predList)){
		stop("predictors must be selected")}

	
	## asks if any predictors are factors
	any.factors<-select.list(c("NO","YES"), title="any predictors catagorical?")
	if(any.factors=="YES"){
		predFactor<-select.list(predList, "Select catagorical predictors", multiple = TRUE)
		if(length(predFactor)==0 || is.null(predFactor)){
			predFactor<-FALSE
		}
	}	
}

if(is.null(model.obj) || predict==TRUE){

	print("About to check if predictors are in dataset")
	print("predList:")
	print(predList)
	print("names(qdata):")
	print(names(qdata))

	pred.not.in.qdata<-!(predList%in%names(qdata))
	if(any(pred.not.in.qdata)){
		stop("Predictors ",paste(predList[pred.not.in.qdata]," ",sep=""),"are not found as column names in the training data and not reconciled by  Look-Up-Table")}
}

#############################################################################################
############################## Select Factored Predictors ###################################
#############################################################################################

print("About to pick factored predictors")
print("predFactor:")
print(predFactor)

if(is.null(model.obj)){
	
	factored.qdata<-sapply(qdata[,match(predList,names(qdata))],is.factor)
	character.data<-sapply(qdata[,match(predList,names(qdata))],is.character)
	factored.qdata<-factored.qdata|character.data

	if(any(predFactor==FALSE)){
		if(any(factored.qdata)){
			stop(	"training or test dataset predictors: ",
				paste(names(factored.qdata)[factored.qdata],collapse=" "),
				" from 'predList' are catagorical predictors (i.e. are non-numeric, such as factors or characters), but are not included in 'predFactor'. Either add these predictors to 'predFactor' or correct the dataset."
				)
		}
	}

	if(!any(predFactor==FALSE)){

		if(any(!names(factored.qdata)[factored.qdata]%in%predFactor)){
			stop(	"training or test dataset predictors: ",
				paste(names(factored.qdata)[factored.qdata][!names(factored.qdata)[factored.qdata]%in%predFactor],collapse=" "),
				" from 'predList' are catagorical predictors (i.e. are non-numeric, such as factors or characters), but are not included in 'predFactor'. Either add these predictors to 'predFactor' or correct the dataset."
				)
		}

		for(i in 1:length(predFactor)){
			qdata[,predFactor[i]]<-factor(qdata[,predFactor[i]])
		}
	}
}

if(!is.null(model.obj)){
		if(is.null(model.obj$levels)){
		if(model.type=="SGB"){
			var.factors<-model.obj$var.type!=0
			if( any(var.factors)){
				model.levels<-as.list(1:sum(var.factors))
				names(model.levels)<-model.obj$var.names[var.factors]
				for(p in 1:sum(var.factors)){
					model.levels[[p]]<-model.obj$var.levels[var.factors][[p]]
				}
				model.obj$levels<-model.levels
			}
		}
		if(model.type=="RF"){
			var.factors<-!sapply(model.obj$forest$xlevels,identical,0)
			if( any(var.factors)){
				model.levels<-as.list(1:sum(var.factors))
				names(model.levels)<-names(var.factors)[var.factors]
				for(p in 1:sum(var.factors)){
					model.levels[[p]]<-model.obj$forest$xlevels[var.factors][[p]]
				}
				model.obj$levels<-model.levels
			}
		}
	}
}



#############################################################################################
############################# Generate Output File Names ####################################
#############################################################################################


print("folder:")
print(folder)

if(is.null(MODELfn)){
	MODELfn<- paste(model.type,"_",response.type,"_",response.name,sep="")}

if(identical(basename(MODELfn),MODELfn)){MODELfn<-paste(folder,"/",MODELfn,sep="")}


if(predict==TRUE){
	if(is.null(MODELpredfn)){
		MODELpredfn<-paste(MODELfn,"_pred.csv",sep="")}
	if(identical(basename(MODELpredfn),MODELpredfn)){MODELpredfn<-paste(folder,"/",MODELpredfn,sep="")}

	if(is.null(DIAGNOSTICfn)){
		DIAGNOSTICfn<-MODELfn}
	if(identical(basename(DIAGNOSTICfn),DIAGNOSTICfn)){DIAGNOSTICfn<-paste(folder,"/",DIAGNOSTICfn,sep="")}
}

if(predict==FALSE){
	if (diagnostics==TRUE){
		warning("predict must equal TRUE to create diagnostic plots")
		diagnostics=FALSE
	}
}

if(map==TRUE){
	if(is.null(asciifn)){
		asciifn<-paste(MODELfn,"_map.txt",sep="")}
	if(identical(basename(asciifn),asciifn)){asciifn<-paste(folder,"/",asciifn,sep="")}

	if(is.null(asciifn.mean)){
		asciifn.mean<-paste(MODELfn,"_mean.txt",sep="")}
	if(identical(basename(asciifn.mean),asciifn.mean)){asciifn.mean<-paste(folder,"/",asciifn,sep="")}

	if(is.null(asciifn.stdev)){
		asciifn.stdev<-paste(MODELfn,"_stdev.txt",sep="")}
	if(identical(basename(asciifn.stdev),asciifn.stdev)){asciifn.stdev<-paste(folder,"/",asciifn.stdev,sep="")}

	if(is.null(asciifn.coefv)){
		asciifn.coefv<-paste(MODELfn,"_coefv.txt",sep="")}
	if(identical(basename(asciifn.coefv),asciifn.coefv)){asciifn.coefv<-paste(folder,"/",asciifn.coefv,sep="")}
}


#############################################################################################
####################################### Build Model #########################################
#############################################################################################

if(!is.null(seed)){
	set.seed(seed)}

if(is.null(model.obj)){
	if (model.type=="RF"){
		model.obj<-create.model(qdata=qdata,
						train=train,
						model.type=model.type,
						folder=FALSE,
						MODELfn=MODELfn,
						predList=predList,
						response.name=response.name,
						response.type=response.type,
						unique.rowname=unique.rowname,
						mtry=mtry,
						seed=NULL,

					# RF arguments:
						ntree=ntree,
					)
	#dev.off()

	}
	if(model.type=="SGB"){
		model.obj<-create.model(qdata=qdata,
						train=train,
						model.type=model.type,		
						folder=FALSE,		
						MODELfn=MODELfn,
						predList=predList,
						response.name=response.name,
						response.type=response.type,			
						unique.rowname=unique.rowname,	
						seed=NULL,
	
					# SGB arguments:
						n.trees=n.trees,                 	
						shrinkage=shrinkage,   	      
                  			interaction.depth=interaction.depth,	
						bag.fraction=bag.fraction,          	
						train.fraction=train.fraction,       	
                  			n.minobsinnode=n.minobsinnode)
	}	
}

#############################################################################################
############################ Make validation predictions ####################################
#############################################################################################

if(predict == TRUE){


	print("Starting validation predictions")
	print(attr(model.obj, "class"))
	print("MODELpredfn:")
	print(MODELpredfn)

	Bpredict.model(	model.obj=model.obj,
				qdata=qdata,
				train=train,
				folder=FALSE,		# No ending slash, to output to working dir = getwd()
				response.name=response.name,
				unique.rowname=unique.rowname,	# Row identifier
				seed=NULL,

			# Model Evaluation Arguments
				diagnostics=diagnostics,		# Diagnostic graphs of predictions (above). Predict must be TRUE
				MODELpredfn=MODELpredfn,
				na.action=na.action,			# also used for mapping
				v.fold=v.fold,

			# SGB arguments
				n.trees=n.trees
				)
	if(diagnostics==TRUE){

		print("starting diagnostics")

		diagnostics.function(	model.obj=model.obj,
					MODELfn=MODELfn,
					MODELpredfn=MODELpredfn,
					response.name=response.name,
					folder=FALSE,
					device.type=device.type,
					DIAGNOSTICfn=DIAGNOSTICfn,
					jpeg.res=jpeg.res,
					device.width=device.width,
					device.height=device.height,
					cex=cex,
					req.sens=req.sens,
					req.spec=req.spec,
					FPC=FPC,
					FNC=FNC)
	#if(.Platform$OS.type=="windows"){bringToTop(which = dev.cur(), stay = FALSE)}
	}


}

#############################################################################################
######################################### Make Map ##########################################
#############################################################################################

## Begin prediction

print("starting production prediction")

if (map==TRUE){
	production.prediction(	model.obj=model.obj,
					rastLUT=rastLUT,
					na.action=na.action,
					folder=FALSE,
					response.name=response.name,
					numrows=numrows,	
					map.sd=map.sd,
					asciifn=asciifn,
					asciifn.mean=asciifn.mean,
					asciifn.stdev=asciifn.stdev,
					asciifn.coefv=asciifn.coefv)

}


#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################


A<-formals(model.map)
A<-mget(names(A),ifnotfound="NULL",envir=as.environment(-1))

ARGfn<-paste(MODELfn,"_arguments.txt",sep="")

if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	A$qdata.trainfn<-"preloaded dataframe"
}

if(is.matrix(qdata.testfn)==TRUE || is.data.frame(qdata.testfn)==TRUE){
	A$qdata.testfn<-"preloaded dataframe"
}

A$datestamp<-Sys.time()
A<-A[c(length(A),1:(length(A)-1))]


print("ARGfn:")
print(ARGfn)

capture.output(print(A),file=ARGfn)




return(model.obj)
}


#############################################################################################
#############################################################################################
############# A Function by any other name would smaell as sweet ############################
#############################################################################################
#############################################################################################

model.run<-model.map