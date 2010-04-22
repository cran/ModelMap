
#############################################################################################
#############################################################################################
################################# Model build Function ######################################
#############################################################################################
#############################################################################################


model.build<-function(	model.type=NULL,	# "RF", "SGB"
				qdata.trainfn=NULL,
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				MODELfn=NULL,
				predList=NULL,
				predFactor=FALSE,
				response.name=NULL,
				response.type=NULL,		# "binary", "continuous",
				seed=NULL,
				na.action="na.omit",
			# RF arguments:
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
			# SGB arguments:
				n.trees=NULL,                 	# number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				train.fraction = 1.0,       	# fraction of data for training,
                  	n.minobsinnode = 10        	# minimum total weight needed in each node
){

## Note: Must have R version 2.8 or greater
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
################################# Select Model Type #########################################
#############################################################################################


if(is.null(model.type)){
	model.type <- select.list(c("RF","SGB"), title="Select model type.")}
if(model.type=="" || is.null(model.type)){
	stop("model.type is required")}

if (model.type == "SGB") {
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package. OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap.")
}

if(!model.type%in%c("RF","SGB")){
	stop("ModelMap currently supports only RF and SGB for model.type")}

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
############################# Generate Output File Names ####################################
#############################################################################################


print("folder:")
print(folder)

## MODELfn
if(is.null(MODELfn)){
	MODELfn<- paste(model.type,"_",response.type,"_",response.name,sep="")}

if(identical(basename(MODELfn),MODELfn)){MODELfn<-paste(folder,"/",MODELfn,sep="")}
	
#############################################################################################
###################################### Load Data ############################################
#############################################################################################


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
	qdata<-qdata.trainfn
}else{
	qdata<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)}


#############################################################################################
#################################### Pick Response ##########################################
#############################################################################################

## If the response variable is NULL, then the user selects variable from pop-up list.
if (is.null(response.name)){
	response.name <- select.list(names(qdata), title="Select response name.")
	if(response.name=="" || is.null(response.name)){
		stop("response.name is required")}	
}

print(paste("response.name =",response.name))

## If model.obj null, ask for response.type

if(is.null(response.type)){
	response.type <- select.list(c("continuous","binary"), title="Select response type.")}
if(response.type=="" || is.null(response.type)){
	stop("response.type is required")}
if(response.type=="categorical"){response.type<-"binary"}	

if(!response.type%in%c("continuous","binary")){
	stop("ModelMap currently supports only continuous and binary for response.type")}

#############################################################################################
##################################### Check Strata ##########################################
#############################################################################################

if(!is.null(strata)){
	if(length(strata)==1){
		if(strata%in%names(qdata)){
			strata<-qdata$strata
		}else{
			stop("'strata' must be either a collumn name from 'qdata' or vector or factor with one element for each row of qdata")
		}
	}
}  


#############################################################################################
################################ Load Libraries #############################################
#############################################################################################

## Loads necessary libraries.

print("loading libraries")

if(model.type=="RF" || na.action=="na.roughfix"){library(randomForest)}
if(model.type=="SGB"){library(gbm)}

#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.



#############################################################################################
################################### Select Predictors #######################################
#############################################################################################

## This section gets the list of predictors from the user (either explicitly or through
##	pop-up window user selection. If the predictor list is NULL, allows user to  
##	select the predictors from collumn names of training data

print("Select predictors")

if(is.null(predList)){

	## Presents list of possible predictors from raster lookup table to user for selection.
	predList = select.list(names(qdata), "Select predictors", multiple = TRUE)
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


#############################################################################################
############################## Select Factored Predictors ###################################
#############################################################################################

print("About to pick factored predictors")
print("predFactor:")
print(predFactor)

factored.qdata<-sapply(qdata[,match(predList,names(qdata))],is.factor)
character.qdata<-sapply(qdata[,match(predList,names(qdata))],is.character)
factored.qdata<-factored.qdata|character.qdata

if(any(predFactor==FALSE)){
	if(any(factored.qdata)){
		stop(	"predictors: ",
			paste(names(factored.qdata)[factored.qdata],collapse=" "),
			" are catagorical predictors (i.e. are non-numeric, such as factors or characters), but are not included in 'predFactor'. Either add these predictors to 'predFactor' or correct the dataset."
			)
	}
}

if(!any(predFactor==FALSE)){

	if(any(!names(factored.qdata)[factored.qdata]%in%predFactor)){
		stop(	"predictors: ",
			paste(names(factored.qdata)[factored.qdata][!names(factored.qdata)[factored.qdata]%in%predFactor],collapse=" "),
			" are catagorical predictors (i.e. are non-numeric, such as factors or characters), but are not included in 'predFactor'. Either add these predictors to 'predFactor' or correct the dataset."
			)
	}

	for(i in 1:length(predFactor)){
		qdata[,predFactor[i]]<-factor(qdata[,predFactor[i]])
	}
}

#############################################################################################
###################################### Deal with NA #########################################
#############################################################################################


if(na.action=="na.roughfix"){
	warning("Replacing NA predictors with median value or most common category")
	qdata[,predList]<-na.roughfix(qdata[,predList])
}



#############################################################################################
####################################### Build Model #########################################
#############################################################################################

if(!is.null(seed)){
	set.seed(seed)}

print("About to create model")

if (model.type=="RF"){
	model.obj<-create.model(qdata=qdata,
					model.type=model.type,
					folder=FALSE,
					predList=predList,
					response.name=response.name,
					response.type=response.type,
					seed=NULL,

				# RF arguments:
					ntree=ntree,
					mtry=mtry,
					replace=replace,
					strata=strata,
					sampsize=sampsize
				)


}
if(model.type=="SGB"){
	model.obj<-create.model(qdata=qdata,
					model.type=model.type,		
					folder=FALSE,		
					predList=predList,
					response.name=response.name,
					response.type=response.type,				
					seed=NULL,
	
				# SGB arguments:
					n.trees=n.trees,                 	
					shrinkage=shrinkage,   	      
                  		interaction.depth=interaction.depth,	
					bag.fraction=bag.fraction,          	
					train.fraction=train.fraction,       	
                  		n.minobsinnode=n.minobsinnode)
}	


#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################


ARGfn<-paste(MODELfn,"_model_building_arguments.txt",sep="")

A<-formals(model.build)
A<-mget(names(A),ifnotfound="NULL",envir=as.environment(-1))

if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	A$qdata.trainfn<-"preloaded dataframe"
}


A$datestamp<-Sys.time()
A<-A[c(length(A),1:(length(A)-1))]

print("ARGfn:")
print(ARGfn)

capture.output(print(A),file=ARGfn)

#############################################################################################
#################################### Return Model Object ####################################
#############################################################################################

return(model.obj)
}


