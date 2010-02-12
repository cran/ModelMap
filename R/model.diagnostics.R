


#############################################################################################
#############################################################################################
############################ Model Diagnostics Function #####################################
#############################################################################################
#############################################################################################



model.diagnostics<-function(	model.obj=NULL,
				qdata.trainfn=NULL,
				qdata.testfn=NULL,
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				MODELfn=NULL,
				response.name=NULL,
				unique.rowname=NULL,	# Row identifier
				seed=NULL,
			# Model Evaluation Arguments
				prediction.type=NULL,
				MODELpredfn=NULL,
				na.action="na.omit",	# also used for mapping
				v.fold=10,
				device.type=NULL,	# options: "windows", "jpeg", "postscript", "win.metafile"
				DIAGNOSTICfn=NULL,
				jpeg.res=72,
				device.width=7,
				device.height=7,
				cex=par()$cex,
				req.sens,req.spec,FPC,FNC
			# SGB arguments
				#n.trees=NULL
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


#####################################################################################
############################# check model.obj #######################################
#####################################################################################

if(is.null(model.obj)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select model object", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("Must provide a model object")}
		}else{stop("Must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("File must contain single model object")}
	assign("model.obj",get(modelname))
}

#####################################################################################
##################### Extract Model Type from model.obj #############################
#####################################################################################


model.type.long<-attr(model.obj,"class")
if(is.null(model.type.long)){model.type.long <- "unknown"}
model.type<-switch(model.type.long,	"randomForest"="RF",
						"gbm"="SGB",
						"unknown")
if(model.type=="unknown"){
	stop("model.obj is of unknown type")}
	
print(paste("model.type =",model.type))

if (model.type == "SGB") {
	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package. OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap.")
}

#############################################################################################
########################### Extract predList from model.obj #################################
#############################################################################################


if(model.type == "RF"){predList<-row.names(model.obj$importance)}

if(model.type == "SGB"){predList<-model.obj$var.names}

#############################################################################################
####################### Extract Factored Predictors from model.obj ##########################
#############################################################################################

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


#############################################################################################
################## Extract response.name and response.type from model.obj ###################
#############################################################################################


if(!is.null(model.obj$response)){
	model.response.name<-model.obj$response
	print(paste("     model.response.name =",model.response.name))
	print(paste("     response.name =",response.name))
	if(is.null(response.name)){
		response.name<-model.response.name
	}else{
		if(response.name!=model.response.name){
			warning(paste("supplied model.obj has response",model.response.name,"not supplied response",response.name))
		}
	}
}


## If the response variable is NULL, then the user selects variable from pop-up list.
if (is.null(response.name)){
	response.name <- select.list(names(qdata), title="Select response name.")
	if(response.name=="" || is.null(response.name)){
		stop("response.name is required")}	
}

print(paste("response.name =",response.name))


## extract response.type from model.obj

if(model.type=="RF"){
	response.type<-switch(model.obj$type,"regression"="continuous","classification"="binary","unknown")}
if(model.type=="SGB"){
	response.type<-switch(model.obj$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")}
if(response.type=="unknown"){stop("supplied model.obj has an unknown response type")}


#############################################################################################
######################### Select and Check Prediction Type ##################################
#############################################################################################

if (is.null(prediction.type)){
	if(model.type=="RF"){
		prediction.type <- select.list(c("TEST","CV","OOB"), title="Select prediction type.")}
	if(model.type=="SGB"){
		prediction.type <- select.list(c("TEST","CV","TRAIN"), title="Select prediction type.")}
}

if(prediction.type=="" || is.null(prediction.type)){
	stop("prediction.type is required")}
if(prediction.type=="TRAIN"){
	warning("Predictions will be made made on the training data, and will yeild unrealistic accuracy estimates with regard to independant data")}

if(model.type=="SGB" && prediction.type=="OOB"){
	stop("'OOB' predictions not available for 'SGB' models")}

if(model.type=="RF" && prediction.type=="TRAIN"){
	stop("'TRAIN' predictions not available for 'RF' models")}


#############################################################################################
################# Extract sampsize and strata from model.obj ################################
#############################################################################################


if (model.type == "RF") {
	if(prediction.type=="CV"){

		if(!is.null(model.obj$call$strata)){
			stop("random forest arguments 'strata' and 'sampsize' not currently supported by MOdelMap for cross validation")
		}
		if(!is.null(model.obj$call$sampsize)){
			stop("random forest arguments 'strata' and 'sampsize' not currently supported by MOdelMap for cross validation")
		}
	}
}




#############################################################################################
############################# SGB + CV: check for n.trees ###################################
#############################################################################################

#if(model.type="SGB" && prediction.type="CV"){
#	if(is.null(n.trees)){
#		n.trees==model.obj$n.trees
#	}
#	if(!is.null(model.obj$iter)){
#		n.trees==NULL}
#}

n.trees<-NULL
if(model.type=="SGB"){
	n.trees<-model.obj$n.trees
	if(!is.null(model.obj$iter)){n.trees<-NULL}
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
############################# Generate Output File Names ####################################
#############################################################################################


print("folder:")
print(folder)

## MODELfn
if(is.null(MODELfn)){
	MODELfn<- paste(model.type,"_",response.type,"_",response.name,sep="")}

if(identical(basename(MODELfn),MODELfn)){MODELfn<-paste(folder,"/",MODELfn,sep="")}

## MODELpredfn
if(is.null(MODELpredfn)){
	MODELpredfn<-paste(MODELfn,"_pred",sep="")}
if(identical(basename(MODELpredfn),MODELpredfn)){MODELpredfn<-paste(folder,"/",MODELpredfn,sep="")}



#############################################################################################
################################ Load Training Data #########################################
#############################################################################################

if (prediction.type %in% c("CV","TRAIN","OOB")){

	## If test data supplied give warning
	if (!is.null(qdata.testfn)){
		warning("if 'prediction.type' is ",prediction.type," then 'qdata.testfn' is ignored")}

	## If training data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.trainfn)){
		if(.Platform$OS.type=="windows"){
			qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
			if(is.null(qdata.trainfn)){stop("")}
		}else{stop("If prediction.type is CV, OOB, or TRAIN, you must provide qdata.trainfn")}
	}

	## Check if file name is full path or basename
	if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
		if(identical(basename(qdata.trainfn),qdata.trainfn)){qdata.trainfn<-paste(folder,"/",qdata.trainfn,sep="")}
	}

	## Read in training data
	if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
		qdata<-qdata.trainfn
		qdata<-data.frame(qdata)
	}else{
		qdata<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)
	}

}

#############################################################################################
################################ Load Test Data #############################################
#############################################################################################

if (prediction.type=="TEST"){

	## If training data supplied give warning
	if (!is.null(qdata.trainfn)){
		warning("if 'prediction.type' is 'TEST' then 'qdata.trainfn' is ignored")}

	## If training data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.testfn)){
		if(.Platform$OS.type=="windows"){
			qdata.testfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
			if(is.null(qdata.testfn)){stop("")}
		}else{stop("If prediction.type is TEST, you must provide qdata.testfn")}
	}

	## Check if file name is full path or basename
	if(is.matrix(qdata.testfn)!=TRUE && is.data.frame(qdata.testfn)!=TRUE){
		if(identical(basename(qdata.testfn),qdata.testfn)){qdata.testfn<-paste(folder,"/",qdata.testfn,sep="")}
	}

	## Read in training data
	if(is.matrix(qdata.testfn)==TRUE || is.data.frame(qdata.testfn)==TRUE){
		qdata<-qdata.testfn
		qdata<-data.frame(qdata)
	}else{
		qdata<-read.table(file=qdata.testfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)}

}

#############################################################################################
###################### Check if response or any predList not in qdata #######################
#############################################################################################

if(!(response.name %in% names(qdata))){
	stop("Model Response: ", response.name, " is missing from the diagnostic dataset")}

missing.predictors<-!(predList %in% names(qdata))
if(any(missing.predictors)){
	stop("Model predictors: ", predList[missing.predictors], " are missing from the diagnostic dataset")}


#############################################################################################
######################## Deal with factored predictors ######################################
#############################################################################################





#############################################################################################
######################## Select unique row identifier #######################################
#############################################################################################

if (is.null(unique.rowname)){
	unique.rowname <- select.list(c(names(qdata),"row_index"), title="Select unique row identifier")	
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

#############################################################################################
############################# Check Device Type #############################################
#############################################################################################


if(is.null(device.type)){
	device.type <- select.list(c("default","jpeg","pdf","postscript","win.metafile"), title="Diagnostic Output?", multiple = TRUE)
	device.type <- c(device.type,"default")
}
if(length(device.type)==0 || is.null(device.type)){
	device.type <- "default"
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
################################ Load Libraries #############################################
#############################################################################################

## Loads necessary libraries.

if(model.type=="RF" || na.action=="na.roughfix"){library(randomForest)}
if(model.type=="SGB"){library(gbm)}

#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.


#############################################################################################
############################ Make validation predictions ####################################
#############################################################################################


print("Starting validation predictions")

if(!is.null(seed)){
	set.seed(seed)}

PRED <- prediction.model(	model.obj=model.obj,
						model.type=model.type,
						qdata=qdata,
						folder=folder,		# No ending slash, to output to working dir = getwd()
						response.name=response.name,
						response.type=response.type,
						unique.rowname=unique.rowname,	# Row identifier

					# Model Evaluation Arguments
						prediction.type=prediction.type,
						MODELpredfn=MODELpredfn,
						na.action=na.action,			# also used for mapping
						v.fold=v.fold,

					# SGB arguments
						n.trees=n.trees
						)


print("starting diagnostics")

diagnostics.function(	model.obj=model.obj,
				model.type=model.type,
				predList=predList,
				PRED=PRED,
				MODELfn=MODELfn,
				MODELpredfn=MODELpredfn,
				response.name=response.name,
				response.type=response.type,
				folder=folder,
				device.type=device.type,
				jpeg.res=jpeg.res,
				device.width=device.width,
				device.height=device.height,
				cex=cex,
				req.sens=req.sens,
				req.spec=req.spec,
				FPC=FPC,
				FNC=FNC)
#if(.Platform$OS.type=="windows"){bringToTop(which = dev.cur(), stay = FALSE)}



#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################

ARGfn<-paste(MODELfn,"_model_diagnostics_arguments.txt",sep="")

A<-formals(model.diagnostics)
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


#############################################################################################
##################################### Returns ###############################################
#############################################################################################


return(PRED)
}

