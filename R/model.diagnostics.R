

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
				diagnostic.flag=NULL,
				seed=NULL,
			# Model Evaluation Arguments
				prediction.type=NULL,
				MODELpredfn=NULL,
				na.action=NULL,	# also used for mapping
				v.fold=10,
				device.type=NULL,	# options: "default", "jpeg", "none","postscript"
				DIAGNOSTICfn=NULL,
				jpeg.res=72,
				device.width=7,
				device.height=7,
				cex=par()$cex,
				req.sens,req.spec,FPC,FNC,
			# SGB arguments
				n.trees=NULL
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

#model.type<-switch(model.type.long,	"randomForest"="RF",
#						"gbm"="SGB",
#						"unknown")

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
	
print(paste("model.type =",model.type))

#if (model.type == "SGB") {
#	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package but OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive however using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap")
#}

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
######################### Extract response.name from model.obj ##############################
#############################################################################################


if(!is.null(model.obj$response)){
	model.response.name<-model.obj$response
	#print(paste("     model.response.name =",model.response.name))
	#print(paste("     response.name =",response.name))
	if(is.null(response.name)){
		response.name<-model.response.name
	}else{
		if(response.name!=model.response.name){
			warning("supplied model.obj has response ",model.response.name," not supplied response ",response.name)
		}
	}
}


## If the response variable is NULL, then the user selects variable from pop-up list.
if (is.null(response.name)){
	response.name <- select.list(names(qdata), title="Select response name.")
	if(response.name=="" || is.null(response.name)){
		stop("response.name is needed")}	
}

print(paste("response.name =",response.name))


#############################################################################################
######################## Extract  response.type from model.obj ##############################
#############################################################################################

if(model.type=="RF"){
	response.type<-switch(model.obj$type,"regression"="continuous","classification"="classification","unknown")
	if(response.type=="classification"){
		if(identical(levels(model.obj$y),c("0","1"))){
			response.type<-"binary"
		}else{
			response.type<-"categorical"}}}
if(model.type=="SGB"){
	response.type<-switch(model.obj$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")}
if(response.type=="unknown"){stop("supplied model.obj has an unknown response type")}

print(paste("response.type=",response.type))



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
	stop("prediction.type is needed")}
if(prediction.type=="TRAIN"){
	warning("predictions will be made made on the training data and will yeild unrealistic accuracy estimates with regard to independant data")}

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
			stop("random forest arguments 'strata' and 'sampsize' not currently supported by ModelMap for cross validation")
		}
		if(!is.null(model.obj$call$sampsize)){
			stop("random forest arguments 'strata' and 'sampsize' not currently supported by ModelMap for cross validation")
		}
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
############################# Generate Output File Names ####################################
#############################################################################################


print(paste("folder =",folder))
print(paste("MODELfn =",MODELfn))
## MODELfn
if(is.null(MODELfn)){
	MODELfn<- paste(model.type,"_",response.type,"_",response.name,sep="")}

if(identical(basename(MODELfn),MODELfn)){
	MODELfn<-file.path(folder,MODELfn)}


print(paste("MODELpredfn =",MODELpredfn))
## MODELpredfn
if(is.null(MODELpredfn)){
	MODELpredfn<-paste(MODELfn,"_pred",sep="")}

if(identical(basename(MODELpredfn),MODELpredfn)){
	MODELpredfn<-file.path(folder,MODELpredfn)}



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
		}else{stop("if prediction.type is CV OOB or TRAIN you must provide 'qdata.trainfn'")}
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
		if(identical(basename(qdata.testfn),qdata.testfn)){
			qdata.testfn<-file.path(folder,qdata.testfn)}
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

datakeep<-c(predList,response.name)

#############################################################################################
######################## Deal with factored response ########################################
#############################################################################################

if(response.type=="categorical"){
	if(!is.factor(qdata[,response.name])){
		qdata[,response.name]<-factor(qdata[,response.name])}

	#levels.qdata<-levels(qdata[,response.name])
	#levels.model<-levels(model.obj$y)

}


#############################################################################################
######################## Select unique row identifier #######################################
#############################################################################################

if (is.null(unique.rowname)){
	unique.rowname <- select.list(c(names(qdata),"row_index"), title="Select unique row identifier")	
	if(unique.rowname!="row_index" && unique.rowname!=""){

		#if(class(RF$na.action)=="omit")
		print(paste("A: rownames length =",length(rownames(qdata)) ))
		print(paste("A: qdata length =",length(qdata[,unique.rowname]) ))

		if(anyDuplicated(qdata[,unique.rowname])){
			stop("'unique.rowname' contains duplicated values")}
		rownames(qdata)<-qdata[,unique.rowname]

	}	
}else{
	if(!(unique.rowname%in%names(qdata))){
		warning("unique.rowname ",unique.rowname," not found in qdata therefore row index numbers will be used instead")
		unique.rowname<-FALSE
		rownames(qdata)<-1:nrow(qdata)
	}
	if(unique.rowname!=FALSE){
		print(paste("B: rownames length =",length(rownames(qdata))))
		print(paste("B: qdata length =",length(qdata[,unique.rowname])))

		if(anyDuplicated(qdata[,unique.rowname])){
			stop("'unique.rowname' contains duplicated values")}
		rownames(qdata)<-qdata[,unique.rowname]
	}
}

#############################################################################################
############################ check 'diagnostic.flag' ########################################
#############################################################################################

print("checking diagnostic.flag")

DFLAG<-FALSE
if(!is.null(diagnostic.flag)){DFLAG<-TRUE}

if(DFLAG){

	if(!(diagnostic.flag%in%names(qdata))){
		if(prediction.type=="TEST"){
			stop("diagnostic.flag ",diagnostic.flag," is not found in 'qdata.testfn'")
		}else{
			stop("diagnostic.flag ",diagnostic.flag," is not found in 'qdata.trainfn'")}
		DFLAG<-FALSE
	}

	if(is.factor(qdata[,diagnostic.flag])){
		qdata[,diagnostic.flag]<-as.logical(qdata[,diagnostic.flag])
	}

	if(is.character(qdata[,diagnostic.flag])){
		qdata[,diagnostic.flag]<-as.logical(as.factor(qdata[,diagnostic.flag]))
	}

	if(is.numeric(qdata[,diagnostic.flag])){
		if(all(qdata[,diagnostic.flag]%in%c(0,1))){
			qdata[,diagnostic.flag]<-as.logical(qdata[,diagnostic.flag])
		}else{
			stop("diagnostic.flag must be a column of either 0 and 1 (0=FALSE, and 1=TRUE), or of TRUE and FALSE")
		}
	} 

	if(any(is.na(qdata[,diagnostic.flag]))){
		stop("unable to interpret all values of diagnostics flag as true or false, NA's generated")
	}

}

if(DFLAG){
	datakeep<-c(datakeep,diagnostic.flag)
}


#############################################################################################
################################## Drop unused columns ######################################
#############################################################################################

print("dropping unused collumns")


print("  qdata dimensions before")
print(dim(qdata))

qdata<-qdata[,datakeep]

print("  qdata dimensions after")
print(dim(qdata))



#############################################################################################
########## Check for factored predictors with levels not found in training data #############
#############################################################################################
	
print("checking factors")

if(!is.null(model.obj$levels)){
	invalid.levels<-model.obj$levels

	for(p in names(model.obj$levels)){
		valid.levels<-c(model.obj$levels[[p]],NA)
		invalid.levels[[p]]<-unique(qdata[,p][!qdata[,p]%in%valid.levels] )
		qdata[,p]<-factor(qdata[,p],levels=model.obj$levels[[p]])
		if(length(invalid.levels[[p]])>0){
			invalid.lev<-paste(invalid.levels[[p]],collapse=", ")
			warning(	"categorical factored predictor ",p," contains levels: ",invalid.lev,
					" not found in training data and these categories will be treated as NA and either omited or replaced")
		}
	}
}

#############################################################################################
######################### Omit rows with NA response ########################################
#############################################################################################

NA.resp<-is.na(qdata[,response.name])

if(any(NA.resp)){
	warning("Omiting ", sum(NA.resp), " datapoints with NA values for response")
	qdata <- qdata[!NA.resp,]
}

#############################################################################################
############################### Determine NA.ACTION #########################################
#############################################################################################

print("Na action")

###create logical vector of datarows containing NA in predictors or reponses###

NA.pred<-apply(qdata[,predList],1,function(x){any(is.na(x))})


#NA.resp<-is.na(qdata[,response.name])
#NA.data<-NA.pred|NA.resp
#print(any(NA.data))


###Check Model Object for na.action###

#NA.ACTION is character string. Choices: "omit", "rough.fix".
#na.action is the actual function, no quotes.
#model.NA.ACTION and model.na.action are similar but extracted from model object,
#   and need to be checked against 'model.diagnostics(na.action)'.

model.NA.ACTION<-NULL
model.na.action<-NULL
NA.ACTION<-NULL

if(any(NA.pred) || (any(var.factors) && prediction.type=="CV")){

	
	###extract na.action from model.obj and check if valid, and if argument is NULL, default to option from model.obj###
	if(!is.null(model.obj$na.action)){
		model.NA.ACTION<-class(model.obj$na.action)
		print(paste("model.NA.ACTION =",model.NA.ACTION))
		if(!model.NA.ACTION%in%c("omit","roughfix")){
			warning("Model Object was built with 'na.action' ",model.NA.ACTION," which is not supported by model map")
			model.NA.ACTION<-NULL
		}else{
			model.na.action<-switch(model.NA.ACTION,omit=na.omit,roughfix=na.roughfix)
			if(is.null(na.action)){			#if model.obj has na.action and function call does not then use action from model.obj
				NA.ACTION<-model.NA.ACTION	#if this step is reached then whole next section not needed
				na.action<-model.na.action
			}
		}
	}



	###Check if na.action from argument is valid###

	NAvalid<-c("na.omit","na.roughfix")
	NAwarn<-"ModelMap currently supports only \"na.omit\" and \"na.roughfix\" for 'na.action'"

	if(is.null(NA.ACTION)){	
		if(is.null(na.action)){
			na.action <- select.list(c("na.omit","na.roughfix"), title="Select na.action")
			if(na.action=="" || is.null(na.action)){
				stop("NA found in data, therefore na.action needed")}
			if(!na.action%in%NAvalid){stop(NAwarn)}
		}else{
			if(is.function(na.action)){
				 stop("quotes needed when specifying the argument 'na.action'")
			}else{
				if(!na.action%in%NAvalid){stop(NAwarn)}
			}
		}

		NA.ACTION<-switch(na.action,na.omit="omit",na.roughfix="roughfix","invalid")
		#turn na.action into function
		na.action<-switch(na.action,na.omit=na.omit,na.roughfix=na.roughfix)
	}

	###Check for impossible combo of model na.action, diagnostics na.action, and prediction.type=OOB
	if(!is.null(model.NA.ACTION)){
		if(model.NA.ACTION=="omit" && NA.ACTION=="roughfix" && prediction.type=="OOB"){
			warning("'model.obj' was created with 'na.action=\"na.omit\"' therefore OOB predictions not available on data rows containing NA and 'na.action' defaulting to \"na.omit\"")
			NA.ACTION<-"omit"
			na.action<-na.omit
		}
		
	}
	
}

				#print("NA.ACTION:")
				#print(NA.ACTION)
				#print("na.action:")
				#print(na.action)
				#print("model.NA.ACTION:")
				#print(model.NA.ACTION)
				#print("model.na.action:")
				#print(model.na.action)
#############################################################################################
################################ Deal with NA's #############################################
#############################################################################################

print("starting dealing with NA")

if(any(NA.pred)){

	###na.omit###

	if(NA.ACTION=="omit"){
		print("Omitting data points with NA predictors")
		warning("Omitting ", sum(NA.pred), " datapoints with NA values for predictors")
		qdata<-na.action(qdata)
	}
	

	###na.roughfix###

	if(NA.ACTION=="roughfix"){
		warning("Rough fixing ", sum(NA.pred), " data points with NA values for predictors by replacing NA with median or most common value")

		qdata<-na.action(qdata)

		na.ac<-(1:nrow(qdata))[NA.pred]
		names(na.ac)<-rownames(qdata)[NA.pred]
		class(na.ac)<-NA.ACTION
		attr(qdata,"na.action")<-na.ac
	}
}

print("done dealing with NA")


##for OOB models, check number of rows is equal to data used in building model
#if(prediction.type=="OOB"){
#	if( nrow(qdata)!= length(model.obj$y)){
#		stop("'prediction.type' is OOB but number of rows in 'qdata.train' does not match dataset used to build model")
#	}
#}


#############################################################################################
############################# Check Device Type #############################################
#############################################################################################

print("check device type")

if(is.null(device.type)){
	device.type <- select.list(c("default","jpeg","none","pdf","postscript"), title="Diagnostic Output?", multiple = TRUE)
	device.type <- c(device.type,"default")
}
if(length(device.type)==0 || is.null(device.type)){
	device.type <- "default"
}

if(!is.null(device.type)){
	device.type[device.type=="windows"]<-"default"
	if(any(!device.type%in%c("default","jpeg","none","pdf","postscript"))){
		stop("illegal 'device.type' device types must be one or more of 'default' 'jpeg' 'pdf' or 'postscript'")
	}
	device.type<-sort(device.type)
	if("default"%in%device.type){
		device.type<-c(device.type[device.type!="default"],"default")
	}
}

if("none"%in%device.type){
	device.type<-"none"
}



#############################################################################################
############################# SGB + CV: check for n.trees ###################################
#############################################################################################

print("check n.trees")

if(model.type=="SGB"){
	if(is.null(n.trees)){
		if(!is.null(model.obj$best.iter)){
			n.trees<-model.obj$best.iter
		}else{
			n.trees<-model.obj$n.trees
		}
	}
	#if(!is.null(model.obj$nTrain)){
	#	nTrain <- model.obj$nTrain
	#}else{
	#	nTrain <- length(model.obj$train.error)
	#}
}


#############################################################################################
############################ Make validation predictions ####################################
#############################################################################################


if(!is.null(seed)){
	set.seed(seed)}


print("starting predictions")

PRED <- prediction.model(		model.obj=model.obj,
						model.type=model.type,
						qdata=qdata,
						folder=folder,		# No ending slash, to output to working dir = getwd()
						response.name=response.name,
						response.type=response.type,
						unique.rowname=unique.rowname,	# Row identifier

					# Model Evaluation Arguments
						prediction.type=prediction.type,
						MODELpredfn=MODELpredfn,
						v.fold=v.fold,

						na.action=na.action,			
						NA.ACTION=NA.ACTION,
						model.na.action=model.na.action,			
						model.NA.ACTION=model.NA.ACTION,

					# SGB arguments
						n.trees=n.trees
						)


###if na.action="na.roughfix" must turn estimated responses back to NA


#############################################################################################
########################### If DFLAG=T remove unused rows ###################################
#############################################################################################

if(DFLAG){
	PREDICTIONfn<-paste(MODELpredfn,".csv",sep="")
	PRED.OUT<-cbind(PRED,qdata[,diagnostic.flag])
	names(PRED.OUT)<-c(names(PRED),diagnostic.flag)
	write.table(PRED.OUT,file=PREDICTIONfn,sep=",",row.names=FALSE)
}





print("removing rows")

#print("PRED dimensions before:")
#print(dim(PRED))

if(DFLAG){
	PRED<-PRED[qdata[,diagnostic.flag],]
}


#print("PRED dimensions after:")
#print(dim(PRED))
#print(PRED)

#print(paste("diagnostic.flag =",diagnostic.flag))

#print("names(qdata)=")
#print(names(qdata))

#print("qdata[,diagnostic.flag]")
#print(qdata[,diagnostic.flag])
#############################################################################################
################################### Run Diagnostics #########################################
#############################################################################################



print("starting diagnostics")

diagnostics.function(	model.obj=model.obj,
				model.type=model.type,
				predList=predList,
				PRED=PRED,
				MODELfn=MODELfn,
				MODELpredfn=MODELpredfn,
				response.name=response.name,
				response.type=response.type,
				prediction.type=prediction.type,
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

print("ending diagnostics")

#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################

print("starting argument record")
ARGfn<-paste(MODELfn,"_model_diagnostics_arguments.txt",sep="")

A<-formals(model.diagnostics)
envir<-as.environment(-1)
A<-mget(names(A),ifnotfound="NULL",envir=envir)

#ARGfn<-paste(MODELfn,"_arguments.txt",sep="")

if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	A$qdata.trainfn<-"preloaded dataframe"
}


if(is.matrix(qdata.testfn)==TRUE || is.data.frame(qdata.testfn)==TRUE){
	A$qdata.testfn<-"preloaded dataframe"
}

A$datestamp<-Sys.time()
A<-A[c(length(A),1:(length(A)-1))]

#print(A)
#capture.output(print(A),file=ARGfn)
print("ending argument record")

#############################################################################################
##################################### Returns ###############################################
#############################################################################################


return(PRED)
}

