
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
				response.type=NULL,	# "binary", "continuous",
				unique.rowname=NULL,	# Row identifier
				seed=NULL,
				na.action=NULL,
				keep.data = TRUE,
			# RF arguments:
				ntree=500,
				mtry=NULL,
				replace=TRUE,
				strata=NULL,
				sampsize = NULL,
				proximity = TRUE,
			# SGB arguments:
				n.trees=NULL,                 	# number of trees
				shrinkage=0.001,   	      # shrinkage or learning rate,
                  	interaction.depth=10,		# 1: additive model, 2: two-way interactions, etc.
				bag.fraction = 0.5,          	# subsampling fraction, 0.5 is probably best
				train.fraction = NULL,       	# fraction of data for training,
				nTrain = NULL,
                  	n.minobsinnode = 10,        	# minimum total weight needed in each node
				var.monotone = NULL

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
	stop("'model.type' needed")}

if(!model.type%in%c("RF","SGB")){
	stop("ModelMap currently supports only RF and SGB for 'model.type'")}


#############################################################################################
################################# Select Response Type ######################################
#############################################################################################

## If model.obj null, ask for response.type

if(is.null(response.type)){
	response.type <- select.list(c("continuous","binary","categorical"), title="Select response type.")}
if(response.type=="" || is.null(response.type)){
	stop("'response.type' needed")}	

if(response.type=="categorical" && model.type=="SGB"){
	stop("categorical response only supported for Random Forest models")}

if(!response.type%in%c("continuous","binary","categorical")){
	stop("ModelMap currently supports only continuous binary and categorical for 'response.type'")}

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


print(paste("folder =", folder))

## MODELfn
if(is.null(MODELfn)){
	MODELfn<- paste(model.type,"_",response.type,"_",response.name,sep="")}

if(identical(basename(MODELfn),MODELfn)){
	MODELfn<-file.path(folder,MODELfn)}
	
#############################################################################################
###################################### Load Data ############################################
#############################################################################################


## If training data is NULL, then the user selects file from pop-up browser.
if (is.null(qdata.trainfn)){
	if(.Platform$OS.type=="windows"){
		qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
		if(is.null(qdata.trainfn)){stop("")}
	}else{stop("to create a model or make validation predictions you must provide 'qdata.trainfn'")}
}

## Check if file name is full path or basename
if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
	if(identical(basename(qdata.trainfn),qdata.trainfn)){
		qdata.trainfn<-file.path(folder,qdata.trainfn)}
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
		stop("'response.name' is needed")}	
}

print(paste("response.name =",response.name))

#############################################################################################
##################################### Check Strata ##########################################
#############################################################################################

if(!is.null(strata)){
	if(length(strata)==1){
		if(strata%in%names(qdata)){
			#print(paste("strata:",strata))
			strata<-qdata[,strata]
			#print(paste("strata:",strata))
		}else{
			stop("'strata' must be either a collumn name from 'qdata' or vector or factor with one element for each row of qdata")
		}
	}
}  



#############################################################################################
################################### Select Predictors #######################################
#############################################################################################

## This section gets the list of predictors from the user (either explicitly or through
##	pop-up window user selection. If the predictor list is NULL, allows user to  
##	select the predictors from collumn names of training data

#print("Select predictors")

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

#print("Select factored predictors")
#print(paste("     predFactor =",predFactor)


factored.qdata<-sapply(qdata[,match(predList,names(qdata))],is.factor)
character.qdata<-sapply(qdata[,match(predList,names(qdata))],is.character)
factored.qdata<-factored.qdata|character.qdata

if(any(predFactor==FALSE)){
	if(any(factored.qdata)){
		fact.q<-paste(names(factored.qdata)[factored.qdata],collapse=" ")
		stop(	"predictors ",fact.q,
			" are catagorical predictors i.e. are non numeric such as factors or characters but are not included in 'predFactor' either add these predictors to 'predFactor' or correct the dataset"
			)
	}
}

if(!any(predFactor==FALSE)){

	if(any(!names(factored.qdata)[factored.qdata]%in%predFactor)){
		fact.q<-paste(names(factored.qdata)[factored.qdata][!names(factored.qdata)[factored.qdata]%in%predFactor],collapse=" ")
		stop(	"predictors ",fact.q,
			" are catagorical predictors i.e. are non-numeric such as factors or characters but are not included in 'predFactor' either add these predictors to 'predFactor' or correct the dataset"
			)
	}

	for(i in 1:length(predFactor)){
		qdata[,predFactor[i]]<-factor(qdata[,predFactor[i]])
	}
}

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
######################### Drop unused columns of qdata ######################################
#############################################################################################

qdata<-qdata[,c(predList,response.name)]


#############################################################################################
############################### Determine NA.ACTION #########################################
#############################################################################################

###create logical vector of datarows containing NA in predictors or reponses###

NA.pred<-apply(qdata[,predList],1,function(x){any(is.na(x))})
NA.resp<-is.na(qdata[,response.name])
NA.data<-NA.pred|NA.resp

### if any NA in data ###

#NA.ACTION is character string. Choices: "omit", "rough.fix".
#na.action is the actual function, no quotes.

if(any(NA.data)){

	###Check if na.action is valid###

	NAvalid<-c("na.omit","na.roughfix")
	NAwarn<-"ModelMap currently supports only \"na.omit\" and \"na.roughfix\" for 'na.action'"

	if(is.null(na.action)){
		na.action <- select.list(c("na.omit","na.roughfix"), title="Select na.action")
		if(na.action=="" || is.null(na.action)){
			stop("NA found in data therefore 'na.action' is needed")}
		if(!na.action%in%NAvalid){stop(NAwarn)}
	}else{
		if(is.function(na.action)){
			 stop("quotes needed when specifying the argument 'na.action'")
		}else{
			if(!na.action%in%NAvalid){stop(NAwarn)}
		}
	}

	#make NA.ACTION characters without "na."
	NA.ACTION<-switch(na.action,na.omit="omit",na.roughfix="roughfix","invalid")
	#turn na.action into function
	na.action<-switch(na.action,na.omit=na.omit,na.roughfix=na.roughfix)

}



#############################################################################################
################################ Deal with NA's #############################################
#############################################################################################

if(any(NA.data)){

	###na.omit###

	if(NA.ACTION=="omit"){
		print("Omiting data points with NA predictors or NA response")
		warning(paste(sum(NA.data), "data point(s) with NA values for prdictor or response omitted"))
		qdata<-na.action(qdata)
	}
	

	###na.roughfix###

	if(NA.ACTION=="roughfix"){
		
		#for(i in 1:ncol(qdata)){
		#	print(names(qdata)[i])
		#	print(mode(qdata[,i]))
		#}



		print("Replacing NA predictors and responses with median value or most common category")
		if(any(NA.resp)){warning(paste(sum(NA.data), "data point(s) with NA values (including", sum(NA.resp),"point(s) with NA response) replaced with median/most common response"))
		}else{warning(paste(sum(NA.data), "data point(s) with NA values replaced with median/most common value"))}

		qdata<-na.roughfix(qdata)

		na.ac<-(1:nrow(qdata))[NA.data]
		names(na.ac)<-rownames(qdata)[NA.data]
		class(na.ac)<-NA.ACTION
		attr(qdata,"na.action")<-na.ac
	}
}

#############################################################################################
##################### SGB models: check nTrain vs train.fraction ############################
#############################################################################################

if(model.type=="SGB"){
    if (!is.null(nTrain) && !is.null(train.fraction)) {
        stop("parameters 'nTrain' and 'train.fraction' cannot both be specified")
    }
    else if (!is.null(train.fraction)) {
        warning("parameter 'train.fraction' of gbm.fit is deprecated please specify 'nTrain' instead")
        nTrain <- floor(train.fraction * nrow(qdata))
    }
    else if (is.null(nTrain)) {
        nTrain <- nrow(qdata)
    }
}

print("nrow(qdata):")
print(nrow(qdata))
print("nTrain:")
print(nTrain)

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
					sampsize=sampsize,
					proximity=proximity
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
					keep.data=keep.data,
	
				# SGB arguments:
					n.trees=n.trees,            	
					shrinkage=shrinkage,   	      
                  		interaction.depth=interaction.depth,	
					bag.fraction=bag.fraction,
					nTrain=nTrain,          	
					#train.fraction=train.fraction,       	
                  		n.minobsinnode=n.minobsinnode,
					var.monotone = var.monotone)	
}	

#############################################################################################
################# add a copy of the predictor data to the model object ######################
#############################################################################################

#if(keep.data){model.obj$predictor.data<-qdata[,predList]}

if(keep.data){model.obj$predictor.data<-qdata}

#############################################################################################
########################## add 'na.action' to the model object ##############################
#############################################################################################

if(sum(NA.data>0)){
	model.obj$na.action<-attr(qdata,"na.action")
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

#print(paste("ARGfn =",ARGfn))


capture.output(print(A),file=ARGfn)

#############################################################################################
#################################### Return Model Object ####################################
#############################################################################################

return(model.obj)
}


