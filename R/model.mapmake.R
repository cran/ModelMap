


#############################################################################################
#############################################################################################
################################# Model Map Making ##########################################
#############################################################################################
#############################################################################################



model.mapmake<-function(model.obj=NULL,
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				MODELfn=NULL,
				rastLUTfn=NULL,
			# Model Evaluation Arguments
				na.action=NULL,	# also used for mapping
			# Mapping arguments
				numrows = 500,		# Only used if mapping. Number of rows predicted at a time.				
				map.sd=FALSE,		# Only used if mapping. Generates 3 additional maps (mean,sd,coefvar)
				asciifn=NULL,
				asciifn.mean=NULL,
				asciifn.stdev=NULL,
				asciifn.coefv=NULL,
				make.img=TRUE,
			# SGB arguments
				n.trees=NULL
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
	if(length(modelname)!=1){
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
#	warning("ModelMap currently uses OOB estimation to determine optimal number of trees in SGB model when calling gbm.perf in the gbm package. OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance but is not yet supported in ModelMap.")
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
###################### Extract response.name from model.obj #################################
#############################################################################################

if(is.null(MODELfn) && is.null(asciifn)){
	if(!is.null(model.obj$response)){
		response.name<-model.obj$response
	}

	## If the response variable is NULL, then the user selects variable from pop-up list.
	if (is.null(response.name)){
		stop("if model is built outside of ModelMap then either MODELfn or asciifn must be provided in argument list")
	}
}

#############################################################################################
###################### Extract response.type from model.obj #################################
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

if(response.type=="categorical"){
	
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
	if(is.null(asciifn)){
		MODELfn<-paste(model.type,"_",response.type,"_",response.name,sep="")
	}else{
		MODELfn<-strsplit(basename(asciifn),".txt")[[1]]
	}
}

if(identical(basename(MODELfn),MODELfn)){MODELfn<-paste(folder,"/",MODELfn,sep="")}

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



#############################################################################################
################################ Load Libraries #############################################
#############################################################################################

## Loads necessary libraries.

library(rgdal)
library(raster)

if(model.type=="RF" ){library(randomForest)}
if(model.type=="SGB"){library(gbm)}
#if(response.type=="categorical"){library(gdata)}

# Note:	'rgdal' is only used to make map
#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.
#	'gdata' is used for mapping levels of categorical reponse variables


#############################################################################################
############################### Determine NA.ACTION #########################################
#############################################################################################


###Check Model Object for na.action###

#NA.ACTION is character string. Choices: "omit", "rough.fix".
#na.action is the actual function, no quotes.
#model.NA.ACTION and model.na.action are similar but extracted from model object,
#   and need to be checked against 'model.diagnostics(na.action)'.

model.NA.ACTION<-NULL
model.na.action<-NULL
NA.ACTION<-NULL


###extract na.action from model.obj and check if valid, and if argument is NULL, default to option from model.obj###
if(!is.null(model.obj$na.action)){
	model.NA.ACTION<-class(model.obj$na.action)
	if(!model.NA.ACTION%in%c("omit","roughfix")){
		warning(paste("Model Object was built with 'na.action'",model.NA.ACTION,"which is not supported by model map"))
		model.NA.ACTION<-NULL
	}else{
		model.na.action<-switch(model.NA.ACTION,omit=na.omit,roughfix=na.roughfix)
		if(is.null(na.action)){             #if model.obj has na.action and function call does not then use action from model.obj
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
			stop("NA found in data, therefore na.action is required")}
		if(!na.action%in%NAvalid){stop(NAwarn)}
	}else{
		if(is.function(na.action)){
			 stop("ModelMap requires the use of quotes when specifying the argument 'na.action'")
		}else{
			if(!na.action%in%NAvalid){stop(NAwarn)}
		}
	}
	print(paste("NA.ACTION =",NA.ACTION))
	print(paste("na.action =",na.action))

	#make NA.ACTION characters without "na."
	NA.ACTION<-switch(na.action,na.omit="omit",na.roughfix="roughfix","invalid")
	#turn na.action into function
	na.action<-switch(na.action,na.omit=na.omit,na.roughfix=na.roughfix)
}

	print(paste("NA.ACTION =",NA.ACTION))
	print(paste("is function na.action =",is.function(na.action)))

if(model.type=="SGB" && NA.ACTION=="roughfix"){library(randomForest)}

#############################################################################################
############################# SGB + CV: check for n.trees ###################################
#############################################################################################

if(model.type=="SGB" && is.null(n.trees)){
	if(!is.null(model.obj$best.iter)){
		n.trees<-model.obj$best.iter
	}else{
		n.trees<-model.obj$n.trees
	}
}


#############################################################################################
################################## Ask for rastLUT ##########################################
#############################################################################################

### If rastLUTfn is NULL, then the user selects file from pop-up browser.

if(is.null(rastLUTfn)){
	if( .Platform$OS.type=="windows"){
		LUT.available<-select.list(c("YES","NO"), title="rastLUT available?")
		if(LUT.available=="YES"){
			rastLUTfn<-choose.files(caption="Select raster look up table", filters = Filters["csv",], multi = FALSE)
		}else{
			build.rastLUT(	predList=predList,
						rastLUTfn=paste(MODELfn,"_rastLUT.csv",sep=""))	
		}
	}else{
		stop("You must provide a raster Look Up Table")
	}	
}

### Check if file name is full path or basename

if(is.matrix(rastLUTfn)!=TRUE && is.data.frame(rastLUTfn)!=TRUE){
	if(identical(basename(rastLUTfn),rastLUTfn))
		{rastLUTfn<-paste(folder,"/",rastLUTfn,sep="")
	}
}

### if rastLUT is filename, read in lookup table

if(is.matrix(rastLUTfn)==TRUE || is.data.frame(rastLUTfn)==TRUE){
	rastLUT<-rastLUTfn
	rastLUT<-data.frame(rastLUT)
}else{
	rastLUT<-read.table(file=rastLUTfn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
}


### Check that collumns of rastLUT are correct format

if(is.factor(rastLUT[,1])){rastLUT[,1]<-as.character(rastLUT[,1])}
if(is.factor(rastLUT[,2])){rastLUT[,2]<-as.character(rastLUT[,2])}
if(is.factor(rastLUT[,3])){rastLUT[,3]<-as.numeric(as.character(rastLUT[,1]))}


### Check that all predictors in predList are in rastLUT

pred.not.in.LUT<-!(predList%in%rastLUT[,2])
if(any(pred.not.in.LUT)){
	stop("Predictors ",paste(predList[pred.not.in.LUT]," ",sep=""),"from predList are not found in rastLUT")}


#############################################################################################
######################################### Make Map ##########################################
#############################################################################################

## Begin prediction

print("starting production prediction")


production.prediction(	model.obj=model.obj,
				model.type=model.type,
				rastLUT=rastLUT,
				na.action=na.action,
				NA.ACTION=NA.ACTION,
				response.type=response.type,
				numrows=numrows,	
				map.sd=map.sd,
				asciifn=asciifn,
				asciifn.mean=asciifn.mean,
				asciifn.stdev=asciifn.stdev,
				asciifn.coefv=asciifn.coefv,
				make.img=make.img,
				n.trees=n.trees)


#############################################################################################
############## If response/type=="categorical Write a key to levels #########################
#############################################################################################

if(response.type=="categorical"){
	#mapkey<-mapLevels(model.obj$y)

	Ylev<-levels(model.obj$y)

	if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
		mapkey<-data.frame(row=1:length(Ylev), category=Ylev,integercode=1:length(Ylev))
	}else{
		mapkey<-data.frame(row=1:length(Ylev), category=Ylev,integercode=as.numeric(Ylev))
	}

	MapKeyfn<-paste(MODELfn,"_map_key.csv",sep="")
	write.table(mapkey,file=MapKeyfn,sep=",",row.names=FALSE)
}
	
#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################


A<-formals(model.mapmake)
A<-mget(names(A),ifnotfound="NULL",envir=as.environment(-1))

ARGfn<-paste(MODELfn,"_mapmake_arguments.txt",sep="")

if(is.matrix(rastLUTfn)==TRUE || is.data.frame(rastLUTfn)==TRUE){
	A$rastLUTfn<-"preloaded dataframe"
}

A$datestamp<-Sys.time()
A<-A[c(length(A),1:(length(A)-1))]

print("ARGfn:")
print(ARGfn)

capture.output(print(A),file=ARGfn)

#############################################################################################
####################################### Return nothing ######################################
#############################################################################################


}

