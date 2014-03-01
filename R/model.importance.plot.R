
########## check model type #################


########################################################################
########################################################################
###########################  importance plot   #########################
########################################################################
########################################################################

model.importance.plot<-function(	model.obj.1=NULL, 
						model.obj.2=NULL, 
						model.name.1="Model 1", 
						model.name.2="Model 2",
						imp.type.1=NULL,
						imp.type.2=NULL,
						type.label=TRUE,
						class.1=NULL,
						class.2=NULL, 
						scale.by="sum",
						sort.by="model.obj.1", 
						predList=NULL,
						folder=NULL,
						PLOTfn=NULL,
						device.type=NULL,	# options: "default", "jpeg", "none","postscript"
						jpeg.res=72,
						device.width=7,
						device.height=7,
						cex=par()$cex,
						...){


#############################################################################################
################################### Add to Filters Table ####################################
#############################################################################################

## Adds to file filters to Cran R Filters table.
if(.Platform$OS.type=="windows"){
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))}

###################################################################################
########################## Check Device Type ######################################
###################################################################################

if(is.null(device.type)){
	device.type <- select.list(c("default","jpeg","none","pdf","postscript"), title="Diagnostic Output?", multiple = TRUE)
	device.type <- c(device.type,"default")
}
if(length(device.type)==0 || is.null(device.type)){
	device.type <- "default"
}

device.type[device.type=="windows"]<-"default"
if(any(!device.type%in%c("default","jpeg","none","pdf","postscript","win.metafile"))){
	stop("illegal 'device.type' device types must be one or more of 'default' 'jpeg' 'pdf' or 'postscript'")
}

device.type<-sort(device.type)
if("default"%in%device.type){
	device.type<-c(device.type[device.type!="default"],"default")
}

if("none"%in%device.type){
	device.type<-"none"
}

#if(.Platform$OS.type!="windows"){
#	if("win.metafile" %in% device.type){
#		stop("'win.metafile' is only a legal 'device.type' in a windows environment")}}

###################################################################################
######################## Select Output Folder #####################################
###################################################################################

if(is.null(folder)){
	if(any(device.type%in%c("jpeg","pdf","postscript"))){
		if(.Platform$OS.type=="windows"){
			folder<-choose.dir(default=getwd(), caption="Select directory")
		}else{
			folder<-getwd()}
	}
}

###################################################################################
######################## check output filename ####################################
###################################################################################

if(is.null(PLOTfn)){PLOTfn<- paste(model.name.1,"_",model.name.2,sep="")}
if(identical(basename(PLOTfn),PLOTfn)){
	PLOTfn<-file.path(folder,PLOTfn)}

#####################################################################################
############################# check model.obj #######################################
#####################################################################################

if(is.null(model.obj.1)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select first model", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("file must contain single model object")}
	assign("model.obj.1",get(modelname))
}

if(is.null(model.obj.2)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select second model", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("file must contain single model object")}
	assign("model.obj.2",get(modelname))
}
#####################################################################################
####################### Check Model Type from model.obj #############################
#####################################################################################

model.type.1<-check.model.type(model.obj.1)
model.type.2<-check.model.type(model.obj.2)


#############################################################################################
######################## Extract  response.type from model.obj ##############################
#############################################################################################

#Note: response type is only used to print definition of importance used for the two models

if(model.type.1=="RF"){
	response.type.1<-switch(model.obj.1$type,"regression"="continuous","classification"="classification","unknown")
	if(response.type.1=="classification"){
		if(identical(levels(model.obj.1$y),c("0","1"))){
			response.type.1<-"binary"
		}else{
			response.type.1<-"categorical"}}}
if(model.type.1=="SGB"){
	response.type.1<-switch(model.obj.1$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")}
if(response.type.1=="unknown"){stop("supplied model.obj.1 has an unknown response type")}


if(model.type.2=="RF"){
	response.type.2<-switch(model.obj.2$type,"regression"="continuous","classification"="classification","unknown")
	if(response.type.2=="classification"){
		if(identical(levels(model.obj.2$y),c("0","1"))){
			response.type.2<-"binary"
		}else{
			response.type.2<-"categorical"}}}
if(model.type.2=="SGB"){
	response.type.2<-switch(model.obj.2$distribution$name,"gaussian"="continuous","bernoulli"="binary","unknown")}
if(response.type.2=="unknown"){stop("supplied model.obj.2 has an unknown response type")}

#print(paste("response.type.1 =",response.type.1))
#print(paste("response.type.2 =",response.type.2))

#####################################################################################
###################### Extract Importace from model.obj #############################
#####################################################################################

if(model.type.1=="RF"){
	if(is.null(imp.type.1)){imp.type.1<-1}
	if(!is.null(class.1) && imp.type.1==2){
		warning("no class specific measure for 'imp.type.1=2' therefore importance type changed to 'imp.type.1=1'")
		imp.type.1<-1}
	IMP.1<-imp.extract.rf(model.obj.1,imp.type=imp.type.1,class=class.1)
}else{
	if(model.type.1=="SGB"){
		if(is.null(imp.type.1)){imp.type.1<-2}
		if(!is.null(class.1)){
			warning("no class specific measure for SGB models therefore 'class.1' ignored")}
		print(paste("imp.type.1:",imp.type.1))
		IMP.1<-imp.extract.sgb(model.obj.1,imp.type=imp.type.1)
	}else{stop("model.obj is of unknown type")}
}

if(model.type.2=="RF"){
	if(is.null(imp.type.2)){imp.type.2<-1}
	if(!is.null(class.2) && imp.type.2==2){
		warning("no class specific measure for 'imp.type.2=2' therefore importance type changed to 'imp.type.2=1'")
		imp.type.2<-1}
	IMP.2<-imp.extract.rf(model.obj.2,imp.type=imp.type.2,class=class.2)
}else{
	if(model.type.2=="SGB"){
		if(is.null(imp.type.2)){imp.type.2<-2}
		if(!is.null(class.2)){
			warning("no class specific measure for SGB models therefore 'class.2' ignored")}
		print(paste("imp.type.2:",imp.type.2))
		IMP.2<-imp.extract.sgb(model.obj.2,imp.type=imp.type.2)
	}else{stop("model.obj is of unknown type")}
}

############################################################################################
############################ Print Importance Types ########################################
############################################################################################

if(is.null(class.1)){
	CLASS.1<-"Overall"
}else{
	CLASS.1<-paste("Class",class.1)
}

if(is.null(class.2)){
	CLASS.2<-"Overall"
}else{
	CLASS.2<-paste("Class",class.2)
}


if(model.type.1=="SGB"){
	if(imp.type.1==1){
		if(response.type.1%in%c("continuous")){
			print("model.obj.1 is a continuous SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.1<-"DecPredPerf"} 
		if(response.type.1%in%c("binary")){
			print("model.obj.1 is a binary SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.1<-"DecPredPerf"}
	}
	if(imp.type.1==2){
		if(response.type.1%in%c("continuous")){
			print("model.obj.1 is a continuous SGB model with importance measured by decrease of squared error")
			IMP.MEASURE.1<-"DecSqError"} 
		if(response.type.1%in%c("binary")){
			print("model.obj.1 is a binary SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.1<-"DecSumSqError"}
	}
}
if(model.type.1=="RF"){
	if(response.type.1%in%c("continuous")){
		if(imp.type.1==1){
			print(paste("model.obj.1 is a continuous RF model with", CLASS.1, "Importance measured by %IncMSE"))
			IMP.MEASURE.1<-"%IncMSE"}
		if(imp.type.1==2){
			print(paste("model.obj.1 is a continuous RF model with", CLASS.1, "Importance measured by IncNodePurity"))
			IMP.MEASURE.1<-"IncNodePurity"}
	}
	if(response.type.1%in%c("binary","categorical")){
		if(imp.type.1==1){
			print(paste("model.obj.1 is a", response.type.1, "RF model with", CLASS.1, "Importance measured by MeanDecreaseAccuracy"))
			IMP.MEASURE.1<-"MeanDecAccuracy"}
		if(imp.type.1==2){
			print(paste("model.obj.1 is a", response.type.1, "RF model with", CLASS.1, "Importance measured by MeanDecreaseGini"))
			IMP.MEASURE.1<-"MeanDecGini"}
	}
}

if(model.type.2=="SGB"){
	if(imp.type.2==1){
		if(response.type.2%in%c("continuous")){
			print("model.obj.2 is a continuous SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.2<-"DecPredPerf"} 
		if(response.type.2%in%c("binary")){
			print("model.obj.2 is a binary SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.2<-"DecPredPerf"}
	}
	if(imp.type.2==2){
		if(response.type.2%in%c("continuous")){
			print("model.obj.2 is a continuous SGB model with importance measured by decrease of squared error")
			IMP.MEASURE.2<-"DecSqError"} 
		if(response.type.2%in%c("binary")){
			print("model.obj.2 is a binary SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.2<-"DecSumSqError"}
	}
}
if(model.type.2=="RF"){
	if(response.type.2%in%c("continuous")){
		if(imp.type.2==1){
			print(paste("model.obj.2 is a continuous RF model with", CLASS.2, "Importance measured by %IncMSE"))
			IMP.MEASURE.2<-"%IncMSE"}
		if(imp.type.2==2){
			print(paste("model.obj.2 is a continuous RF model with", CLASS.2, "Importance measured by IncNodePurity"))
			IMP.MEASURE.2<-"IncNodePurity"}
	}
	if(response.type.2%in%c("binary","categorical")){
		if(imp.type.2==1){
			print(paste("model.obj.2 is a", response.type.2, "RF model with", CLASS.2, "Importance measured by MeanDecreaseAccuracy"))
			IMP.MEASURE.2<-"MeanDecAccuracy"}
		if(imp.type.2==2){
			print(paste("model.obj.2 is a", response.type.2, "RF model with", CLASS.2, "Importance measured by MeanDecreaseGini"))
			IMP.MEASURE.2<-"MeanDecGini"}
	}
}

#####################################################################################
############################ Scale Importace ########################################
#####################################################################################

if(!scale.by %in% c("max","sum")){
	stop("scale.by must be either max or sum")}


IMP.1<-imp.scale(IMP.1,scale.by=scale.by)
IMP.2<-imp.scale(IMP.2,scale.by=scale.by)

if(scale.by=="sum"){
	MAX.IMP<-max(IMP.1$imp,IMP.2$imp)

	IMP.1$imp<-IMP.1$imp/MAX.IMP
	IMP.2$imp<-IMP.2$imp/MAX.IMP
}

#####################################################################################
######################### Check predList from model.obj #############################
#####################################################################################

if(!identical(sort(as.character(IMP.1$pred)),sort(as.character(IMP.2$pred)))){
	stop("models contain different predictors")
}

if(!is.null(predList)){
	if(!identical(sort(as.character(IMP.1$pred)),sort(predList))){
		stop("predList contains different predictors than models")
	}

	if(!identical(sort(as.character(IMP.2$pred)),sort(predList))){
		stop("predList contains different predictors than models")
	}
}

#####################################################################################
#################################### Sort Imp #######################################
#####################################################################################


if(!sort.by %in% c("predList","model.obj.1","model.obj.2")){
	stop("sort.by must be 'model.obj.1' or 'model.obj.2' or 'predList'")}



if(sort.by=="model.obj.1"){
	sort.by=IMP.1$pred
}else{
	if(sort.by=="model.obj.2"){
		sort.by=IMP.2$pred
	}else{
		if(sort.by=="predList"){
			sort.by=rev(predList)
		}
	}
}

IMP.1<-IMP.1[match(sort.by,IMP.1$pred),]
IMP.2<-IMP.2[match(sort.by,IMP.2$pred),]

#####################################################################################
################################### Make plot #######################################
#####################################################################################

	names.arg=IMP.1$pred
	NCHAR<-max(nchar(as.character(names.arg)))


#if(!"none"%in%device.type){
for(i in 1:length(device.type)){

###################################################################
### Output filenames ###

if(device.type[i] == "jpeg"){
	IMPORTANCEfn<-paste(PLOTfn,".jpg",sep="")
}

if(device.type[i] == "pdf"){
	IMPORTANCEfn<-paste(PLOTfn,".pdf",sep="")
}

if(device.type[i] == "postscript"){
	IMPORTANCEfn<-paste(PLOTfn,".ps",sep="")
}


if(device.type[i]=="default"){dev.new(width = device.width, height = device.height,  record = TRUE)}
if(device.type[i]=="jpeg"){jpeg(filename=IMPORTANCEfn,width = device.width, height = device.height, res=jpeg.res, units="in")}
if(device.type[i]=="postscript"){postscript(file=IMPORTANCEfn,width = device.width, height = device.height)}
if(device.type[i]=="pdf"){pdf(file=IMPORTANCEfn,width = device.width, height = device.height)}


	op<-par(mar=par()$mar+c(0,(2*NCHAR/3)-2,0,0),cex=cex)

	barplot(-IMP.1$imp,horiz=TRUE,xlim=c(-1,1),las=1,axes=F,names.arg=names.arg,...)
	barplot(IMP.2$imp,horiz=TRUE,xlim=c(-1,1),las=1,col="black",add=TRUE,axes=F,...)

	IMP.TEXT<-paste(IMP.MEASURE.1,"  ",IMP.MEASURE.2)

	#mtext(IMP.TEXT,side=1, line=0, adj=.5, cex=cex)

	if(type.label){
		mtext(IMP.MEASURE.1,side=1, line=0, adj=0, cex=1*cex)
		mtext(IMP.MEASURE.2,side=1, line=0, adj=1, cex=1*cex)}

	mtext(model.name.1,side=1, line=1.5, adj=0, cex=1.3*cex)
	mtext(model.name.2,side=1, line=1.5, adj=1, cex=1.3*cex)

	par(op)
if(!device.type[i]%in%c("default","none")){dev.off()}
}
#}

#############################################################################
}


