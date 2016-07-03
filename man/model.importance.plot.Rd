\name{model.importance.plot}

\alias{model.importance.plot}

\title{Compares the variable importance of two models with a back to back barchart.}

\description{Takes two models and produces a back to back bar chart to compare the importance of the predictor variables. Models can be any combination of Random Forest or Stochastic Gradient Boosting, as long as both models have the same predictor variables. }

\usage{
model.importance.plot(model.obj.1 = NULL, model.obj.2 = NULL, 
model.name.1 = "Model 1", model.name.2 = "Model 2", imp.type.1 = NULL, 
imp.type.2 = NULL, type.label=TRUE, class.1 = NULL, class.2 = NULL, 
quantile.1=NULL, quantile.2=NULL,
col.1="grey", col.2="black", scale.by = "sum", sort.by = "model.obj.1", 
cf.mincriterion.1 = 0, cf.conditional.1 = FALSE, cf.threshold.1 = 0.2, 
cf.nperm.1 = 1, cf.mincriterion.2 = 0, cf.conditional.2 = FALSE, 
cf.threshold.2 = 0.2, cf.nperm.2 = 1, predList = NULL, folder = NULL, 
PLOTfn = NULL, device.type = NULL, res=NULL, jpeg.res = 72, 
device.width = 7,  device.height = 7, units="in", pointsize=12, 
cex=par()$cex,...)
}

\arguments{

  \item{model.obj.1}{\code{R} model object. The model object to use for left side of barchart. The model object must be of type \code{"RF"} (random forest), \code{"QRF"} (quantile random forest), \code{"CF"} (conditional forest), or \code{"SGB"} (stochastic gradient boosting). The model object must have the exact same predictors as \code{model.obj.2}.}

  \item{model.obj.2}{\code{R} model object. The model object to use for right side of barchart. The model object must be of type \code{"RF"} (random forest), \code{"QRF"} (quantile random forest), \code{"CF"} (conditional forest), or \code{"SGB"} (stochastic gradient boosting).The model object must have the exact same predictors as \code{model.obj.1}.}

  \item{model.name.1}{String. Label for left side of barchart.}

  \item{model.name.2}{String. Label for right side of barchart.}
 
  \item{imp.type.1}{Number. Type of importance to use for model 1. Importance type 1 is permutation based, as described in Breiman (2001).  Importance type 2 is model based. For RF models is the decrease in node impurities attributable to each predictor variable. For SGB models, it is the reduction attributable to each variable in predicting the gradient on each iteration. Default for random forest models is \code{imp.type.1 = 1}. Default for SGB models is \code{imp.type.1 = 2}.}

  \item{imp.type.2}{Number. Type of importance to use for model 2. Importance type 1 is permutation based, as described in Breiman (2001).  Importance type 2 is model based. For RF models is the decrease in node impurities attributable to each predictor variable. For SGB models, it is the reduction attributable to each variable in predicting the gradient on each iteration. Default for random forest models is \code{imp.type.2 = 1}. Default for SGB models is \code{imp.type.2 = 2}.}

  \item{type.label}{Logical. Should axis labels include importance type for each side of plot.}

  \item{class.1}{String. For binary and categorical random forest models. If the name a class is specified, the class-specific relative influence is used for plot. If \code{class.1 = NULL} overall relative influence used for plot.}

  \item{class.2}{String. For binary and categorical random forest models. If the name a class is specified, the class-specific relative influence is used for plot. If \code{class.2 = NULL} overall relative influence used for plot.}

  \item{quantile.1}{Numeric. QRF models. Quantile to use for model 1. Must be one of the quantiles used in building the QRF model.}
 
  \item{quantile.2}{Numeric. QRF models. Quantile to use for model 2. Must be one of the quantiles used in building the QRF model.}

  \item{col.1}{String. For binary and categorical random forest models. Color to use for bars for model 1. Defaults to grey.}

  \item{col.2}{String. For binary and categorical random forest models. Color to use for bars for model 2. Defaults to black.}

  \item{scale.by}{String. Scale by: \code{"max"} or \code{"sum"}. When \code{scale.by="max"} the importance are scaled for each model so that the maximum importance for each model fills the graph. When \code{scale.by="sum"}, the importance for each model are scaled to sum to 100.}

  \item{sort.by}{String. Sort by: \code{"model.obj.1"}, \code{"model.obj.2"}, \code{"predList"}. Gives the order to draw the bars for the predictor variables. When \code{sort.by="model.obj.1"} the predictors are sorted largest to smallest based on importance from model 1. When \code{sort.by="model.obj.2"} the predictors are sorted largest to smallest based on importance from model 2. When \code{sort.by="predList"} the predictors are sorted to match the order given in \code{"predList"}.}

 \item{cf.mincriterion.1}{Number. CF models. The value of the test statistic or 1 - p-value that must be exceeded in order to include a split in the computation of the importance. The default \code{cf.mincriterion.1 = 0} guarantees that all splits are included.}

 \item{cf.conditional.1}{Logical. CF models. A logical determining whether unconditional or conditional computation of the importance is performed for \code{model.obj.1}.}

\item{cf.threshold.1}{Number. CF models. The value of the test statistic or 1 - p-value of the association between the variable of interest and a covariate that must be exceeded inorder to include the covariate in the conditioning scheme for the variable of interest (only relevant if \code{conditional = TRUE}).}
					
\item{cf.nperm.1}{Number. CF models. The number of permutations performed.}				 
						
 \item{cf.mincriterion.2}{Number. CF models. The value of the test statistic or 1 - p-value that must be exceeded in order to include a split in the computation of the importance. The default \code{cf.mincriterion.2 = 0} guarantees that all splits are included.}

 \item{cf.conditional.2}{Logical. CF models. A logical determining whether unconditional or conditional computation of the importance is performed for \code{model.obj.2}.}

\item{cf.threshold.2}{Number. CF models. The value of the test statistic or 1 - p-value of the association between the variable of interest and a covariate that must be exceeded inorder to include the covariate in the conditioning scheme for the variable of interest (only relevant if \code{conditional = TRUE}).}
					
\item{cf.nperm.2}{Number. CF models. The number of permutations performed.}					

  \item{predList}{String.  A character vector of the predictor short names used to build the models. If \code{sort.by="predList"}, then \code{predList} is used to specify the order to draw the predictors in the barchart.}

  \item{folder}{ String.  The folder used for all output.  Do not add ending slash to path string.  If \code{folder = NULL} (default), a GUI interface prompts user to browse to a folder.  To use the working directory, specify \code{folder = getwd()}.}

  \item{PLOTfn}{ String.  The file name to use to save the generated graphical plots. If \code{PLOTfn = NULL} a default name is generated by pasting \code{model.name.1_model.name.2}. The filename can be the full path, or it can be the simple basename, in which case the output will be to the folder specified by \code{folder}.}

 \item{device.type}{ String or vector of strings.  Model validation.  One or more device types for graphical output from model validation diagnostics. 

Current choices:

\tabular{lllll}{
	  \tab \tab \tab \code{"default"} \tab default graphics device\cr
	  \tab \tab \tab \code{"jpeg"} \tab *.jpg files\cr
	  \tab \tab \tab \code{"none"} \tab no graphics device generated\cr	
	  \tab \tab \tab \code{"pdf"} \tab *.pdf files\cr
	  \tab \tab \tab \code{"png"} \tab *.png files\cr
	  \tab \tab \tab \code{"postscript"} \tab *.ps files\cr
	  \tab \tab \tab \code{"tiff"} \tab *.tif files }

 }

  \item{res}{ Integer.  Model validation.  Pixels per inch for jpeg, png, and tiff plots.  The default is 72dpi, good for on screen viewing. For printing, suggested setting is 300dpi. }

  \item{jpeg.res}{ Integer.  Model validation.  Deprecated. Ignored unless \code{res} not provided. }

  \item{device.width}{ Integer.  Model validation.  The device width for diagnostic plots in inches. }

  \item{device.height}{ Integer.  Model validation.  The device height for diagnostic plots in inches. }

  \item{units}{ Model validation.  The units in which \code{device.height} and \code{device.width} are given. Can be \code{"px"} (pixels), \code{"in"} (inches, the default), \code{"cm"} or \code{"mm"}. }

  \item{pointsize}{ Integer.  Model validation.  The default pointsize of plotted text, interpreted as big points (1/72 inch) at \code{res} ppi}

  \item{cex}{ Integer.  Model validation.  The cex for diagnostic plots. }

  \item{\dots}{Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).}

}
\details{
The importance measures used in this plot depend on the model type (RF verses SGB) and the response type (continuous, categorical, or binary). 

Importance type 1 is permutation based, as described in Breiman (2001). Importance is calculated by randomly permuting each predictor variable and computing the associated reduction in predictive performance using Out Of Bag error for RF models and training error for SGB models. Note that for SGB models permutation based importance measures are still considered experimental.  Importance type 2 is model based. For RF models, importance type 2 is calculated by the decrease in node impurities attributable to each predictor variable. For SGB models, importance type 2 is the reduction attributable to each variable in predicting the gradient on each iteration as described in described in Friedman (2001). 


For RF models:
\tabular{lllclll}{
	  \tab \tab  \code{response type} \tab \code{type}	\tab			\tab \tab Importance Measure \cr	
	  \tab \tab  \code{"continuous"}  \tab \code{1}   	\tab permutation 	\tab \tab \%IncMSE \cr
	  \tab \tab  \code{"binary"}      \tab \code{1}    	\tab permutation 	\tab \tab Mean Decrease Accuracy \cr
	  \tab \tab  \code{"categorical"} \tab \code{1}    	\tab permutation 	\tab \tab Mean Decrease Accuracy \cr
	  \tab \tab  \code{"continuous"}  \tab \code{2}    	\tab node impurity	\tab \tab Residual sum of squares \cr
	  \tab \tab  \code{"binary"}      \tab \code{2}    	\tab node impurity	\tab \tab Mean Decrease Gini \cr
	  \tab \tab  \code{"categorical"} \tab \code{2}    	\tab node impurity	\tab \tab Mean Decrease Gini   }
  
For Random Forest models, if \code{imp.type} not specified, importance type defaults to \code{imp.type} of \code{1} - permutation importance. For SGB models, permutation importance is considered experimental so importance defaults to \code{imp.type} of \code{2} - reduction of gradient of the loss function.

Also, for binary and categorical Random Forest models, class specific importance plots can be generated by the use of the \code{class} argument. Note that class specific importance is only available for Random Forest models with importance type 1.

For CF models:
\tabular{lllclll}{
	  \tab \tab  \code{response type} \tab \code{type}	\tab			\tab \tab Importance Measure \cr	
	  \tab \tab  \code{"continuous"}  \tab \code{1}   	\tab permutation 	\tab \tab Mean Decrease Accuracy \cr
	  \tab \tab  \code{"binary"}      \tab \code{1}    	\tab permutation 	\tab \tab Mean Decrease Accuracy \cr
	  \tab \tab  \code{"categorical"} \tab \code{1}    	\tab permutation 	\tab \tab Mean Decrease Accuracy \cr
	  \tab \tab  \code{"continuous"}  \tab \code{2}    	\tab node impurity	\tab \tab Not Available \cr
	  \tab \tab  \code{"binary"}      \tab \code{2}    	\tab node impurity	\tab \tab Mean Decrease in AUC \cr
	  \tab \tab  \code{"categorical"} \tab \code{2}    	\tab node impurity	\tab \tab Not Available   } 

For binary CF models, if{importance.type = 2}, function uses AUC-based variables importances as described by Janitza et al. (2012). Here, the area under the curve instead of the accuracy is used to calculate the importance of each variable. This AUC-based variable importance measure is more robust towards class imbalance. 

Also, for CF models, if \code{cf.conditional = TRUE}, the importance of each variable is computed by permuting within a grid defined by the covariates that are associated (with 1 - p-value greater than threshold) to the variable of interest. The resulting variable importance score is conditional in the sense of beta coefficients in regression models, but represents the effect of a variable in both main effects and interactions. See Strobl et al. (2008) for details. Conditional improtance can be slow for large datasets.

For SGB models:
\tabular{lllclll}{
	  \tab \tab  \code{response type} \tab \code{type}	\tab				\tab \tab Importance Measure \cr
	  \tab \tab  \code{"continuous"}  \tab \code{1}		\tab permutation 		\tab \tab reduction predictive performance \cr
	  \tab \tab  \code{"binary"}      \tab \code{1}		\tab permutation		\tab \tab reduction predictive performance \cr
	  \tab \tab  \code{"continuous"}  \tab \code{2}		\tab gradient of loss function	\tab \tab reduction of squared error \cr
	  \tab \tab  \code{"binary"}      \tab \code{2}		\tab gradient of loss function	\tab \tab reduction in sum of squared error }
}

\value{
The function returns a two element list: \code{IMP1} is the variable importance for \code{model.obj.1}; and, \code{IMP2} is the variable importance for \code{model.obj.2}. This is mostly intended for CF models, where calculating the conditional importance can represent a considerable time investment. For other model types it would be just as easy to recalcuate importances on the fly as needed.
}

\note{
Importance currently unavailable for QRF models.}

\references{

Breiman, L. (2001) Random Forests. Machine Learning, 45:5-32.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. Ann. Stat., 29(5):1189-1232.

Alexander Hapfelmeier, Torsten Hothorn, Kurt Ulm, and Carolin Strobl (2012). A New Variable Importance Measure for Random Forests with Missing Data. Statistics and Computing, http://dx.doi.org/10.1007/s11222-012-9349-1 

Torsten Hothorn, Kurt Hornik, and Achim Zeileis (2006b). Unbiased Recursive Partitioning: A Conditional Inference Framework. Journal of Computational and Graphical Statistics, 15 (3), 651-674. Preprint available from http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf 

Silke Janitza, Carolin Strobl and Anne-Laure Boulesteix (2013). An AUC-based Permutation Variable Importance Measure for Random Forests. BMC Bioinformatics.2013, 14 119. http://www.biomedcentral.com/1471-2105/14/119 

Carolin Strobl, Anne-Laure Boulesteix, Thomas Kneib, Thomas Augustin, and Achim Zeileis (2008). Conditional Variable Importance for Random Forests. BMC Bioinformatics, 9, 307. http://www.biomedcentral.com/1471-2105/9/307

}

\author{Elizabeth Freeman}


\seealso{\code{\link{model.build}}}

\examples{

###########################################################################
############################# Run this set up code: #######################
###########################################################################

# set seed:
seed=38

# Define training and test files:

qdata.trainfn = system.file("extdata", "helpexamples","DATATRAIN.csv", package = "ModelMap")

# Define folder for all output:
folder=getwd()	


#identifier for individual training and test data points

unique.rowname="ID"

##################################################################
########## Continuous Response, Continuous Predictors ############
##################################################################

#file names:
MODELfn.RF="RF_Bio_TC"				
MODELfn.SGB="SGB_Bio_TC"

#predictors:
predList=c("TCB","TCG","TCW")	

#define which predictors are categorical:
predFactor=FALSE	

# Response name and type:
response.name="BIO"
response.type="continuous"

########## Build Models #################################

model.obj.RF = model.build( model.type="RF",
                       qdata.trainfn=qdata.trainfn,
                       folder=folder,		
                       unique.rowname=unique.rowname,		
                       MODELfn=MODELfn.RF,
                       predList=predList,
                       predFactor=predFactor,
                       response.name=response.name,
                       response.type=response.type,
                       seed=seed
)

model.obj.SGB = model.build( model.type="SGB",
                       qdata.trainfn=qdata.trainfn,
                       folder=folder,		
                       unique.rowname=unique.rowname,		
                       MODELfn=MODELfn.SGB,
                       predList=predList,
                       predFactor=predFactor,
                       response.name=response.name,
                       response.type=response.type,
                       seed=seed+1
)

############## Make Imortance Plot - RF vs. SGB ###################

model.importance.plot(	model.obj.1=model.obj.RF, 
			model.obj.2=model.obj.SGB, 
			model.name.1="RF Model", 
			model.name.2="SGB Model", 
			scale.by="sum",
			sort.by="predList", 
			predList=predList,
			main="RF verses SGB",
			device.type="default")

########## Make Imortance Plot - RF Importance type 1 vs 2 #######

model.importance.plot(	model.obj.1=model.obj.RF, 
			model.obj.2=model.obj.RF, 
			model.name.1="PercentIncMSE", 
			model.name.2="IncNodePurity",
			imp.type.1=1,
			imp.type.2=2,
			scale.by="sum",
			sort.by="predList", 
			predList=predList,
			main="Imp type 1 vs Imp type 2",
			device.type="default")


\dontrun{
##################################################################
########## Categorical Response, Continuous Predictors ###########
##################################################################

file name:
MODELfn="RF_NLCD_TC"				

predictors:
predList=c("TCB","TCG","TCW")	

define which predictors are categorical:
predFactor=FALSE	

 Response name and type:
response.name="NLCD"
response.type="categorical"

########## Build Model #################################

model.obj.NLCD = model.build( model.type="RF",
                       qdata.trainfn=qdata.trainfn,
                       folder=folder,		
                       unique.rowname=unique.rowname,	
                       MODELfn=MODELfn,
                       predList=predList,
                       predFactor=predFactor,
                       response.name=response.name,
                       response.type=response.type,
                       seed=seed)

############## Make Imortance Plot ###################

model.importance.plot(	model.obj.1=model.obj.NLCD, 
			model.obj.2=model.obj.NLCD, 
			model.name.1="NLCD=41", 
			model.name.2="NLCD=42",
			class.1="41",
			class.2="42",
			scale.by="sum",
			sort.by="predList", 
			predList=predList,
			main="Class 41 vs. Class 42",
			device.type="default")


##################################################################
############## Conditonal inference forest models ################
##################################################################

#predictors:
predList=c("TCB","TCG","TCW","NLCD")

#define which predictors are categorical:
predFactor=c("NLCD")

#binary response
response.name="CONIFTYP"
response.type="binary"
MODELfn.CF="CF_CONIFTYP_TCandNLCD"

####################### Build Model ##############################

model.obj.CF = model.build( model.type="CF",
                       qdata.trainfn=qdata.trainfn,
                       folder=folder,		
                       unique.rowname=unique.rowname,		
                       MODELfn=MODELfn.CF,
                       predList=predList,
                       predFactor=predFactor,
                       response.name=response.name,
                       response.type=response.type,
                       seed=seed
)

################## Make Imortance Plot ##########################

#Conditional vs. Unconditional importance#

model.importance.plot(	model.obj.1=model.obj.CF, 
			model.obj.2=model.obj.CF, 
			model.name.1="conditional",
			model.name.2="unconditional",
			imp.type.1=1,
			imp.type.2=1, 
			cf.conditional.1=TRUE,
			cf.conditional.2=FALSE,
			scale.by="sum",
			sort.by="predList", 
			predList=predList,
			main="Conditional verses Unconditional",
			device.type="default"
)

} # end dontRun
}

\keyword{ models }

