\name{model.explore}
\alias{model.explore}

\title{
Exploratory data analysis
}
\description{
Graphically explores the relationships between the training data and the predictor rasters.
}
\usage{
model.explore(qdata.trainfn = NULL, folder = NULL, predList = NULL, 
predFactor = FALSE, response.name = NULL, response.type = NULL, 
response.colors = NULL, unique.rowname = NULL, OUTPUTfn = NULL, 
device.type = NULL, allow.default.graphics=FALSE, res=NULL, jpeg.res = 72, 
MAXCELL=100000, device.width = NULL, device.height = NULL, units="in", 
pointsize=12, cex=1, rastLUTfn = NULL, create.extrapolation.masks = FALSE, 
na.value = -9999, col.ramp = rainbow(101, start = 0, end = 0.5), 
col.cat = palette()[-1])
}

\arguments{
  \item{qdata.trainfn}{String.  The name (full path or base name with path specified by \code{folder}) of the training data file used for building the model (file should include columns for both response and predictor variables).  The file must be a comma-delimited file \code{*.csv} with column headings. \code{qdata.trainfn} can also be an \code{R} dataframe. If predictions will be made (\code{predict = TRUE} or \code{map=TRUE}) the predictor column headers must match the names of the raster layer files, or a \code{rastLUT} must be provided to match predictor columns to the appropriate raster and band. If \code{qdata.trainfn = NULL} (the default), a GUI interface prompts user to browse to the training data file.}

  \item{folder}{String.  The folder used for all output from predictions and/or maps.  Do not add ending slash to path string.  If \code{folder = NULL} (default), a GUI interface prompts user to browse to a folder.  To use the working directory, specify \code{folder = getwd()}.}

  \item{predList}{ String.  A character vector of the predictor short names used to build the model.  These names must match the column names in the training/test data files and the names in column two of the \code{rastLUT}.  If \code{predList = NULL} (the default), a GUI interface prompts user to select predictors from column 2 of \code{rastLUT}.}

  \item{predFactor}{ String.  A character vector of predictor short names of the predictors from \code{predList} that are factors (i.e categorical predictors).  These must be a subset of the predictor names given in \code{predList} Categorical predictors may have multiple categories.}

  \item{response.name}{ String.  The name of the response variable used to build the model.  If \code{response.name = NULL}, a GUI interface prompts user to select a variable from the list of column names from training data file.  \code{response.name} must be column name from the training/test data files. }

  \item{response.type}{ String.  Response type: \code{"binary"}, \code{"categorical"} or \code{"continuous"}.  Binary response must be binary 0/1 variable with only 2 categories.  All zeros will be treated as one category, and everything else will be treated as the second category.}

  \item{response.colors}{ Data frame. A two column data frame. Column names must be:\code{category}, the response categories; and, \code{color}, the colors associated with each category.}

  \item{unique.rowname}{ String.  The name of the unique identifier used to identify each row in the training data. If \code{unique.rowname = NULL}, a GUI interface prompts user to select a variable from the list of column names from the training data file.  If \code{unique.rowname = FALSE}, a variable is generated of numbers from \code{1} to \code{nrow(qdata)} to index each row. }

  \item{OUTPUTfn}{ String. Filename that ouput file names will be based on.} 

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

Note that the \code{"default"} device is disabled unless \code{allow.default.graphics=TRUE}. This is because these graphics are slow to produce, and if the onscreen graphics window is moved or closed while the function is in progress there is a risk of crashing the entire R session. }

  \item{allow.default.graphics}{ Logical.  Should the default on-screen graphics device be allowed. USE WITH CAUTION! These graphics are complicated and slow to produce. If the on-screen default graphics device is moved or closed before the plot is completed it can crash the entire R session.}

  \item{res}{ Integer.  Model validation.  Pixels per inch for jpeg, png, and tiff plots.  The default is 72dpi, good for on screen viewing. For printing, suggested setting is 300dpi. }

  \item{jpeg.res}{ Integer. Graphical output. Deprecated. Ignored unless \code{res} not provided. }

  \item{MAXCELL}{ Integer.  Graphical output.  The maximum number of raster cells used to create the graphical output. Rasters larger than this value will be subsampled for the graphical maps and figures. The default value of \code{MAXCELL=100000} is generally a good resolution for onscreen viewing with the default jpeg resolution of 72dpi. Publication quality qraphics may require higher \code{MAXCELL}. Higher values require more memory and are slower to process.

Note: \code{MAXCELL} only affects graphical figures. Output rasters generated when \code{create.extrapolation.masks=TRUE} are always done on full resolution rasters. }

  \item{device.width}{ Integer.  Model validation.  The device width for diagnostic plots in inches. }

  \item{device.height}{ Integer.  Model validation.  The device height for diagnostic plots in inches. }

  \item{units}{ Model validation.  The units in which \code{device.height} and \code{device.width} are given. Can be \code{"px"} (pixels), \code{"in"} (inches, the default), \code{"cm"} or \code{"mm"}. }

  \item{pointsize}{ Integer.  Model validation.  The default pointsize of plotted text, interpreted as big points (1/72 inch) at \code{res} ppi}

  \item{cex}{ Integer.  Model validation.  The cex for diagnostic plots. }

 \item{rastLUTfn}{ String.  The file name (full path or base name with path specified by \code{folder}) of a \code{.csv} file for a \code{rastLUT}. Alternatively, a dataframe containing the same information. The \code{rastLUT} must include 3 columns: (1) the full path and name of the raster file; (2) the shortname of each predictor / raster layer (band); (3) the layer (band) number.  The shortname (column 2) must match the names \code{predList}, the predictor column names in training/test data set (\code{qdata.trainfn} and \code{qdata.testfn}, and the predictor names in \code{model.obj}. 

Example of comma-delimited file:

\tabular{llllll}{
	  \tab \tab \tab \code{C:/button_test/tc99_2727subset.img,} \tab \code{tc99_2727subsetb1,} \tab \code{1}\cr
	  \tab \tab \tab \code{C:/button_test/tc99_2727subset.img,} \tab \code{tc99_2727subsetb2,} \tab \code{2}\cr
	  \tab \tab \tab \code{C:/button_test/tc99_2727subset.img,} \tab \code{tc99_2727subsetb3,} \tab \code{3}}}

  \item{create.extrapolation.masks}{Logical.  If \code{TRUE} then the raster brick containing the masks for all predictors from \code{predList} is saved as image file. The layers in this file will be in the same order as the predictors in \code{predList}   }

  \item{na.value}{Value used in rasters to indicate \code{NA}. Note this value is only used  for \code{NA} values in the predictor rasters.  Note: all predictor rasters must use the same value for \code{NA}. \code{NA} values in the training data should be indicated with \code{NA}. }

  \item{col.ramp}{Color ramp to use for continuous predictors}

  \item{col.cat}{Vector. Vector of colors to use for categorical predictors. }

}
\details{
The \code{model.explore} function is intended to aid with preliminary data exploration before model building. It includes graphical tools to explore the relationships between the training data (both predictors and responses) as well as the predictor rasters. It uses the \code{corrplot} package to create a correlation plot of the continuous predictor. This can aid in interpreting the \code{\link{model.importance.plot}} output from the models, as Random Forest models divide importance between correlated predictors, while Stochastic Gradient Boosting models assing the majority of the importance to the correlated predictor that is used earlies in the model.

The \code{model.explore} function also can aid in identifying if additional training data is needed. For example, the maps of the extrapolation masks for the predictor rasters help spot areas of the map where the pixels lie outside the range of the training data, and therefore any model predictions will be extrapolations, and possibly unreliable. The user can decide to either collect additional training data, or mask out these areas of the final prediction output of \code{\link{model.mapmake}}.

To increase speed, the default behavior for large predictor rasters is to create the graphics from subsampled rasters. (Note: for categorical predictors, the full raster is always used to identify all categories found in the map area.) If \code{create.extrapolation.masks=TRUE}, then the full rasters are used for the extrapolation masks, regardless of size of the reasters. This option runs much slower, as large rasters need to be read into R a block at a time.
}

\value{
Function does not return a value, but does create files. 

Graphical files are created for each predictor variable, with file type determined by \code{device.type}. In addition, if \code{create.extrapolation.masks}, an extrapolation mask raster is created for each predictor as well as an overall extrapolation mask, with the value \code{1} for pixels with predictor values within the range of the training data, or categories found in the training data, and the value \code{0} for pixels outside the range of the training data, categories not found in the training data, or NA value. The overall extrapolation mask has \code{0} if any of the predictors for that pixel are extrapolated. Note that this option is much slower to run.

}

\author{Elizabeth Freeman}

\note{
The default graphics device is disabled unless \code{allow.default.graphics} is set to \code{TRUE}. These graphics can be slow to produce, and if the on screen graphics device is moved or closed while the graphic is in progress, it can crash R. It is recomended that graphics be written to a file by using jpeg, pdf, etc... \code{device.type}.
}


\examples{

\dontrun{

###########################################################################
############################# Run this set up code: #######################
###########################################################################

###Define training and test files:
qdata.trainfn = system.file("extdata", "helpexamples","DATATRAIN.csv", package = "ModelMap")

###Define folder for all output:
folder=getwd()	

###identifier for individual training and test data points
unique.rowname="ID"			

###predictors:
predList=c("TCB","TCG","TCW","NLCD")

###define which predictors are categorical:
predFactor=c("NLCD")

###Create a the filename (including path) for the rast Look up Tables ###
rastLUTfn.2001 <- system.file(	"extdata",
				"helpexamples",
				"LUT_2001.csv",
				package="ModelMap")

###Load rast LUT table, and add path to the predictor raster filenames in column 1 ###
rastLUT.2001 <- read.table(rastLUTfn.2001,header=FALSE,sep=",",stringsAsFactors=FALSE)

for(i in 1:nrow(rastLUT.2001)){
	rastLUT.2001[i,1] <- system.file("extdata",
					"helpexamples",
					rastLUT.2001[i,1],
					package="ModelMap")
}                                 

#################Continuous Response###################

###Response name and type:
response.name="BIO"
response.type="continuous"

###file name to store model:
OUTPUTfn="BIO_TCandNLCD.img"

###run model.explore

model.explore(	qdata.trainfn=qdata.trainfn,
		folder=folder,		
		predList=predList,
		predFactor=predFactor,

		response.name=response.name,
		response.type=response.type,
	
		unique.rowname=unique.rowname,

		OUTPUTfn=OUTPUTfn,
		device.type="jpeg",
		jpeg.res=144,


		# Raster arguments
		rastLUTfn=rastLUT.2001,
		na.value=-9999,

		# colors for continuous predictors
		col.ramp=rainbow(101,start=0,end=.5), 
		# colors for categorical predictors
		col.cat=c("wheat1","springgreen2","darkolivegreen4",
			  "darkolivegreen2","yellow","thistle2",
			  "brown2","brown4")
)
} # end dontrun



}

\keyword{ models }

