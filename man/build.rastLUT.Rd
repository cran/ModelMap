\name{build.rastLUT}

\alias{build.rastLUT}

\title{Build a raster Look-UP-Table for training dataset}

\description{
  GUI prompts will help the user build a Look-Up-Table to associat predictor variable with their corresponding spatial rasters.
}

\usage{
build.rastLUT(imageList=NULL,predList=NULL,qdata.trainfn=NULL,rastLUTfn=NULL,folder=NULL)
}

\arguments{
 \item{imageList}{ Vector. A vector of character strings giving names and full paths to all raster data files used in model. }
  \item{predList}{ Vector. A vector of character strings giving the predictor names used as headers in the model training data.}
 \item{qdata.trainfn}{String.  The name (full path or base name with path specified by \code{folder}) of the training data file used for building the model. The file must be a comma-delimited file \code{*.csv} with column headings. \code{qdata.trainfn} can also be an \code{R} dataframe. The column headers from \code{qdata.trainfn} are used to generate a list of possible predictors for the raster Look-UP-Table.}
  \item{rastLUTfn}{ String. The name of the file output for the Look-Up-Table. By default, if a file name is provided by the \code{"qdatatrainfn"} argument \code{"_rastLUT.csv"} appended after \code{"qdatatrainfn"}. Otherwise, default filename for look-up-table is\code{"rastLUT.csv"}}
  \item{folder}{ String.  The folder used for output.  Do not add ending slash to path string.  If \code{folder = NULL} (default), a GUI interface prompts user to browse to a folder.  To use the working directory, specify \code{folder = getwd()}. }
}
\details{
This function helps the user create a raster Look-Up-Table to be used later by \code{model.mapmake()}. Currently this function only works in a Windows environment. 

First, if \code{"folder"} is not given, the user selects the output folder for the Look-UP-Table.

Second, if \code{"predList"} or \code{"qdatatrainfn"}  are not given, the user selects the file containing the training data. The header of the file is used to generate a selection list of possible predictor variables.

Third, if \code{"imageList"} is not provided, the user selects the rasters.

Finally, the function steps through each band of each raster, and the user selects the apropriate predictor.

}
\value{
Returns a data frame containing the raster Look-Up-Table. Also Writes a \code{.csv} file containing the raster Look-Up-Table.
}
\references{}
\author{ Elizabeth Freeman }
\note{  

}
\seealso{ }
\examples{

folder<-system.file("external", "helpexamples", package = "ModelMap")
qdata.trainfn = paste(folder,"/DATATRAIN.csv",sep="")

#build.rastLUT(	qdata.trainfn=qdata.trainfn,
#		folder=folder)


}

\keyword{ models }

