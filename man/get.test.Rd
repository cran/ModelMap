\name{get.test}
\alias{get.test}

\title{Randomly Divide Data into Training and Test Sets}
\description{
  Uses random selection to split a dataset into training and test data sets
}
\usage{
get.test(proportion.test, qdatafn = NULL, seed = NULL, folder=NULL, 
qdata.trainfn = paste(strsplit(qdatafn, split = ".csv")[[1]], "_train.csv", sep = ""), 
qdata.testfn = paste(strsplit(qdatafn, split = ".csv")[[1]], "_test.csv", sep = ""))
}

\arguments{
  \item{proportion.test}{ Number. The proportion of the training data that will be randomly extracted for use as a test set. Value between 0 and 1. }
  \item{qdatafn}{ String. The name (basename or full path) of the data file to be split into training and test data. This data should include both response and predictor variables. The file must be a comma-delimited file \code{*.csv)} with column headings and the predictor names in the file must match the raster layer files, if applying predictions (\code{predict = TRUE}). If \code{NULL} (the default), a GUI interface prompts user to browse to the data file. }
  \item{seed}{ Integer. The number used to initialize randomization to randomly select rows for a test data set. If you want to produce the same model later, use the same seed. If \code{seed = NULL} (the default), a new one is created each time. }
  \item{folder}{ String.  The folder used for all output from predictions and/or maps.  Do not add ending slash to path string.  If \code{folder = NULL} (default), a GUI interface prompts user to browse to a folder.  To use the working directory, specify \code{folder = getwd()}. }
  \item{qdata.trainfn}{ String. The name of the file output of training data. By default, \code{_train} appended after \code{qdatafn}.}
  \item{qdata.testfn}{ String. The name of the file output of test data. By default, \code{_test} appended after \code{qdatafn}. }
}
\details{
This function should be run once, before starting analysis to create training and test sets. If the cross validation option is to be used with RF or SGB models, or if the OOB option is to be used for RF models, then this step is unnecessary. 
}
\value{
Outputs a training data file and test data file. Unless \code{qdata.trainfn} or \code{qdata.testfn} are specified, the output will be located in \code{folder}. The output will have the same rows and columns as the original data.
}

\author{ Elizabeth Freeman }

\examples{

\dontrun{
qdatafn<-system.file("extdata", "helpexamples","DATATRAIN.csv", package = "ModelMap")

qdata<-read.table(file=qdatafn,sep=",",header=TRUE,check.names=FALSE)

get.test(	proportion.test=0.2,
		qdatafn=qdatafn,
		seed=42,
		folder=getwd(),
		qdata.trainfn="example.train.csv",
		qdata.testfn="example.test.csv")
} # end dontrun

}

\keyword{ models }

