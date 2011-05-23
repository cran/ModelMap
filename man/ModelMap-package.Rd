\name{ModelMap-package}
\alias{ModelMap-package}
\alias{ModelMap}
\docType{package}
\title{
Modeling and Map production using Random Forest and Stochastic Gradient Boosting
}
\description{
The \code{ModelMap} package creates sophisticated models of training data and validates the models with an independant test set, cross validation, or in the case of Random Forest Models, with Out OF Bag (OOB) predictions on the training data. The package creates graphs and tables of the model validation results. It applies these models to GIS .img files of predictors to create detailed prediction surfaces. The package handles large predictor files for map making by reading in the predictor \code{.img} files in chunks, and writing the predictions for each data chunk to the \code{.txt} file, before reading the next chunk of predictor data.
}
\details{
\tabular{ll}{
Package: \tab ModelMap\cr
Type: \tab Package\cr
Version: \tab 2.3.1\cr
Date: \tab 2011-05-23
\cr
License: \tab Unlimited. This code was written and prepared by a U.S. Government employee on official time, 
and therefore it is in the public domain and not subject to copyright.\cr
}

This package provides a push button appraoch to complex model building and production mapping. It contains five functions. 

It contains two simple functions that can be used before beginning the model building process: \link{get.test} that can be used to radomly divide a training dataset into training and test/validation sets; and \link{build.rastLUT} that uses GUI prompts to walk a user through the process of setting up a Raster look up table to link predictors from the training data with the rasters used for map contruction.

It also contains four functions used for the model building and map contruction process itself: \link{model.build}, \link{model.diagnostics},\link{model.interaction.plot} and \link{model.mapmake}.

\code{ModelMap} can be run in a traditional R command mode, where all arguments are specified in the function call. However it can also be used in a full push button mode, where you type in the simple command such as \link{model.build}, and GUI pop up windows ask questions about the type of model, the file locations of the data, etc...

Random Forest is implemented through the \code{randomForest} package within \code{R}. Random Forest is more user friendly than Stochastic Gragient Boosting, as it has fewer parameters to be set by the user, and is less sensitive to tuning of these parameters. A Random Forest model consists of multiple trees that vote on predictions. For each tree a random subset of the training data is used to construct the tree, with the remaining data points used to construct out-of-bag (OOB) error estimates. At each node of the tree a random selection of predictors is chosen to determine the split. The number of predictors used to select the splits is the primary user specified parameter that can affect model performance, and this parameter can be automatically optimized using the \code{randomForest} function \code{tuneRF()}. Random Forest will not over fit data, therefore the only penalty of increasing the number of trees is computation time. Random Forest can compute variable importance, an advantage over some "black box" modeling techniques if it is important to understand the ecological relationships underlying a model (Brieman, 2001).

Stochastic gradient boosting (Friedman 2001, 2002), is related to both boosting and bagging. Many small classification or regression trees are built sequentially from "pseudo"-residuals (the gradient of the loss function of the previous tree). At each iteration, a tree is built from a random sub-sample of the dataset (selected without replacement) and an incremental improvement in the model. Using only a fraction of the training data increases both the computation speed and the prediction accuracy, while also helping to avoid over-fitting the data. An advantage of stochastic gradient boosting is that it is not necessary to pre-select or transform predictor variables. It is also resistant to outliers, as the steepest gradient algorithm emphasizes points that are close to their correct classification. Stochastic gradient boosting is implemented through the \code{gbm} package within \code{R}.  One disadvantege of Stochastic Gradient Boosting, compared to Random Forest, is increased number of user specified parameters, and the SGB models tend tobe more sensitive to these parameters. Model fitting parameters include distribution, interaction depth, bagging fraction, shrinkage rate, and training fraction. These parameters can be set in the argument list when calling \code{model.map()}. Values for these parameters other than the defaults can not be set by point and click in the GUI pop up windows. Friedman (2001, 2002) and Ridgeway (1999) provide guidelines on appropriate settings for model fitting options.

For Presence-Absence data, the package \code{PresenceAbsence} is used for model validation.

For map making, the \code{rgdal} is used to read \code{.img} files.

For interaction plots, the \code{fields} package is used to produce image plots.
}
\author{
Author: Elizabeth Freeman and Tracey Frescino

Maintainer: Elizabeth Freeman <eafreeman@fs.fed.us>

}
\references{

Breiman, L. (2001) Random Forests. Machine Learning, 45:5-32.

Elith, J., Leathwick, J. R. and Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology. 77:802-813.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. Ann. Stat., 29(5):1189-1232.

Friedman, J.H. (2002). Stochastic gradient boosting. Comput. Stat. Data An., 38(4):367-378.

Liaw, A. and  Wiener, M. (2002). Classification and Regression by randomForest. R News 2(3), 18--22.

Ridgeway, G., (1999). The state of boosting. Comp. Sci. Stat. 31:172-181
}

\keyword{ package }

