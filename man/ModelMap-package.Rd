\name{ModelMap-package}
\alias{ModelMap-package}
\alias{ModelMap}
\docType{package}
\title{
Modeling and Map Production using Random Forest and Related Stochastic Models
}
\description{
Creates sophisticated models of training data and validates the models with an independent test set, cross validation, or with Out Of Bag (OOB) predictions on the training data. Create graphs and tables of the model validation results. Applies these models to GIS .img files of predictors to create detailed prediction surfaces. Handles large predictor files for map making, by reading in the .img files in chunks, and output to the .txt file the prediction for each data chunk, before reading the next chunk of data.
}
\details{
\tabular{ll}{
Package: \tab ModelMap\cr
Type: \tab Package\cr
Version: \tab 3.4.0.4\cr
Date: \tab 2023-04-04
\cr
License: \tab Unlimited. This code was written and prepared by a U.S. Government employee on official time, 
and therefore it is in the public domain and not subject to copyright.\cr
}

This package provides a push button approach to complex model building and production mapping.  It contains three main functions: \code{\link{model.build}},\code{\link{model.diagnostics}}, and \code{\link{model.mapmake}}. 

In addition it contains a simple function \code{\link{get.test}} that can be used to randomly divide a training dataset into training and test/validation sets; \code{\link{build.rastLUT}} that uses GUI prompts to walk a user through the process of setting up a Raster look up table to link predictors from the training data with the rasters used for map contruction; \code{\link{model.explore}}, for preliminary data exploration; and, \code{\link{model.importance.plot}} and \code{\link{model.interaction.plot}} for interpreting the effects of individual model predictors. 

\code{ModelMap} can be run in a traditional R command mode, where all arguments are specified in the function call. However it can also be used in a full push button mode, where you type in the simple command such as \link{model.build}, and GUI pop-up windows ask questions about the type of model, the file locations of the data, etc...

Random Forest is implemented through the \code{randomForest} package within \code{R}. Random Forest is more user friendly than Stochastic Gradient Boosting, as it has fewer parameters to be set by the user, and is less sensitive to tuning of these parameters. A Random Forest model consists of multiple trees that vote on predictions. For each tree a random subset of the training data is used to construct the tree, with the remaining data points used to construct out-of-bag (OOB) error estimates. At each node of the tree a random selection of predictors is chosen to determine the split. The number of predictors used to select the splits is the primary user specified parameter that can affect model performance, and this parameter can be automatically optimized using the \code{randomForest} function \code{tuneRF()}. Random Forest will not over fit data, therefore the only penalty of increasing the number of trees is computation time. Random Forest can compute variable importance, an advantage over some "black box" modeling techniques if it is important to understand the ecological relationships underlying a model (Brieman, 2001).

Quantile Regression Forests is implemented through the \code{quantregForest} package.

Conditional Forests is implemented with the \code{cforest()} function in the \code{party} package. As stated in the \code{party} package, ensembles of conditional inference trees have not yet been extensively tested, so this routine is meant for the expert user only and its current state is rather experimental. 

For Presence-Absence data, the package \code{PresenceAbsence} is used for model validation.

For model diagnostics the package \code{corrplot} is used to plot the correlation between predictor variables.

For map making, the \code{raster} is used to read and write \code{.img} files.

For interaction plots, the \code{fields} package is used to produce image plots.
}
\author{
Author: Elizabeth Freeman and Tracey Frescino

Maintainer: Elizabeth Freeman <elizabeth.a.freeman@usda.gov>

}
\references{
Breiman, L. (2001) Random Forests. Machine Learning, 45:5-32.

Elith, J., Leathwick, J. R. and Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology. 77:802-813.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. Ann. Stat., 29(5):1189-1232.

Friedman, J.H. (2002). Stochastic gradient boosting. Comput. Stat. Data An., 38(4):367-378.

Liaw, A. and  Wiener, M. (2002). Classification and Regression by randomForest. R News 2(3), 18--22.

N. Meinshausen (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7, 983-999 http://jmlr.csail.mit.edu/papers/v7/

Ridgeway, G., (1999). The state of boosting. Comp. Sci. Stat. 31:172-181

Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007). Bias in Random Forest variable Importance Measures: Illustrations, Sources and a Solution. BMC Bioinformatics, 8, 25. http://www.biomedcentral.co,/1471-2105/8/25

Carolin Strobl, James Malley and Gerhard Tutz (2009). An Introduction to Recursive Partitioning: Rationale, Application, and Characteristics of Classification and Regression Trees, Bagging, and Random forests. Phsycological Methods, 14(4), 323-348. 

Torsten Hothorn, Berthold Lausen, Axel Benner and Martin Radespiel-Troeger (2004). Bagging Survival Trees. Statistics in Medicine, 23(1), 77-91.

Torsten Hothorn, Peter Buhlmann, Sandrine Dudoit, Annette Molinaro and Mark J. ven der Laan (2006a). Survival Ensembles. Biostatistics, 7(3), 355-373. 

Torston Hothorn, Kurt Hornik and Achim Zeileis (2006b). Unbiased Recursive Partitioning: A Conditional Inference Framework. JOurnal of Computational and Graphical Statistics, 15(3), 651-674. Preprint available from http://statmath.wu-wein.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf

}

\keyword{ package }

