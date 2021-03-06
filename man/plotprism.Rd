% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotprism.R
\name{plotprism}
\alias{plotprism}
\title{Partially Recursively Induced Structured Moderation}
\usage{
plotprism(l, color = 2, extra = 1, type = 1, digit = 4)
}
\arguments{
\item{l}{An rpart output used to plot prism or hprism tree}

\item{color}{Color of the boxes around the text. Defalt value in 2.}

\item{extra}{Display extra information at the nodes. Default value is 1.}

\item{type}{Type of plot. Defual value is 1.}

\item{digit}{number of significant digits in displayed numbers. Default value is 4.}
}
\description{
Plot the prism or hprism model
}
\details{
This function calls the prp R function in the rpart.plot package which provides prettier more customizable tree plots. Details about how to customize the plot by specifying different values of parameters extra and type can be found in \code{\link{prp}} help file in rpart.plot package.
}
\examples{
\dontrun{
plotprism(l.prism)
plotprism(l.hprism)
plotprism(l.hprism1)
}
}
