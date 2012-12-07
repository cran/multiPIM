\name{summary.multiPIM}
\alias{summary.multiPIM}
\alias{print.summary.multiPIM}
\title{Summary methods for class multiPIM}
\description{Generate and print summaries of \code{"multiPIM"} objects (which result from calling either the \code{\link{multiPIM}} or the \code{\link{multiPIMboot}} function). Summaries may be of type "statistical", "time" or "both" (default). Statistical summaries contain, for each exposure-outcome pair, the parameter estimate, the standard error, the test statistic, the unadjusted p-value, and the Bonferroni-adjusted p-value. Time summaries contain a breakdown by g vs. Q modeling, and (if super learning was used to generate the \code{"multiPIM"} object) by super learner candidate, of the time taken to run \code{multiPIM}.}

\usage{
\method{summary}{multiPIM}(object,
                type = c("both", "statistical", "time"),
                use.plug.in.se = is.null(object$boot.param.array),
                alternative.se.matrix = NULL,
                two.sided.p.vals = TRUE,
                bf.multiplier = object$num.exp * object$num.out,
                by.exposure = TRUE,
                digits = 4,
                \dots)

\method{print}{summary.multiPIM}(x, by.exposure, digits, ...)
}

\arguments{
  \item{object}{an object of class \code{"multiPIM"} (the result of a call to \code{\link{multiPIM}} or \code{\link{multiPIMboot}}) to be summarized.}
  \item{type}{The type of summary required. Default is to include both statistical and timing information in the summary.}
  \item{use.plug.in.se}{logical value indicating whether the plug-in standard errors are to be used. Defaults to \code{TRUE} when \code{object} is the result of a call to \code{\link{multiPIM}}, and to \code{FALSE} when \code{object} is the result of a call to \code{\link{multiPIMboot}}, in which case the bootstrap standard errors are used to calculate test statistics and p-values. This argument is ignored when \code{alternative.se.matrix} is non-null.}
  \item{alternative.se.matrix}{matrix of standard errors which were obtained by the user through some method other than the normal plug-in (from \code{\link{multiPIM}})  or bootstrap (from \code{\link{multiPIMboot}}). Must have \code{dim} attribute equal to that of \code{object$param.estimates} and corresponding elements will be assumed to match up.}
  \item{two.sided.p.vals}{logical value. Should p-values be two-sided (for \code{TRUE}), or one-sided (for \code{FALSE}).}
  \item{bf.multiplier}{what number should p-values be multiplied by in doing the Bonferroni Correction. Defaults to the number of exposure-outcome pairs.}
  \item{by.exposure}{logical value. If the summary is printed, and if there is more than one exposure AND more than one outcome, should each table correspond to the exposure-outcome pairs involving each exposure (for \code{TRUE}, the default), or each outcome (for \code{FALSE}). The \code{by.exposure} argument to the print method, if given, will override the one specified when the summary was generated with the summary method.}
  \item{digits}{passed to \code{\link{print.default}} when and if the summary object is printed. The \code{digits} argument to the print method, if given, will override the one specified when the summary was generated with the summary method.}
  \item{x}{an object of class \code{"summary.multiPIM"} (the result of a call to the summary method) to be printed.}
  \item{\dots}{for the summary method: currently ignored. For the print method: passed to \code{\link{print.default}}.}
}

\author{Stephan Ritter, with design contributions from Alan Hubbard and Nicholas Jewell.}

\value{

For the summary method:

An object of class \code{"summary.multiPIM"}, which will have different slots depending on the value of the \code{type} argument.

The object will always have the following slots (regardless of the value of \code{type}):

\item{type}{the value of the type argument to summary.} 
\item{digits}{integer value which will be the default value passed to \code{\link{print.default}} when this object is printed.}
\item{call}{a copy of the call which was used to generate the \code{"multiPIM"} object on which the summary method was called.}

If \code{type} is \code{"statistical"} or \code{"both"}, the result of summary will in addition have the following slots: 

\item{summary.array}{a three-dimensional array containing the information which can be used to build summary tables (see below).}
\item{two.sided.p.vals}{logical value indicating whether p-values used are two-sided (for \code{TRUE}) or one-sided (for \code{FALSE}).}
\item{stand.err.type}{the type of standard error which has been used to generate this summary object: either \code{"plug.in"}, \code{"bootstrap"}, or \code{"alternative"}.}
\item{bf.multiplier}{the value of the bf.multiplier argument}
\item{by.exposure}{logical value which will be used by default when this object is printed to decide whether the tables should be arranged by exposure (for \code{TRUE}) or by outcome (for \code{FALSE}).}

Details for the \code{summary.array} slot: the first dimension corresponds to the exposures (columns of A from the \code{"multiPIM"} object for which the summary is being generated), the second dimension to the outcomes (columns of Y) and the third dimension has length 5 and corresponds to the 5 relevant attributes for each exposure-outcome pair (i.e. the parameter estimate, the standard error of that estimate, the test statistic, the unadjusted p-value and the Bonferroni-adjusted p-value, in that order). Thus, \code{summary.array[1,2,3]} would be the test statistic for the pair consisting of the first exposure (first column of A) and the 2nd outcome (2nd column of Y), while \code{summary.array[3,2,1]} would be the parameter estimate for the pair consisting of the 3rd exposure and the 2nd outcome. To access the matrix containing all unadjusted p-values, use \code{summary.array[,,4]}, to access the matrix consisting of everything that involves the fourth outcome use \code{summary.array[,4,]}.

If \code{type} is \code{"time"} or \code{"both"}, the result of summary will in addition have the following slots: 

\item{main.time}{total time (in seconds) taken to generate the multiPIM result which is being summarized.}
\item{g.time}{time in seconds taken for running g models.}
\item{Q.time}{time in seconds taken for running Q models.}
\item{g.Q.time.frame}{data frame containing breakdown of total time by g vs. Q modeling, with seconds and percentages}
\item{g.sl.time}{if g.method is "sl", time in seconds taken for running cross-validation of g models.}
\item{Q.sl.time}{if Q.method is "sl", time in seconds taken for running cross-validation of Q models.}
\item{g.sl.xval.time.mat}{if super learning was used for g, a matrix containing a breakdown by super learner candidate of the time taken for cross-validation of g models, with seconds and percentages.}
\item{Q.sl.xval.time.mat}{if super learning was used for Q, a matrix containing a breakdown by super learner candidate of the time taken for cross-validation of Q models, with seconds and percentages.}

Print method:

The print method returns its first argument (\code{x}, which should be an object of class \code{"summary.multiPIM"}) invisibly.
}

\seealso{
\code{\link{multiPIM}} and \code{\link{multiPIMboot}}
}

\examples{

## load example from multiPIM help file

example(multiPIM)

## The results can also be displayed by outcome instead of by exposure:

summary(result, by.exposure = FALSE)
## now each table corresponds to all the pairs involving a single outcome

## may be best to store the summary object

sum.obj <- summary(result, by.exposure = FALSE)
sum.obj

## now the print method can be used to overide the values for
## by.exposure and digits (but not the other arguments):

print(sum.obj, by.exposure = TRUE, digits =  3)

## also can hand pick the info that we want from the summary.array slot
## e.g. let's say we are interested in all of the standard errors:

sum.obj$summary.array[,,2]

## or we are only interested in the exposure1-outcome2 pair:

sum.obj$summary.array[1,2,]

## or by name

sum.obj$summary.array["A1","Y2",]
}
