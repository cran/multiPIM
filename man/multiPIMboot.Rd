\name{multiPIMboot}
\alias{multiPIMboot}
\title{Bootstrap the multiPIM Function}
\description{This function will run \code{\link{multiPIM}} once on the actual data, then sample with replacement from the rows of the data and run \code{\link{multiPIM}} again (with the same options) the desired number of \code{times}.}
\usage{multiPIMboot(Y, A, W = NULL,
             times = 5000,
             id = 1:nrow(Y),
             multicore = FALSE,
             mc.num.jobs,
             mc.seed = 123,
             estimator = c("TMLE", "DR-IPCW", "IPCW", "G-COMP"),
             g.method = "main.terms.logistic", g.sl.cands = NULL,
             g.num.folds = NULL, g.num.splits = NULL,
             Q.method = "sl", Q.sl.cands = "default",
             Q.num.folds = 5, Q.num.splits = 1,
             Q.type = NULL,
             adjust.for.other.As = TRUE,
             truncate = 0.05,
             return.final.models = TRUE,
             na.action,
             verbose = FALSE,
             extra.cands = NULL,
             standardize = TRUE,
             \dots)}
\arguments{
  \item{Y}{a data frame of outcomes containing only numeric (integer or double) values. See details section of \code{\link{multiPIM}} for the default method of determining, based on the values in \code{Y}, which regression types to allow for modelling Q. Must have unique names.}
  \item{A}{a data frame containing binary exposure variables. \emph{Binary} means that all values must be either 0 (indicating unexposed, or part of target group) or 1 (indicating exposed or not part of target group). Must have unique names.}
  \item{W}{an optional data frame containing possible confounders of the effects of the variables in \code{A} on the variables in \code{Y}. No effect measures will be calculated for these variables. May contain numeric (integer or double), or factor values. Must be left as \code{NULL} if not required. If not NULL, must have unique names.}
  \item{times}{single integer greater than or equal to 2. The number of bootstrap replicates of \code{Y}, \code{A} and \code{W} to generate and pass to \code{\link{multiPIM}}.}
  \item{id}{vector which identifies clusters. If obervations i and j are in the same cluster, then \code{id[i]} should be equal to \code{id[j]}. Bootstrapping will be carried out by sampling with replacement from the clusters. Keeping the default value will result in sampling with replacement from the observations (i.e. no clustering).}
  \item{multicore}{logical value indicting whether bootstrapping should be done using multiple simultaneous jobs (as of multiPIM version 1.3-1 this requires the \emph{\link[parallel]{parallel}} package, which is distributed with R version 2.14.0 or later. For earlier versions of multiPIM, this feature relied on  CRAN packages \emph{multicore} and \emph{rlecuyer}.}
  \item{mc.num.jobs}{number of simultaneous multicore jobs, e.g. if you want to use a quad core processor with hyperthreading, use \code{mc.num.jobs = 8}. This must be specified whenever \code{multicore} is true. Automatic detection of the number of cores is no longer available.}
  \item{mc.seed}{integer value with which to seed the RNG when using parallel processing (internally, \code{\link[base]{RNGkind}} will be called to set the RNG to \code{"L'Ecuyer-CMRG"}). Will be ignored if \code{multicore} is \code{FALSE}. If \code{mulicore} is \code{FALSE}, one \dQuote{should} (depending on the candidates used) be able to get reprodicible results by setting the seed normally (with \code{\link[base]{set.seed}}) prior to running multiPIMboot.}
  \item{estimator}{the estimator to be used. The default is \code{"TMLE"}, for the targeted maximum likelihood estimator. Alternatively, one may specify \code{"DR-IPCW"}, for the Double-Robust Inverse Probability of Censoring-Weighted estimator, or \code{"IPCW"}, for the regular IPCW estimator. If the regular IPCW estimator is selected, all arguments which begin with the letter Q are ignored, since only g (the regression of each exposure on possible confounders) needs to be modeled in this case.}
  \item{g.method}{a length one character vector indicating the regression method to use in modelling g. The default value, \code{"main.terms.logistic"}, is meant to be used with the default TMLE estimator. If a different estimator is used, it is recommended to use super learning by specifying \code{"sl"}. In this case, the arguments \code{g.sl.cands}, \code{g.num.folds} and \code{g.num.splits} must also be specified. Other possible values for the \code{g.method} argument are: one of the elements of the vector \code{\link{all.bin.cands}}, or, if \code{extra.cands} is supplied, one of the names of the \code{extra.cands} list of functions. Ignored if \code{estimator} is \code{"G-COMP"}.}
  \item{g.sl.cands}{character vector of length \eqn{\geq 2}{>= 2} indicating the candidate algorithms that the super learner fits for g should use. The possible values may be taken from the vector \code{\link{all.bin.cands}}, or from the names of the \code{extra.cands} list of functions, if it is supplied. Ignored if \code{estimator} is \code{"G-COMP"}. or if \code{g.method} is not \code{"sl"}. NOTE: The TMLE estimator is recommended, but if one is using either of the IPCW estimators, a reasonable choice is to specify \code{g.method = "sl"} and \code{g.sl.cands = default.bin.cands}.}
  \item{g.num.folds}{the number of folds to use in cross-validating the super learner fit for g (i.e. the v for v-fold cross-validation). Ignored if \code{estimator} is \code{"G-COMP"}, or if \code{g.method} is not \code{"sl"}.}
  \item{g.num.splits}{the number of times to randomly split the data into \code{g.num.folds} folds in cross-validating the super learner fit for g. Cross-validation results will be averaged over all splits. Ignored if \code{estimator} is \code{"G-COMP"}, or if \code{g.method} is not \code{"sl"}.}
  \item{Q.method}{character vector of length 1. The regression method to use in modelling Q. See details to find out which values are allowed. The default value, \code{"sl"}, indicates that super learning should be used for modelling Q. Ignored if \code{estimator} is \code{"IPCW"}.}
  \item{Q.sl.cands}{either of the length 1 character values \code{"default"} or \code{"all"} or a character vector of length \eqn{\geq 2}{>= 2} containing elements of either \code{all.bin.cands} or of \code{all.cont.cands}, or of the names of the \code{extra.cands} list of functions, if it is supplied. See details. Ignored if \code{estimator} is \code{"IPCW"} or if \code{Q.method} is not \code{"sl"}.}
  \item{Q.num.folds}{the number of folds to use in cross-validating the super learner fit for Q (i.e. the v for v-fold cross-validation). Ignored if \code{estimator} is \code{"IPCW"} or if \code{Q.method} is not \code{"sl"}.}
  \item{Q.num.splits}{the number of times to randomly split the data into \code{Q.num.folds} folds in cross-validating the super learner fit for Q. Ignored if \code{estimator} is \code{"IPCW"} or if \code{Q.method} is not \code{"sl"}.}
  \item{Q.type}{either \code{NULL} or a length 1 character vector (which must be either \code{"binary.outcome"} or \code{"continuous.outcome"}). This provides a way to override the default mechanism for deciding which candidates will be allowed for modeling Q (see details). Ignored if \code{estimator} is \code{"IPCW"}.}
  \item{adjust.for.other.As}{a single logical value indicating whether the other columns of \code{A} should be included (for \code{TRUE}) or not (for \code{FALSE}) in the g and Q models used to calculate the effect of each column of \code{A} on each column of \code{Y}. See details. Ignored if \code{A} has only one column.}
  \item{truncate}{either \code{FALSE}, or a single number greater than 0 and less than 0.5 at which the values of g(0, W) should be truncated in order to avoid instability of the estimator. Ignored if \code{estimator} is \code{"G-COMP"}.}
  \item{return.final.models}{single logical value indicating whether final g and Q models should be returned by the function (in the slots \code{g.final.models} and \code{Q.final.models}). Default is \code{TRUE}. If memory is a concern, you will probably want to set this to FALSE. Note that only g and Q models for the main multiPIM run will be returned, not for each of the bootstrap runs.}
  \item{na.action}{currently ignored. If any of \code{Y}, \code{A} or (a non-null) \code{W} has missing values, \code{multiPIMboot} will throw an error.}
  \item{verbose}{single logical value. Should messages about the progress of the evaluation be printed out. Some of the candidate algorithms may print messages even when \code{verbose} is set to \code{FALSE}.}
  \item{extra.cands}{a named list of functions. This argument provides a way for the user to specify his or her own functions to use either as stand-alone regression methods, or as candidates for a super learner. See details section of \code{\link{multiPIM}}.}
  \item{standardize}{should all predictor variables be standardized before certain regression methods are run. Passed to all candidates, but only used by some (at this point, lars, penalized.bin and penalized.cont)}
  \item{\dots}{currently ignored.}
}

\details{
Bootstrap standard errors can be calculated by running the \code{summary} function on the \code{multiPIMboot} result (see \code{link{summary.multiPIM}}).

As of \emph{multiPIM} version 1.3-1, support for multicore processing is through R's \emph{parallel} package (distributed with R as of version 2.14.0).

For more details on how to use the arguments, see the details section for \code{\link{multiPIM}}.
}

\value{
Returns an object of class \code{"multiPIM"} which is identical to the object resulting from running the \code{\link{multiPIM}} function in the original data, except for two slots which are slightly different: the \code{call} slot contains a copy of the original call to \code{multiPIMboot}, and the \code{boot.param.array} slot now contains the bootstrap distribution of the parameter estimates gotten by running \code{\link{multiPIM}} on the bootstrap replicates of the original data. Thus the object returned has the following slots:
  \item{param.estimates}{a matrix of dimensions \code{ncol(A)} by \code{ncol(Y)} with \code{rownames} equal to \code{names(A)} and \code{colnames} equal to \code{names(Y)}, with each element being the estimated causal attributable risk for the exposure given by its row name vs. the outcome given by its column name.}
  \item{plug.in.stand.errs}{a matrix with the same dimensions as \code{param.estimates} containing the corresponding plug-in standard errors of the parameter estimates. These are obtained from the influence curve. Note: plug-in standard errors are not available for \code{estimator = "G-COMP"}. This field will be set to \code{NA} in this case.}
  \item{call}{a copy of the call to \code{multiPIMboot} which generated this object.}
  \item{num.exposures}{this will be set to \code{ncol(A)}.}
  \item{num.outcomes}{this will be set to \code{ncol(Y)}.}
  \item{W.names}{the names attribute of the \code{W} data frame, if one was supplied. If no \code{W} was supplied, this will be \code{NA}.}
  \item{estimator}{the estimator used.}
  \item{g.method}{the method used for modelling g.}
  \item{g.sl.cands}{in case super learning was used for g, the candidates used in the super learner. Will be \code{NA} if \code{g.method} was not \code{"sl"}.}
  \item{g.winning.cands}{if super learning was used for g, this will be a named character vector with \code{ncol(A)} elements. The ith element will be the name of the candidate which "won" the cross validation in the g model for the ith column of \code{A}.}
  \item{g.cv.risk.array}{array with dim attribute \code{c(ncol(A), g.num.splits, length(g.sl.cands))} containing cross-validated risks from super learner modeling for g for each exposure-split-candidate triple. Has informative dimnames attribute. Note: the values are technically not risks, but log likelihoods (i.e. winning candidate is the one for which this is a \emph{max}, not a min).}
  \item{g.final.models}{a list of length \code{nrow(A)} containing the objects returned by the candidate functions used in the final g models (see \link{Candidates}).}
  \item{g.num.folds}{the number of folds used for cross validation in the super learner for g. Will be \code{NA} if \code{g.method} was not \code{"sl"}.}
  \item{g.num.splits}{the number of splits used for cross validation in the super learner for g. Will be \code{NA} if \code{g.method} was not \code{"sl"}.}
  \item{Q.method}{the method used for modeling Q. Will be \code{NA} if \code{double.robust} was \code{FALSE}.}
  \item{Q.sl.cands}{in case super learning was used for Q, the candidates used in the super learner. Will be \code{NA} if \code{double.robust} was \code{FALSE} or if \code{Q.method} was not \code{"sl"}.}
  \item{Q.winning.cands}{if super learning was used for Q, this will be a named character vector with \code{ncol(Y)} elements. The ith element is the name of the candidate which "won" the cross validation in the super learner for the Q model for the ith column of \code{Y}.}
  \item{Q.cv.risk.array}{array with dim attribute \code{c(ncol(A), ncol(Y), Q.num.splits, length(Q.sl.cands))} containing cross-validated risks from super learner modeling for Q. Has informative dimnames attribute. Note: the values will be log likelihoods when \code{Q.type} is \code{"binary.outcome"} (see note above for \code{g.cv.risk.array}), and they will be mean squared errors when \code{Q.type} is \code{"continuous.outcome"}.}
  \item{Q.final.models}{a list of length \code{ncol(A)}, each element of which is another list of length \code{ncol(Y)} containing the objects returned by the candidate functions used for the Q models. I.e. \code{Q.final.models[[i]][[j]]} contains the Q model information for exposure i and outcome j.}
  \item{Q.num.folds}{the number of folds used for cross validation in the super learner for Q. Will be \code{NA} if \code{double.robust} was \code{FALSE} or if \code{Q.method} was not \code{"sl"}.}
  \item{Q.num.splits}{the number of splits used for cross validation in the super learner for Q.  Will be \code{NA} if \code{double.robust} was \code{FALSE} or if \code{Q.method} was not \code{"sl"}.}
  \item{Q.type}{either \code{"continuous.outcome"} or \code{"binary.outcome"}, depending on the contents of \code{Y} or on the value of the \code{Q.type} argument, if supplied.}
  \item{adjust.for.other.As}{logical value indicating whether the other columns of \code{A} were included in models used to calculate the effect of each column of \code{A} on each column of \code{Y}. Will be set to \code{NA} when \code{A} has only one column.}
  \item{truncate}{the value of the \code{truncate} argument. Will be set to NA if estimator was \code{"G-COMP"}.}
  \item{truncation.occured}{logical value indicating whether it was necessary to trunctate. \code{FALSE} when \code{truncate} is \code{FALSE}. Will be set to NA if estimator was \code{"G-COMP"}.}
  \item{standardize}{the value of the \code{standardize} argument.}
  \item{boot.param.array}{a three dimensional array with \code{dim} attribute equal to \code{c(times, ncol(A), ncol(Y))} containing the corresponding parameter estimate for each bootstrap replicatate-exposure-outcome trio. Also has an informative \code{dimnames} attribute for easy printing.}
  \item{main.time}{time (in seconds) taken for main run of multiPIM on the original data.}
  \item{g.time}{time in seconds taken for running g models.}
  \item{Q.time}{time in seconds taken for running Q models.}
  \item{g.sl.time}{if g.method is "sl", time in seconds taken for running cross-validation of g models.}
  \item{Q.sl.time}{if Q.method is "sl", time in seconds taken for running cross-validation of Q models.}
  \item{g.sl.cand.times}{if g.method is "sl", named vector containing time taken, with each element corresponding to a super learner candidate for g.}
  \item{Q.sl.cand.times}{if Q.method is "sl", named vector containing time taken, with each element corresponding to a super learner candidate for Q.}

Note that all timing results apply only to the first run of code{link{multiPIM}} on the original data, not the subsequent bootstrap runs.
}

\author{Stephan Ritter, with design contributions from Alan Hubbard and Nicholas Jewell.}

\references{
Ritter, Stephan J.,  Jewell, Nicholas P. and Hubbard, Alan E. (2011) \dQuote{Variable Importance Analysis with the multiPIM R Package.} \emph{U.C. Berkeley Division of Biostatistics Working Paper Series}, Working Paper 286.
\url{http://www.bepress.com/ucbbiostat/paper286}

Hubbard, Alan E. and van der Laan, Mark J. (2008) \dQuote{Population Intervention Models in Causal Inference.} \emph{Biometrika} \bold{95}, 1: 35--47.

Young, Jessica G., Hubbard, Alan E., Eskenazi, Brenda, and Jewell, Nicholas P. (2009) \dQuote{A Machine-Learning Algorithm for Estimating and Ranking the Impact of Environmental Risk Factors in Exploratory Epidemiological Studies.} \emph{U.C. Berkeley Division of Biostatistics Working Paper Series}, Working Paper 250.
\url{http://www.bepress.com/ucbbiostat/paper250}

van der Laan, Mark J. and Rose, Sherri (2011) \emph{Targeted Learning}, Springer, New York.
ISBN: 978-1441997814

Sinisi, Sandra E., Polley, Eric C., Petersen, Maya L, Rhee, Soo-Yon and van der Laan, Mark J. (2007) \dQuote{Super learning: An Application to the Prediction of HIV-1 Drug Resistance.} \emph{Statistical Applications in Genetics and Molecular Biology} \bold{6}, 1: article 7.
\url{http://www.bepress.com/sagmb/vol6/iss1/art7}

van der Laan, Mark J., Polley, Eric C. and Hubbard, Alan E. (2007) \dQuote{Super learner.} \emph{Statistical applications in genetics and molecular biology} \bold{6}, 1: article 25.
\url{http://www.bepress.com/sagmb/vol6/iss1/art25}
}

\seealso{
\code{\link{multiPIM}} for the main function which is called by \code{multiPIMboot}.

\code{\link{summary.multiPIM}} for printing summaries of the results.

\code{\link{Candidates}} to see which candidates are currently available, and for information on writing user-defined super learner candidates and regression methods.
}

\examples{
## Warning: This would take a very long time to run!
\dontrun{
## load example from multiPIM help file

example(multiPIM)

## this would run 5000 bootstrap replicates:

boot.result <- multiPIMboot(Y, A)

summary(boot.result)}}
