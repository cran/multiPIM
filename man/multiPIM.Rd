\name{multiPIM}
\alias{multiPIM}
\title{Estimate Variable Importances for Multiple Exposures and Outcomes}
\description{The parameter of interest is a type of causal attributable risk.  One effect measure (and a corresponding plug-in standard error) will be calculated for each exposure-outcome pair. The default is to use a Targeted Maximum Likelihood Estimator (TMLE). The other available estimators are Inverse Probability of Censoring Weighted (IPCW), Double-Robust IPCW (DR-IPCW), and Graphical Computation (G-COMP) estimators. PIM stands for Population Intervention Model.}
\usage{multiPIM(Y, A, W = NULL,
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
         check.input = TRUE,
         verbose = FALSE,
         extra.cands = NULL,
         standardize = TRUE,
         \dots)}
\arguments{
  \item{Y}{a data frame of outcomes containing only numeric (integer or double) values. See details for the default method of determining, based on the values in \code{Y}, which regression types to allow for modelling Q. Must have unique names.}
  \item{A}{a data frame containing binary exposure variables. \emph{Binary} means that all values must be either 0 (indicating unexposed, or part of target group) or 1 (indicating exposed or not part of target group). Must have unique names.}
  \item{W}{an optional data frame containing possible confounders of the effects of the variables in \code{A} on the variables in \code{Y}. No effect measures will be calculated for these variables. May contain numeric (integer or double), or factor values. Must be left as \code{NULL} if not required. See details.}
  \item{estimator}{the estimator to be used. The default is \code{"TMLE"}, for the targeted maximum likelihood estimator. Alternatively, one may specify \code{"DR-IPCW"}, for the Double-Robust Inverse Probability of Censoring-Weighted estimator, or \code{"IPCW"}, for the regular IPCW estimator, or \code{"G-COMP"} for the Graphical Computation estimator. If the regular IPCW estimator is selected, all arguments which begin with the letter Q are ignored, since only g (the regression of each exposure on possible confounders) needs to be modeled in this case. Similarly, if the G-COMP estimator is selected, all arguments which begin with the letter g, as well as the truncate argument, will be ignored, since only Q needs to be modeled in this case. Note: an additional characteristic of the G-COMP estimator is that there are no plug-in standard errors available. If you want to use G-COMP and you need standard errors, the \code{multiPIMboot} function is available and will provide bootstrap standard errors.}
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
  \item{return.final.models}{single logical value indicating whether final g and Q models should be returned by the function (in the slots \code{g.final.models} and \code{Q.final.models}). Default is \code{TRUE}. If memory is a concern, you will probably want to set this to FALSE.}
  \item{na.action}{currently ignored. If any of \code{Y}, \code{A} or (a non-null) \code{W} has missing values, \code{multiPIM} will throw an error.}
  \item{check.input}{a single logical value indicating whether all of the input to the function should be subjected to strict error checking. \code{FALSE} is not recommended.}
  \item{verbose}{a single logical value indicating whether messages about the progress of the evaluation should be printed out. Some of the candidate algorithms may print messages even when \code{verbose} is set to \code{FALSE}.}
  \item{extra.cands}{a named list of functions. This argument provides a way for the user to specify his or her own functions to use either as stand-alone regression methods, or as candidates for a super learner. See details.}
  \item{standardize}{should all predictor variables be standardized before certain regression methods are run. Passed to all candidates, but only used by some (at this point, lars, penalized.bin and penalized.cont)}
  \item{\dots}{currently ignored.}
}

\details{
The parameter of interest is a type of attributable risk.  This means
that it is a measure (adjusted for known confounders) of the
\emph{difference} between the mean value of \code{Y} for the units in
the target (or unexposed) group and the \emph{overall} mean value of
\code{Y}. Units which are in the target (or unexposed) group with
respect to one of the variables in \code{A} are characterized as such by
having the value 0 in the respective column of \code{A}. Members of the
the non-target (or exposed) group should have a 1 in that column of
\code{A}. \emph{Assuming all causal assumptions hold} (see the paper), each
parameter estimate can be thought of as estimating the hypothetical
effect on the respective outcome of totally eliminating the respective
exposure from the population (i.e. setting everyone to 0 for that
exposure). For example, in the case of a binary outcome, a parameter
estimate for exposure x and outcome y of -0.03 could be interpreted
as follows: the effect of an intervention in which the entire population
was set to exposure x = 0 would be to reduce the level of outcome y
by 3 percentage points.

If \code{check.input} is \code{TRUE} (which is the default and is highly
recommended), all of the arguments will be checked to make sure they
have permissible values. Many of the arguments, especially those for
which a single logical value (\code{TRUE} or \code{FALSE}) or a single
character value (such as, for example, \code{"all"}) is expected, are
checked using the \code{identical} function, which means that if any of
these arguments has any extraneous attributes (such as names), this may
cause \code{multiPIM} to throw an error.

On the other hand, the arguments \code{Y} and \code{A} (and \code{W} if
it is non-null) \emph{must} have valid names attributes. \code{multiPIM}
will throw an error if there is any overlap between the names of the
columns of these data frames, or if any of the names cannot be used in a
\code{formula} (for example, because it begins with a number and not a
letter).

By default, the regression methods which will be allowed
for fitting models for Q will be determined from the contents of
\code{Y} as follows: if all values in Y are either 0 or 1 (i.e. all
outcomes are binary), then \dQuote{logistic}-type regression methods
will be used (and only these methods will be allowed in the arguments
\code{Q.method} and \code{Q.sl.cands}); however, if there are any values
in \code{Y} which are not equal to 0 or 1 then it will be assumed that
all outcomes are continuous, \dQuote{linear}-type regression will be
used, and the values allowed for \code{Q.method} and \code{Q.sl.cands}
will change accordingly.  This behavior can be overriden by specifying
\code{Q.type} as either \code{"binary.outcome"} (for logistic-type
regression), or as \code{"continuous.outcome"} (for linear-type
regression). If \code{Q.type} is specified, \code{Y} will not be checked for
binaryness.

The values allowed for \code{Q.method} (which should have length 1) are:
either \code{"sl"} if one would like to use super learning, or one of the
elements of the vector \code{\link{all.bin.cands}} (for the binary outcome case),
or of \code{\link{all.cont.cands}} (for the continuous outcome case), if one would
like to use only a
particular regression method for all modelling of Q. If \code{Q.method}
is given as \code{"sl"}, then the candidates used by the super learner
will be determined from the value of \code{Q.sl.cands}. If the value of
\code{Q.sl.cands} is \code{"default"}, then the candidates listed in either
\code{\link{default.bin.cands}} or \code{\link{default.cont.cands}} will
be used. If the value of \code{Q.sl.cands} is \code{"all"}, then the candidates
listed in either \code{\link{all.bin.cands}} or \code{\link{all.cont.cands}}
will be used. The function will automatically choose the candidates which
correspond to the correct outcome type (binary or continuous). Alternatively,
one may specify \code{Q.sl.cands} explicitly as a vector of names of the
candidates to be used.

If \code{A} has more than one column, the \code{adjust.for.other.As}
argument can be used to specify whether the other
columns of \code{A} should possibly be included in the g and Q models
which will be used in calculating the effect of a
certain column of \code{A} on each column of \code{Y}.

With the argument \code{extra.cands}, one may supply alternative R
functions to be used as stand-alone regression methods, or as super
learner candidates, within the \code{multiPIM}
function. \code{extra.cands} should be given as a named list of
functions. See \link{Candidates} for the form (e.g. arguments) that the
functions in this list should have. In order to supply your own stand
alone regression method for g or Q, simply specify \code{g.method} or
\code{Q.method} as the name of the function you want to use (i.e. the
corresponding element of the names attribute of \code{extra.cands}). To
add candidates to a super learner, simply use the corresponding names of
your functions (from the names attribute of \code{extra.cands}) when you
supply the \code{g.sl.cands} or \code{Q.sl.cands} arguments. Note that
you may mix and match between your own extra candidates and the built-in
candidates given in the \code{\link{all.bin.cands}} and
\code{\link{all.cont.cands}} vectors. Note
also that extra candidates must be explicitly specified as
\code{g.method}, \code{Q.method}, or as elements of \code{g.sl.cands} or
\code{Q.sl.cands} -- Specifying \code{Q.sl.cands} as \code{"all"} will not
cause any extra candidates to be used.
}

\value{
Returns an object of class \code{"multiPIM"} with the following elements:
  \item{param.estimates}{a matrix of dimensions \code{ncol(A)} by \code{ncol(Y)} with \code{rownames} equal to \code{names(A)} and \code{colnames} equal to \code{names(Y)}, with each element being the estimated causal attributable risk for the exposure given by its row name vs. the outcome given by its column name.}
  \item{plug.in.stand.errs}{a matrix with the same dimensions as \code{param.estimates} containing the corresponding plug-in standard errors of the parameter estimates. These are obtained from the influence curve. Note: plug-in standard errors are not available for \code{estimator = "G-COMP"}. This field will be set to \code{NA} in this case.}
  \item{call}{a copy of the call to \code{multiPIM} which generated this object.}
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
  \item{boot.param.array}{this slot will be \code{NULL} for objects returned by the \code{multiPIM} function. See \code{\link{multiPIMboot}} for details on what this slot is actually used for.}
  \item{main.time}{total time (in seconds) taken to generate this multiPIM result.}
  \item{g.time}{time in seconds taken for running g models.}
  \item{Q.time}{time in seconds taken for running Q models.}
  \item{g.sl.time}{if g.method is "sl", time in seconds taken for running cross-validation of g models.}
  \item{Q.sl.time}{if Q.method is "sl", time in seconds taken for running cross-validation of Q models.}
  \item{g.sl.cand.times}{if g.method is "sl", named vector containing time taken, with each element corresponding to a super learner candidate for g.}
  \item{Q.sl.cand.times}{if Q.method is "sl", named vector containing time taken, with each element corresponding to a super learner candidate for Q.}
}

\author{Stephan Ritter, with design contributions from Alan Hubbard and Nicholas Jewell.}

\references{
Ritter, Stephan J.,  Jewell, Nicholas P. and Hubbard, Alan E. (2014) \dQuote{R Package multiPIM: A Causal Inference Approach to Variable Importance Analysis} \emph{Journal of Statistical Software} \bold{57}, 8: 1--29. 
\url{http://www.jstatsoft.org/v57/i08/.}

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
\code{\link{multiPIMboot}} for running \code{multiPIM} with automatic bootstrapping to get standard errors.

\code{\link{summary.multiPIM}} for printing summaries of the results.

\code{\link{Candidates}} to see which candidates are currently available, and for information on writing user-defined super learner candidates and regression methods.
}

\examples{
num.columns <- 3
num.obs <- 250

set.seed(23)

## use rbinom with size = 1 to make a data frame of binary data

A <- as.data.frame(matrix(rbinom(num.columns*num.obs, 1, .5),
                          nrow = num.obs, ncol = num.columns))

## let Y[,i] depend only on A[,i] plus some noise
## (start with the noise then add a multiple of A[,i] to Y[,i])

Y <- as.data.frame(matrix(rnorm(num.columns*num.obs),
                          nrow = num.obs, ncol = num.columns))
for(i in 1:num.columns)
  Y[,i] <- Y[,i] + i * A[,i]

## make sure the names are unique

names(A) <- paste("A", 1:num.columns, sep = "")
names(Y) <- paste("Y", 1:num.columns, sep = "")

result <- multiPIM(Y, A)

summary(result)}
