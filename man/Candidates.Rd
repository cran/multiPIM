\name{Candidates}
\alias{Candidates}
\alias{all.bin.cands}
\alias{all.cont.cands}
\alias{default.bin.cands}
\alias{default.cont.cands}
\title{Super learner candidates (regression methods) available for use with the multiPIM and multiPIMboot functions}
\description{When the \bold{multiPIM} package is loaded, four character vectors are made available to the user. They are defined as follows:

\preformatted{all.bin.cands <- c("polyclass", "penalized.bin", "main.terms.logistic",
                   "rpart.bin")}
\preformatted{all.cont.cands <- c("polymars", "lars", "main.terms.linear", "penalized.cont",
                    "rpart.cont")}
\code{default.bin.cands <- c("polyclass", "penalized.bin", "main.terms.logistic")}

\code{default.cont.cands <- c("polymars", "lars", "main.terms.linear")}

These vectors (or subsets thereof) can be supplied as arguments to the \code{multiPIM} and the \code{multiPIMboot} functions, in order to specify which regression methods should be used to estimate the nuisance parameters g(0, W) and Q(0|W). The user may also supply custom written regression methods or super learner candidates. The mechanism for this is described in a section below.}

  \section{Candidates}{
    \subsection{polyclass and polymars}{

These candidate use the functions \code{\link[polspline]{polyclass}} and \code{\link[polspline]{polymars}} from the package \bold{polspline}.

    }
    \subsection{penalized.bin and penalized.cont}{

These candidates perform L1 penalized logistic (penalized.bin) or linear (penalized.cont) regression using the function \code{\link[penalized]{penalized}} from the package \bold{penalized}. The  value of the L1 penalty is selected by cross validation (using the \code{\link[penalized]{profL1}} function).

    }
    \subsection{lars}{

This candidate uses the function \code{\link[lars]{lars}} from the package \bold{lars}. Cross validation is performed using the function \code{\link[lars]{cv.lars}.}

    }
    \subsection{main.terms.logistic and main.terms.linear}{

These candidates perform standard main terms logistic or linear regression, using the functions \code{\link[stats]{glm}} and \code{\link[stats]{lm}}.

    }
%     \subsection{randomForest.bin and randomForest.cont}{

% These candidates use the function \code{\link[randomForest]{randomForest}} from the package \bold{randomForest}. Note that these candidates have not been included among the default candidates due to several issues which may lead to poor performance of the PIM estimator, such as e.g. the possibility of overfitting.

%     }
    \subsection{rpart.bin and rpart.cont}{

These candidates use the function \code{\link[rpart]{rpart}} from the package \bold{rpart}. They are not included as default candidates since their performance is very poor when predictor variables are continuous. It is not recommended to use these candidates if \code{W} contains any continuous variables.

    }
}

\section{Forcing of Variables into Q Models}{
Since some of the available candidates (such as polyclass/polymars, penalized) will sometimes completely drop an input variable from the model, it is necessary to have some mechanism to make sure that the relevant exposure variable stays in the model. How this is done for each candidate is described in greater detail in the technical report referenced below (Ritter, Jewell and Hubbard, 2011).
}

\section{User-Defined Regression Methods and Super Learner Candidates}{

Below is the code which defines the \code{main.terms.logistic} candidate function. This is an example of the form that functions which are passed as elements of the \code{extra.cands} argument should have. (See code for the multiPIM function in the file multiPIM/R/functions.R to see other examples of how candidates are defined.)

\preformatted{candidate.functions$main.terms.logistic <- function(X, Y, newX, force,
                                                    standardize) \{

  result <- vector("list", length = 4)
  names(result) <- c("preds", "fitted", "main.model", "force.model")
  class(result) <- "main.terms.logistic.result"

  formula.string <- paste(names(Y), "~", paste(names(X), collapse = "+"),
                          sep = "")

  result$main.model <- glm(formula.string, data = cbind(Y, X),
                           family = binomial, model = FALSE, y = FALSE)

  result$preds <- predict(result$main.model, newdata = newX,
                          type = "response")
  result$fitted <- predict(result$main.model, type = "response")
  return(result)
\}
}


The functions muse take these four arguments: \code{X} will be a data frame of predictors, \code{Y} will be a single-column data frame containing the outcome, and \code{newX} will be a data frame with columns corresponding to the columns of \code{X}, containing data on which to predict the outcome based on the model fit by the function. \code{force} will be an integer specifying the column number (of \code{X} and \code{newX}) corresponding to the variable which should be forced into the model in case the function is being used to fit a Q model, and standardize will just be a logical value indicating whether or not to standardize the input variables before running the fitting algorithm (this is meant for algorithms like penalized, where the scale of the predictors will make a difference).

For g models, the force argument will be missing when the function is called, so if the function must do something differently in order to force in a variable (unlike the \code{main.terms.logistic} function above), then one can use a conditional such as:

\code{if(missing(force)) ...}

in order to differentiate between g and Q models.

The list returned by the functions must have a slot named \code{preds} containing the predictions (or predicted probabilities for the binary outcome case), based on \code{newX}. Also, for the TMLE estimator it is necessary to have the fitted values, i.e. the predictions on \code{X}. This should be returned in the \code{fitted} slot. Note that for binary outcomes, these predictions and fitted values should be predicted probabilities that the outcome is equal to 1. Thus, if the candidate is being used for a g model, where the outcome is an exposure variable, the returned values will be estimated probabilities that the exposure variable is equal to 1. The probabilities will be converted as necessary elsewhere in the multiPIM function.

The other slots (main.model and force.model), and setting the class of the object returned, are not necessary for the multiPIM function to work correctly, but may be useful if one would like to inspect the final g and Q models after running the function (see the \code{return.final.models} argument).
}

\author{Stephan Ritter, with design contributions from Alan Hubbard and Nicholas Jewell.}

\references{

multiPIM:

Ritter, Stephan J.,  Jewell, Nicholas P. and Hubbard, Alan E. (2011) \dQuote{Variable Importance Analysis with the multiPIM R Package.} \emph{U.C. Berkeley Division of Biostatistics Working Paper Series}, Working Paper 286.
\url{http://www.bepress.com/ucbbiostat/paper286}

lars:

Efron, B et al. (2004). \dQuote{Least angle regression}.
\emph{The Annals of statistics}, \bold{32}(2):407-499.

penalized:

Goeman, J. J. (2010). \dQuote{L1 penalized estimation in the cox proportional hazards model}.
\emph{Biometrical Journal}, \bold{52}(1):70-84.

polyclass and polymars:

Friedman, J. H. (1991). \dQuote{Multivariate adaptive regression splines (with discussion)}. \emph{The Annals of Statistics}, \bold{19}:1-141.

Kooperberg, C. et al. (1997). \dQuote{Polychotomous regression}.
\emph{Journal of the American Statistical Association}, \bold{92}(437):117-127.

Stone, C. J. et al. (1997). \dQuote{The use of polynomial splines and their tensor products in extended linear modeling (with discussion)} .
\emph{Annals of Statistics}, \bold{25}:1371-1470.

% randomForest:

% Liaw, A. and Wiener, M. (2002). \dQuote{Classification and Regression by randomForest}.
% \emph{R News}, \bold{2}(3):18-22.

rpart:

Breiman, L. et al. (1984). \emph{Classification and regression trees}. Wadsworth International Group, Belmont, CA.
ISBN: 0534980538.
}

\seealso{
\code{\link{multiPIM}}, \code{\link{multiPIMboot}}
}
