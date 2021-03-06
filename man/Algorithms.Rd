% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AATtools.R
\name{Algorithms}
\alias{Algorithms}
\alias{aat_doublemeandiff}
\alias{aat_doublemediandiff}
\alias{aat_dscore}
\alias{aat_dscore_multiblock}
\alias{aat_regression}
\alias{aat_standardregression}
\alias{aat_doublemedianquotient}
\alias{aat_doublemeanquotient}
\alias{aat_singlemeandiff}
\alias{aat_singlemediandiff}
\title{AAT score computation algorithms}
\usage{
aat_doublemeandiff(ds, subjvar, pullvar, targetvar, rtvar, ...)

aat_doublemediandiff(ds, subjvar, pullvar, targetvar, rtvar, ...)

aat_dscore(ds, subjvar, pullvar, targetvar, rtvar, ...)

aat_dscore_multiblock(ds, subjvar, pullvar, targetvar, rtvar, blockvar, ...)

aat_regression(ds, subjvar, formula, aatterm, ...)

aat_standardregression(ds, subjvar, formula, aatterm, ...)

aat_doublemedianquotient(ds, subjvar, pullvar, targetvar, rtvar, ...)

aat_doublemeanquotient(ds, subjvar, pullvar, targetvar, rtvar, ...)

aat_singlemeandiff(ds, subjvar, pullvar, rtvar, ...)

aat_singlemediandiff(ds, subjvar, pullvar, rtvar, ...)
}
\arguments{
\item{ds}{A long-format data.frame}

\item{subjvar}{Column name of the participant identifier variable}

\item{pullvar}{Column name of the movement variable (0: avoid; 1: approach)}

\item{targetvar}{Column name of the stimulus category variable (0: control stimulus; 1: target stimulus)}

\item{rtvar}{Column name of the reaction time variable}

\item{...}{Other arguments passed on by functions (ignored)}

\item{blockvar}{name of the variable indicating block number}

\item{formula}{A regression formula to fit to the data to compute an AAT score}

\item{aatterm}{A character naming the formula term representing the approach bias.
Usually this is the interaction of the movement-direction and stimulus-category terms.}
}
\value{
A data.frame containing participant number and computed AAT score.
}
\description{
\itemize{
\item \code{aat_doublemeandiff} computes a mean-based double-difference score:

\code{(mean(push_target) - mean(pull_target)) - (mean(push_control) - mean(pull_control))}
\item \code{aat_doublemediandiff} computes a median-based double-difference score:

\code{(median(push_target) - median(pull_target)) - (median(push_control) - median(pull_control))}
\item \code{aat_dscore} computes D-scores for a 2-block design (see Greenwald, Nosek, and Banaji, 2003):

\code{((mean(push_target) - mean(pull_target)) - (mean(push_control) - mean(pull_control))) / sd(participant_reaction_times)}
\item \code{aat_dscore_multiblock} computes D-scores for pairs of sequential blocks
and averages the resulting score (see Greenwald, Nosek, and Banaji, 2003).
Requires extra \code{blockvar} argument, indicating the name of the block variable.
\item \code{aat_regression} and \code{aat_standardregression} fit regression models to participants' reaction times and extract a term that serves as AAT score.
\code{aat_regression} extracts the raw coefficient, equivalent to a mean difference score.
\code{aat_standardregression} extracts the t-score of the coefficient, standardized on the basis of the variability of the participant's reaction times.
These algorithms can be used to regress nuisance variables out of the data before computing AAT scores.
When using these functions, additional arguments must be provided:
\itemize{
\item \code{formula} - a formula to fit to the data
\item \code{aatterm} - the term within the formula that indicates the approach bias; this is usually the interaction of the pull and target terms.
}
\item \code{aat_doublemeanquotient} and \code{aat_doublemedianquotient} compute a log-transformed ratio of approach to avoidance for both stimulus categories and subtract these ratios:

\code{log(mean(pull_target) / mean(push_target)) - log(mean(pull_control) / mean(push_control))}
\item \code{aat_singlemeandiff} and \code{aat_singlemediandiff} subtract the mean or median approach reaction time from the mean or median avoidance reaction time.
These algorithms are only sensible if the supplied data contain a single stimulus category.
}
}
