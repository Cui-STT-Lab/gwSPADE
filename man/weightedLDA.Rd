% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keyATM.R
\name{weightedLDA}
\alias{weightedLDA}
\title{Weighted LDA main function}
\usage{
weightedLDA(
  docs,
  Weight_Mat,
  model,
  number_of_topics,
  model_settings = list(),
  priors = list(),
  options = list(),
  keep = c()
)
}
\arguments{
\item{docs}{texts read via \code{\link[=keyATM_read]{keyATM_read()}}.}

\item{model}{Weighted LDA model: \code{base}, \code{covariates}, and \code{dynamic}.}

\item{number_of_topics}{the number of regular topics.}

\item{model_settings}{a list of model specific settings (details are in the online documentation).}

\item{priors}{a list of priors of parameters.}

\item{options}{a list of options (details are in the documentation of \code{\link[=keyATM]{keyATM()}}).}

\item{keep}{a vector of the names of elements you want to keep in output.}
}
\value{
A \code{keyATM_output} object containing:
\describe{
\item{V}{number of terms (number of unique words)}
\item{N}{number of documents}
\item{model}{the name of the model}
\item{theta}{topic proportions for each document (document-topic distribution)}
\item{phi}{topic specific word generation probabilities (topic-word distribution)}
\item{topic_counts}{number of tokens assigned to each topic}
\item{word_counts}{number of times each word type appears}
\item{doc_lens}{length of each document in tokens}
\item{vocab}{words in the vocabulary (a vector of unique words)}
\item{priors}{priors}
\item{options}{options}
\item{keywords_raw}{\code{NULL} for LDA models}
\item{model_fit}{perplexity and log-likelihood}
\item{pi}{estimated pi for the last iteration (\code{NULL} for LDA models)}
\item{values_iter}{values stored during iterations}
\item{number_of_topics}{number of topics}
\item{kept_values}{outputs you specified to store in \code{keep} option}
\item{information}{information about the fitting}
}
}
\description{
Fit weighted LDA models.
}
\examples{
\dontrun{
  library(keyATM)
  library(quanteda)
  data(keyATM_data_bills)
  bills_dfm <- keyATM_data_bills$doc_dfm  # quanteda dfm object
  keyATM_docs <- keyATM_read(bills_dfm)

  # Weighted LDA
  out <- weightedLDA(docs = keyATM_docs, model = "base",
                     number_of_topics = 5)

  # Visit our website for full examples: https://keyatm.github.io/keyATM/
}

}
\seealso{
\code{\link[=save.keyATM_output]{save.keyATM_output()}}, \url{https://keyatm.github.io/keyATM/articles/pkgdown_files/Options.html}
}
