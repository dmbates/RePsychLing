% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/createParamMx.R
\name{createParamMx}
\alias{createParamMx}
\title{Create Population Parameter Matrix}
\usage{
createParamMx(nexp = 100000L, simparam.env = getParamRanges(),
  firstseed = NULL, h0 = TRUE, outfile = NULL)
}
\arguments{
\item{nexp}{number of parameter-value sets to generate (default 100000)}

\item{simparam.env}{an environment containing ranges of population parameters}

\item{firstseed}{initial seed}

\item{h0}{logical value indicating if h0 is true or false}

\item{outfile}{name of save file for parameter-value matrix}
}
\value{
a matrix of generated population parameter values
}
\description{
Create a matrix of parameter values for simulations
}
\details{
Generates a matrix of parameters.  Each row are the population parameters
used to generate data from a single hypothetical "experiment."
}

