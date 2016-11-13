\name{vietnamese}
\alias{vietnamese}
\docType{data}
\title{vietnamese visual lexical decision from Pham and Baayen (2015)}
\description{
  Data described in Pham and Baayen (2015)
}
\format{
  A data frame with 15021 observations on the following 6 variables.
  \describe{
    \item{\code{HeadWord}}{a factor specifying the word stimuli}
    \item{\code{RTinv}}{a numeric vector with transformed reaction times (-1000/RT)}
    \item{\code{MidLevelTone}}{a factor specifying whether the first syllabeme carries mid level tone (TRUE/FALSE)}
    \item{\code{LogFreq}}{the frequency of the compound}
    \item{\code{LogFreqSyl1}}{the frequency of the left syllabeme}
    \item{\code{LogFreySyl2}}{the frequency of the right syllabeme}
  }
}
\source{
  Pham, H. and Baayen, R. H. (2015) Vietnamese compounds show an anti-frequency effect in visual lexical decision.
  \emph{Language, Cognition, and Neuroscience}, 
  \bold{30.9}, pp. 1077-1095.
}
\references{
  Pham, H. and Baayen, R. H. (2015) Vietnamese compounds show an anti-frequency effect in visual lexical decision.
  \emph{Language, Cognition, and Neuroscience}, 
  \bold{30.9}, pp. 1077-1095.
}
\examples{
data(vietnamese)
library(mgcv)
vietnamese.gam = bam(RTinv ~ MidLevelTone + s(LogFreq) + te(LogFreqSyl1, LogFreqSyl2), data=vietnamese, method="ML")
summary(vietnamese.gam)
plot(vietnamese.gam, pages=1)
}
\keyword{datasets}
