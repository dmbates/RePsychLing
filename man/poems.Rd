\name{poems}
\alias{poems}
\docType{data}
\title{Baayen and Milin }
\description{
  Data described in Baayen and Milin (2010). 
}
\format{
  A data frame with 275996 observations on the following 24 variables.
  \describe{
    \item{\code{ReadingTime}}{a numeric vector of self-paced reading times}
    \item{\code{Subject}}{a factor with participant identifiers}
    \item{\code{Sex}}{a factor with levels \code{m} (male) and \code{f} (female)}
    \item{\code{Age}}{a numeric vector specifying the participant's age}
    \item{\code{NPoems}}{a numeric vector of the self-reported maximum number of poems read annually, according to a four-choice question}
    \item{\code{MultipleChoiceRT}}{a numeric vector with the response latency to the four-choice question}
    \item{\code{Trial}}{a numeric vector specifying the rank of the item in the subject's experimental list}
    \item{\code{NumberOfWordsIntoLine}}{a numeric vector specifying the position of the
      item in the line of poetry being read}
    \item{\code{PositionBegMidEnd}}{a factor specifying whether the word was initial
      \code{beg}, medial \code{mid} or final \code{end} in the sentence}
    \item{\code{SentenceLength}}{a numeric vector specifying sentence length}
    \item{\code{Poem}}{a factor with as levels identifiers for the poems}
    \item{\code{Word}}{a factor with as levels identifiers for the words}
    \item{\code{WordFrequencyInPoem}}{a numeric vector specifying the frequency of the word in the poem}
    \item{\code{RhymeFreqInPoem}}{a numeric vector specifying the frequency of the word's rhyme in the poem}
    \item{\code{OnsetFreqInPoem}}{a numeric vector specifying the frequency of the word's onset in the poem}
    \item{\code{WordLength}}{a numeric vector specifying the length of the word in letters}
    \item{\code{FamilySize}}{a numeric vector specifying the count of morphological family members}
    \item{\code{InflectionalEntropy}}{a numeric vector specifying Shannon's entropy calculated over the probability distribution of a word's inflected variants}
    \item{\code{LemmaFrequency}}{a numeric vector specifying the frequency of occurrence of the word in the lemma subsection of the CELEX lexical database}
    \item{\code{WordFormFrequency}}{a numeric vector specifying the frequency of occurrence of the word's inflected form in the word form subsection of the CELEX lexical database}
    \item{\code{NumberOfMeanings}}{a numeric vector specifying the number of synsets in WordNet in which the word is listed}
    \item{\code{IsFunctionWord}}{a factor specifying whether the word is a function word \code{TRUE} or not \code{FALSE}}
    \item{\code{HasPunctuationMark}}{a factor specifying whether the word is followed by a punctuation mark, levels \cite{FALSE} (absent) and \cite{TRUE} (present)}
    \item{\code{NumberOfMorphemes}}{a numeric vector specifying the scaled number of morphemes in a word}
  }
}
\source{
  Baayen, R. H. and Milin, P (2010) Analyzing reaction times. 
  \emph{International Journal of Psychological Research}, 
  \bold{3.2}, pp. 12-28.
}
\references{
  Baayen, R. H. and Milin, P (2010) Analyzing reaction times. 
  \emph{International Journal of Psychological Research}, 
  \bold{3.2}, pp. 12-28.
}
\examples{
data(poems)
par(mfrow=c(2,4))
qqnorm(poems$ReadingTime)
qqnorm(poems$WordFormFrequency)
qqnorm(poems$LemmaFrequency)
qqnorm(poems$FamilySize)
qqnorm(poems$MultipleChoiceRT)
qqnorm(poems$NPoems)
qqnorm(poems$NumberOfMeanings)
poems$LogReadingTime        = log(poems$ReadingTime)
poems$LogWordFormFrequency  = log(poems$WordFormFrequency+1)
poems$LogLemmaFrequency     = log(poems$LemmaFrequency+1)
poems$RecFamilySize         = -100/(poems$FamilySize+1)
poems$LogMultipleChoiceRT   = log(poems$MultipleChoiceRT)
poems$LogNPoems             = log(poems$NPoems)
poems$LogNumberOfMeanings   = log(poems$NumberOfMeanings+1)

\dontrun{

p = poems[,c("Age", "LogNPoems", "LogMultipleChoiceRT", "NumberOfWordsIntoLine", "SentenceLength",
                     "WordFrequencyInPoem", "RhymeFreqInPoem", "OnsetFreqInPoem", "WordLength", 
                     "NumberOfMorphemes",
                     "RecFamilySize", "InflectionalEntropy", "LogLemmaFrequency", "LogWordFormFrequency",
                     "LogNumberOfMeanings")]
pc = prcomp(p,center=TRUE, scale=TRUE)
round(pc$rotation[,1:7],2)
#                        PC1   PC2   PC3   PC4   PC5   PC6   PC7
#Age                    0.00  0.01  0.00  0.03  0.61  0.49 -0.01
#LogNPoems              0.00 -0.01  0.01 -0.01 -0.70 -0.02  0.00
#LogMultipleChoiceRT    0.00  0.00  0.00  0.01 -0.37  0.87 -0.02
#NumberOfWordsIntoLine  0.03 -0.19 -0.39 -0.56  0.01  0.02 -0.05
#SentenceLength        -0.09 -0.20 -0.40 -0.52  0.01  0.01 -0.11
#WordFrequencyInPoem   -0.30 -0.36  0.14  0.11  0.00 -0.01 -0.06
#RhymeFreqInPoem       -0.24 -0.54  0.15  0.07  0.01  0.00  0.11
#OnsetFreqInPoem       -0.20 -0.56  0.14  0.06  0.01  0.00  0.13
#WordLength             0.41 -0.16  0.18 -0.08  0.00  0.00  0.15
#NumberOfMorphemes      0.17 -0.13  0.24 -0.03  0.01 -0.01 -0.83
#RecFamilySize         -0.35  0.20 -0.02 -0.11  0.00  0.01  0.34
#InflectionalEntropy    0.30 -0.19 -0.42  0.36 -0.01 -0.01 -0.02
#LogLemmaFrequency     -0.43  0.13 -0.21  0.18 -0.01 -0.01 -0.27
#LogWordFormFrequency  -0.45  0.16 -0.12  0.10  0.00 -0.01 -0.25
#LogNumberOfMeanings    0.11 -0.15 -0.55  0.44 -0.01 -0.01  0.01


poems$PC1 = pc$x[,1]
poems$PC2 = pc$x[,2]
poems$PC3 = pc$x[,3]
poems$PC4 = pc$x[,4]
poems$PC5 = pc$x[,5]
poems$PC6 = pc$x[,6]
poems$PC7 = pc$x[,7]

library(lme4)
poems.lmer = lmer(LogReadingTime ~ 
  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
  HasPunctuationMark*Sex + Trial + PositionBegMidEnd +
  (1|Poem) + (1|Word) + (1|Subject),
  #(1+LogWordFormFrequency+NumberOfMorphemes|Subject) ,
  data=poems, REML=FALSE)
print(summary(poems.lmer), corr=FALSE)

chf <- diag(c(diag(
  getME(poems.lmer, "Tlist")[[2]]), 
  getME(poems.lmer, "Tlist")[[1]], 
  getME(poems.lmer, "Tlist")[[3]]))
chf[1:3, 1:3] <- getME(poems.lmer, "Tlist")[[2]]             

sv <- svd(chf)
round(sv$d^2/sum(sv$d^2)*100, 1)
}
}
\keyword{datasets}
