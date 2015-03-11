## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE------------------
library(lme4)
library(knitr)
library(RePsychLing)
opts_chunk$set(comment=NA)
options(width=92,show.signif.stars=FALSE)

## ----strKWDYZ------------------------------------------------------------
str(KWDYZ)

## ----m0------------------------------------------------------------------
summary(m0 <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2+c3|subj), KWDYZ, REML=FALSE))
summary(rePCA(m0))

## ----chf0----------------------------------------------------------------
chf0 <- getME(m0,"Tlist")[[1]]
zapsmall(chf0)

## ----svd0----------------------------------------------------------------
sv0 <- svd(chf0)
sv0$d

## ----linalg--------------------------------------------------------------
(xx<-tcrossprod(chf0))
sum(diag(xx)) 
diag(xx)
str(sv0)
sv0$v
zapsmall(sv0$v)
sv0$u   # last column is the singular combination of random effects

## ----m1------------------------------------------------------------------
m1 <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2+c3||subj), KWDYZ, REML=FALSE)
VarCorr(m1)
summary(rePCA(m1))
anova(m1, m0)  # significant: too much of a reduction

## ----ci_p1, eval=FALSE---------------------------------------------------
## (m1.ci_profile <- confint(m1, method="profile"))

## ----ci_b1, eval=FALSE---------------------------------------------------
## (m1.ci_boot <- confint(m1, method="boot"))

## ----m2.2----------------------------------------------------------------
m2d <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2|subj), KWDYZ, REML=FALSE)
VarCorr(m2d)
summary(rePCA(m2d))

## ----m2.4----------------------------------------------------------------
print(summary(m2i <- lmer(lrt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | subj), 
                          REML=FALSE, data=KWDYZ)), corr=FALSE)  
summary(rePCA(m2i))
print(summary(m2j <- lmer(prt ~ 1 + c1 + c2 + c3 + (1 + c1 + c2 + c3 | subj),
                          REML=FALSE, data=KWDYZ)), corr=FALSE)  
summary(rePCA(m2j))

## ----versions------------------------------------------------------------
sessionInfo()

