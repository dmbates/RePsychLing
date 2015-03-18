## ----preliminaries,include=FALSE,cache=FALSE-----------------------------
library(lme4)
library(RePsychLing)
library(knitr)
opts_chunk$set(cache=FALSE,comment=NA)
options(show.signif.stars=FALSE,width=92,digits = 5)

## ----m0------------------------------------------------------------------
summary(m0 <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2+c3|subj), KWDYZ, REML=FALSE))

## ----theta---------------------------------------------------------------
zapsmall(getME(m0,"theta"))

## ----chf-----------------------------------------------------------------
zapsmall(chf <- getME(m0,"Tlist")[[1]])

## ----relvcov-------------------------------------------------------------
tcrossprod(chf)

## ----vcov----------------------------------------------------------------
tcrossprod(getME(m0,"sigma")*chf)

## ----rowlengths----------------------------------------------------------
(rowlengths <- sqrt(rowSums(chf*chf)))

## ----corr----------------------------------------------------------------
tcrossprod(chf/rowlengths)

## ----chfsvd--------------------------------------------------------------
(svd0 <- svd(chf,nv=0))

## ----orthogonal----------------------------------------------------------
zapsmall(crossprod(svd0$u))

## ----stddevs-------------------------------------------------------------
zapsmall(getME(m0,"sigma")*svd0$d)

## ----propvar-------------------------------------------------------------
vc <- svd0$d^2   # variances of principal components
zapsmall(vc/sum(vc))

## ----cumsum--------------------------------------------------------------
cumsum(vc/sum(vc))

## ----rePCAm0-------------------------------------------------------------
prc <- rePCA(m0)
class(prc)
length(prc)
names(prc)
class(prc$subj)
prc$subj
summary(prc$subj)

## ----m1------------------------------------------------------------------
VarCorr(m1 <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2+c3||subj), KWDYZ, REML=FALSE))

## ----m1Tlist-------------------------------------------------------------
getME(m1,"Tlist")

## ----bdiag---------------------------------------------------------------
bdiag(getME(m1,"Tlist"))

## ----svdbdiag------------------------------------------------------------
svd(bdiag(getME(m1,"Tlist")),nv=0)

## ----m1rePCA-------------------------------------------------------------
summary(rePCA(m1))

## ----m2------------------------------------------------------------------
VarCorr(m2 <- lmer(rt ~ 1+c1+c2+c3 + (1+c1+c2|subj) + (0+c3|subj), KWDYZ, REML=FALSE))

## ----m2chf---------------------------------------------------------------
(chf <- bdiag(getME(m2,"Tlist")))
svd(chf,nu=0,nv=0)
summary(rePCA(m2))

## ----m3------------------------------------------------------------------
m3 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1|subj) + (1+P|item), kb07, REML=FALSE)
print(summary(m3),corr=FALSE)

## ----m3rePCA-------------------------------------------------------------
summary(rePCA(m3))

## ----versions------------------------------------------------------------
sessionInfo()


