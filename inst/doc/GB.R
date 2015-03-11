## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE-----------------------------------
library(RePsychLing)
library(lme4)
library(knitr)
options(width=92,show.signif.stars = FALSE)

## ----gbstr--------------------------------------------------------------------------------
str(gb12)
summary(gb12)

## ----m0-----------------------------------------------------------------------------------
print(summary(m0 <- lmer(
    sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF+(1+T+F+TF|session)+(1+T+P+TP|item), gb12,REML=FALSE)), corr=FALSE)

## ----m0PCA--------------------------------------------------------------------------------
summary(rePCA(m0))

## ----gbm02--------------------------------------------------------------------------------
m1 <- 
  lmer(sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF + (1+T+F+TF||session) + (1+T+P+TP||item),
       gb12, REML=FALSE)
VarCorr(m1)
anova(m1, m0)

## ----m2-----------------------------------------------------------------------------------
m2 <-
  lmer(sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF + (1+T+F||session) + (1+T||item), 
       gb12, REML=FALSE)
VarCorr(m2)
anova(m2, m1, m0)

## ----m3-----------------------------------------------------------------------------------
m3 <-
  lmer(sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF + (1+T+F||session) + (1|item),
       gb12, REML=FALSE)
VarCorr(m3)
anova(m3, m2)

## ----m4-----------------------------------------------------------------------------------
m4 <- 
  lmer(sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF + (1+T+F|session) + (1|item), 
       gb12, REML=FALSE)
VarCorr(m4)
anova(m3, m4)

## ----m5-----------------------------------------------------------------------------------
m5 <-
  lmer(sottrunc2 ~ 1+T+P+F+TP+TF+PF+TPF + (1+F|session) + (0+T|session) + (1|item),
       gb12, REML=FALSE)
VarCorr(m5)
anova(m5, m4)

## ----versions-----------------------------------------------------------------------------
sessionInfo()

