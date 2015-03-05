## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE-----------------------------------
library(lme4)
library(RePsychLing)
library(knitr)
#opts_chunk$set(cache=TRUE)
options(width=92,show.signif.stars = FALSE)

## ----kb07str------------------------------------------------------------------------------
str(kb07)

## ----m0, eval=FALSE-----------------------------------------------------------------------
#  m0 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC+(1+S+P+C+SP+SC+PC+SPC|subj)+(1+S+P+C+SP+SC+PC+SPC|item),kb07,start=kb07maxLMMtheta,control=lmerControl(optimizer="none"))

## ----load_m0,include=FALSE----------------------------------------------------------------
load("kb07_m0.rda")  # must be in working directory (or provide path)
m0@theta <- Re(kb07maxLMMtheta)

## ----showm0-------------------------------------------------------------------------------
summary(m0)

## ----chf0---------------------------------------------------------------------------------
summary(rePCA(m0))

## ----m1-----------------------------------------------------------------------------------
m1 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1+S+P+C+SP+SC+PC+SPC||subj) +
             (1+S+P+C+SP+SC+PC+SPC||item), kb07, REML=FALSE)
print(summary(m1),corr=FALSE)

anova(m1, m0)  

## ----rePCAm1------------------------------------------------------------------------------
summary(rePCA(m1))

## ----m2-----------------------------------------------------------------------------------
m2 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1+S+P+C+SP+PC||subj) +
             (1+P+C+SC+PC+SPC||item), kb07, REML=FALSE)
anova(m2, m1)  # not significant: prefer m2 over m1

## ----m3-----------------------------------------------------------------------------------
m3 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1+C+SP+PC||subj) +  (1+P+C+PC||item),kb07,REML=FALSE)
anova(m3,  m2)  # not significant: prefer m3 over m2

## ----m4-----------------------------------------------------------------------------------
m4 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1|subj) + (1|item) + (0+P|item),kb07,REML=FALSE)
anova(m4, m3)  # not significant: prefer m4 over m3
anova(m4, m1)  # not significant: prefer m4 over m1 (no accumulation)

## ----m5-----------------------------------------------------------------------------------
m5 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1|subj) + (1|item), data=kb07, REML=FALSE)
anova(m5, m4)  # significant: prefer m4 over m5

## ----m6-----------------------------------------------------------------------------------
m6 <- lmer(RTtrunc ~ 1+S+P+C+SP+SC+PC+SPC + (1|subj) + (1+P|item), kb07, REML=FALSE)
print(summary(m6), corr=FALSE)

anova(m4, m6)  # significant: prefer m6 over m4
anova(m6, m0)  # not significant: prefer m6 over m0 (no accumulation)

## ----versions-----------------------------------------------------------------------------
sessionInfo()

