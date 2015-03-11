## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE------------------
library(lme4)
library(RePsychLing)
library(knitr)
opts_chunk$set(comment=NA)
options(width=92,show.signif.stars = FALSE)

## ----bs10str-------------------------------------------------------------
str(bs10)

## ----m0,warning=FALSE----------------------------------------------------
m0 <- lmer(dif ~ 1+S+F+SF + (1+S+F+SF|subj) + (1+S+F+SF|item), bs10, REML=FALSE, start=thcvg$bs10$m0,
           control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1L),
                               check.conv.grad="ignore",check.conv.hess="ignore"))
summary(m0)

## ----sv_m0_s-------------------------------------------------------------
summary(rePCA(m0))

## ----m1------------------------------------------------------------------
print(summary(m1 <- lmer(dif ~ 1+S+F+SF + (1+S+F+SF||subj) + (1+S+F+SF||item), bs10, REML=FALSE)), corr=FALSE)
anova(m1, m0)

## ----pca1----------------------------------------------------------------
summary(rePCA(m1))

## ----m2------------------------------------------------------------------
print(summary(m2 <- lmer(dif ~ 1+S+F+SF +  (1+S+SF||subj), bs10,  REML=FALSE)), corr=FALSE)

## ----m21anova------------------------------------------------------------
anova(m2,m1)

## ----m3------------------------------------------------------------------
print(summary(m3 <- lmer(dif ~ 1+S+F+SF +  (1+S||subj), bs10, REML=FALSE)), corr=FALSE)
anova(m3, m2)

## ----m4------------------------------------------------------------------
print(summary(m4 <- lmer(dif ~ 1+S+F+SF + (1+S|subj), bs10,  REML=FALSE)), corr=FALSE)
anova(m3, m4, m0)

## ----versions------------------------------------------------------------
sessionInfo()

