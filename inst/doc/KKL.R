## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE-----------------------------------
library(lme4)
library(knitr)
library(RePsychLing)

opts_chunk$set(cache=FALSE)
options(show.signif.stars=FALSE,width=92)

## ----m0-----------------------------------------------------------------------------------
mm <- model.matrix(~ sze*(spt+obj+grv)*orn, data=KKL)
KKL$sze_spt <- mm[,  7]
KKL$sze_obj <- mm[,  8]
KKL$sze_grv <- mm[,  9]
KKL$sze_orn <- mm[, 10]
KKL$spt_orn <- mm[, 11]
KKL$obj_orn <- mm[, 12]
KKL$grv_orn <- mm[, 13]
KKL$sze_spt_orn <- mm[, 14]
KKL$sze_obj_orn <- mm[, 15]
KKL$sze_grv_orn <- mm[, 16]

m0 <- lmer(lrt ~ sze*(spt+obj+grv)*orn + 
                             (spt + obj + grv + orn + spt_orn + obj_orn + grv_orn | subj),
                         data=KKL, REML=FALSE, control=lmerControl(optCtrl=list(maxfun=10000L)))
#print(summary(m0), corr=FALSE)   

## ----sv0----------------------------------------------------------------------------------
chf0 <- getME(m0, "Tlist")[[1]]
zapsmall(chf0, digits=4)           

sv0 <- svd(chf0)

round(sv0$d^2/sum(sv0$d^2)*100, 1)  

## ----m1-----------------------------------------------------------------------------------
print(summary(m1 <- lmer(lrt ~ sze*(spt+obj+grv)*orn + 
                             (spt + obj + grv + orn + spt_orn + obj_orn + grv_orn || subj),
                         data=KKL, REML=FALSE)), corr=FALSE)

sv1 <- svd(diag(getME(m1, "theta")) )
sv1$d                                         

round(sv1$d^2/sum(sv1$d^2)*100, 1) 

anova(m1, m0)                 

## ----m2-----------------------------------------------------------------------------------
print(summary(m2 <- lmer(lrt ~ sze*(spt+obj+grv)*orn + 
                             (spt + obj + grv + orn + spt_orn || subj),
                         data=KKL, REML=FALSE)), corr=FALSE)

sv2 <- svd(diag(getME(m2, "theta")) )
sv2$d                                     

round(sv2$d^2/sum(sv2$d^2)*100, 3) 

anova(m2, m1)   # not significant: prefer m2 to m1
anova(m2, m0)   # significant: m2 is "reduced" too much               

## ----m3-----------------------------------------------------------------------------------
print(summary(m3 <- lmer(lrt ~ sze*(spt+obj+grv)*orn + 
                             (spt + obj + grv + orn + spt_orn | subj),
                         data=KKL, REML=FALSE)), corr=FALSE)

chf3 <- getME(m3, "Tlist")[[1]]
zapsmall(chf3, digits=4)           

sv3 <- svd(chf3)

round(sv3$d^2/sum(sv3$d^2)*100, 1)

anova(m2, m3)  # significant: prefer m3 to  m2 
anova(m3, m0)  # not significant: prefer m3 to m0

## ----m4-----------------------------------------------------------------------------------
print(summary(m4 <- lmer(lrt ~ sze*(spt+obj+grv)*orn +
                             (spt + grv | subj) + (0 + obj | subj) + (0 + orn | subj) + (0 + spt_orn | subj), 
                         data=KKL, REML=FALSE)), corr=FALSE)

getME(m4, "Tlist")      # Variance components look ok
chf4 <- diag(c(diag(getME(m4, "Tlist")[[1]]), getME(m4, "Tlist")[[2]], getME(m4, "Tlist")[[3]], getME(m4, "Tlist")[[4]]))
chf4[1:3, 1:3] <- getME(m4, "Tlist")[[1]]             

sv4 <- svd(chf4)
sv4$d                               # singular value decomposition: ok              

round(sv4$d^2/sum(sv4$d^2)*100, 1)  # percentages of variance accounted: ok

anova(m4, m3)    # not significant: prefer m4 to m3
anova(m4, m0)    # not significant: prefer m4 to m0    

## ----prof4, eval=FALSE--------------------------------------------------------------------
#  # Profiled 9 parameters individually; saved all results in as list in "p_m4.rda"
#  p_m4_1 = profile(m4_1, which=1)   # increment 1 to 9
#  confint(p_m4_1)

## ----versions-----------------------------------------------------------------------------
sessionInfo()

