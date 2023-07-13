#Set working directory
setwd("//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/1. SCHILDKLIER onderzoek/PROJECT Predictiemodel [Sam]/R")
library(broom)
library(moments)
library(caret)
library(calibrator)
library(sampling)
library(DescTools)
library(CVST)
library(boot)
library(dplyr,lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(tidyr, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(ggplot2, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(cowplot, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") #version 0.9.4
library(survival, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(survminer, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(tableone, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(reshape2, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(rms, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(cmprsk, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(gridExtra, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(psfmi, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(PRROC, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(ROCR, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(foreign, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") # Reading external (spps) data
library(rms, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") # Harrell's regression library
library(mice, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") # Multiple imputation
library(mitools, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(ggplot2, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(survminer, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(haven, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") #package needed to open spss files
library(Hmisc, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(ResourceSelection, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R") #package needed for hosmerlemeshow test
library(remotes, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(pROC, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(lazyeval, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(rex, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(riskmetric, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")
library(ValidationMetrics, lib.loc = "//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/Persoonlijk vDijk/R")


#Variable selection
cvar <- c("Record_ID", "Center","Readmission","HypoP", "surgery_type", "Sex", "Indicatie",
          "BSKgezien",  "CHKD", "BSKinPA", "PostPA", "PostOPCa", "Validatie", "Validatie1")
evar <- c("")
nvar <- c("BaselinePTH", "PTH24u", "BaseCa", "Basealbu", "Ca24u", "Albu24u",
           "Age_Years")
dvar <- c("")

#import dataset
Data <- read.csv("//storage/x/avcl05/HEEL/Exchange_X/Research Endocrien/1. SCHILDKLIER onderzoek/PROJECT Predictiemodel [Sam]/R/Data - Copy.csv",
                 stringsAsFactors=TRUE,
                 header = T, sep = ";", dec = ".", na.strings = c("", " ", "NA", "999", 999, "999.00", 999.00, "Missing", "missing"))
data <- as_tibble(Data) # TODO: "data" wordt niet gebruikt?

#Mutate date variables
Data <- Data %>% mutate_at(dvar, ~as.Date(.,"%d-%m-%Y")) # DISCUSS: dvar is leeg?
DAT <- Data[c(cvar, nvar)]
str(DAT)
summary(DAT$Validatie1) #variabele voor de leave-one-out validatie

#Data cleaning
DAT$HypoP <- ifelse(DAT$HypoP %in% c("No ", "No"), "No", "Yes")
DAT$HypoP <- as.factor(DAT$HypoP)
DAT$PostOPCa <- ifelse(DAT$PostOPCa %in% c("No ", "No"), "No", "Yes")
DAT$PostOPCa <- as.factor(DAT$PostOPCa)
DAT$BaselinePTH <- gsub(",", ".", DAT$BaselinePTH)
DAT$BaselinePTH <- as.numeric(DAT$BaselinePTH)
DAT$PTH24u <- gsub(",", ".", DAT$PTH24u)
DAT$PTH24u <- as.numeric(DAT$PTH24u)
DAT$BaseCa <- gsub(",", ".", DAT$BaseCa)
DAT$BaseCa <- as.numeric(DAT$BaseCa)
DAT$Basealbu <- gsub(",", ".", DAT$Basealbu)
DAT$Basealbu <- as.numeric(DAT$Basealbu)
DAT$Ca24u <- gsub(",", ".", DAT$Ca24u)
DAT$Ca24u <- as.numeric(DAT$Ca24u)
DAT$Albu24u <- gsub(",", ".", DAT$Albu24u)
DAT$Albu24u <- as.numeric(DAT$Albu24u)
DAT$Age_Years <- as.numeric(DAT$Age_Years)

#Imputation---------------
md.pattern(DAT)
dd<-datadist(DAT)
options(datadist='dd')
options(digits=8)
miceHypoP <- mice(DAT, maxit = 0)
miceHypoPmeth <- miceHypoP$meth
miceHypoPpred <- miceHypoP$pred
miceHypoP

## will not be used as predictor in model
miceHypoPpred[, "Record_ID"] <- 0
miceHypoPpred[, "Validatie"] <- 0
miceHypoPpred[, "Validatie1"] <- 0
# miceHypoPpred[, "surgery_type"] <- 0
# miceHypoPpred[, "Age_Years"] <- 0
# miceHypoPpred[, "Sex"] <- 0
miceHypoPpred[, "Indicatie"] <- 0
miceHypoPpred[, "PostOPCa"] <- 0

## will not be imputed
miceHypoPmeth["Record_ID"] <- ""
miceHypoPmeth["Center"] <- ""
miceHypoPmeth["Validatie"] <- ""
miceHypoPmeth["Validatie1"] <- ""
miceHypoPmeth["Readmission"] <- ""
miceHypoPmeth["HypoP"] <- ""
miceHypoPmeth["surgery_type"] <- ""
miceHypoPmeth["Age_Years"] <- ""
miceHypoPmeth["Sex"] <- ""
miceHypoPmeth["Indicatie"] <- ""
miceHypoPmeth["PostOPCa"] <- ""

mi <- mice(DAT, maxit = 50, m = 5, seed = 123,
           method=miceHypoPmeth,
           predictorMatrix=miceHypoPpred,
           vis="monotone", print=FALSE)

#Berekenen PTH daling, gecorrigeerd calcium en gecorrigeerd calcium daling
mi_long <- complete(mi, action = "long", include = T)
mi_long$dPTH <- ((mi_long$BaselinePTH - mi_long$PTH24u)/mi_long$BaselinePTH*100)
mi_long$CorrCaBaseline <- mi_long$BaseCa +((34-mi_long$Basealbu)*0.016)
mi_long$CorrCa24u <- mi_long$Ca24u +((34-mi_long$Albu24u)*0.016)
mi_long$dCorrCa24u <- ((mi_long$CorrCaBaseline - mi_long$CorrCa24u)/mi_long$CorrCaBaseline*100)
mi_long$dCa24u <- ((mi_long$BaseCa - mi_long$Ca24u)/mi_long$BaseCa*100)
mi_long$BSKnietgezien <- ifelse(mi_long$BSKgezien=="Yes", "No", ifelse(mi_long$BSKgezien=="No", "Yes", NA))
mi_long$BSKnietgezien <- as.factor(mi_long$BSKnietgezien)

mi_test <- as.mids(mi_long)
minew <- filter(mi_long, .imp > 0) #code om originele dataset uit de gehele dataset te halen
dat.complete <- complete(mi_test, 1) # Eerste geimputeerde dataset

# DISCUSS: only set datadist for original data, not for imputed data
dd <- datadist(minew)
options(datadist = "dd")

#linearity assumption visual inspection------------
# Fit logistic regression model
#Hieronder alle mogelijke predictors. Note (!): PTH24u en dPTH lijken veel op elkaar en kunnen dus eigenlijk niet in hetzelfde model worden meegenomen. Hetzelfde geldt voor Ca24u + CorrCa24u + dCa24u + dCorrCa24u.
formula <- HypoP ~ dPTH + surgery_type + PTH24u + Ca24u + CorrCa24u + dCa24u + dCorrCa24u + Age_Years + Sex + Indicatie + BSK gezien + CHKD
full.model <- glm(formula, data=dat.complete, family="binomial")

## backward selection

model <- glm(HypoP ~ dPTH + dCorrCa24u + BSKnietgezien, data = dat.complete, family = binomial)
# Predict the probability (p) of diabete positivity
probabilities <- predict(model, type = "response")
mydata <- dat.complete %>%
  dplyr::select_if(is.numeric)
predictors <- colnames(mydata)
# Bind the logit and tidying the data for plot
mydata <- dat.complete %>%
  dplyr::select(dPTH, dCorrCa24u) %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") +
  theme_bw() +
  facet_wrap(~predictors, scales = "free_y")

#dPTH linear, dCorrCa24u does not seem to be linear. Hoe aanpakken? Splines? Transformation?

#outliers
plot(model, which = 4, id.n = 3)
model.data <- augment(model) %>%
  mutate(index = 1:n())
model.data %>% top_n(3, .cooksd)
ggplot(model.data, aes(index, .std.resid)) +
  geom_point(aes(color = HypoP), alpha = .5) +
  theme_bw()
model.data %>%
  filter(abs(.std.resid) > 3.0)
#There is no influential observations in our data.

#multicolinearity. As a rule of thumb, a VIF value that exceeds 5 or 10 indicates a problematic amount of collinearity.
car::vif(model)
#In our example, there is no collinearity: all variables have a value of VIF well below 5.


#univariate univariable logistic regression model-------------
rAGE <- fit.mult.impute(HypoP ~ Age_Years,
                        data = mi_test, xtrans = mi_test, fitter = lrm)
rAGE
exp(coef(rAGE))
lower_rAGE <- exp(0.0150-1.96*0.0124)
higher_rAGE <- exp(0.0150+1.96*0.0124)
lower_rAGE
higher_rAGE

rSEX <- fit.mult.impute(HypoP ~ Sex,
                        data = mi_test, xtrans = mi_test, fitter = lrm)
rSEX
exp(coef(rSEX))
lower_rSEX <- exp(0.7342-1.96*0.3931)
higher_rSEX <- exp(0.7342+1.96*0.3931)
lower_rSEX
higher_rSEX

rSurg <- fit.mult.impute(HypoP ~ surgery_type,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
rSurg
exp(coef(rSurg))
lower_rSurg <- exp(0.4520-1.96*0.5120)
higher_rSurg <- exp(0.4520+1.96*0.5120)
lower_rSurg
higher_rSurg

rCHKD <- fit.mult.impute(HypoP ~ CHKD,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
rCHKD
exp(coef(rCHKD))
lower_rCHKD <- exp(1.0063-1.96*0.3899)
higher_rCHKD <- exp(1.0063+1.96*0.3899)
lower_rCHKD
higher_rCHKD

rBSKnietgezien <- fit.mult.impute(HypoP ~ BSKnietgezien,
                                  data = mi_test, xtrans = mi_test, fitter = lrm, x=TRUE, y=TRUE)
rBSKnietgezien
exp(coef(rBSKnietgezien))
lower_rBSKnietgezien <- exp(1.2182-1.96*0.4210)
higher_rBSKnietgezien <- exp(1.2182+1.96*0.4210)
lower_rBSKnietgezien
higher_rBSKnietgezien


rdPTH <- fit.mult.impute(HypoP ~ dPTH,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
rdPTH
exp(coef(rdPTH))
lower_rdPTH <- exp(0.1025-1.96*0.0273)
higher_rdPTH <- exp(0.1025+1.96*0.0273)
lower_rdPTH
higher_rdPTH

rdcCa <- fit.mult.impute(HypoP ~ dCorrCa24u,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
rdcCa
exp(coef(rdcCa))
lower_rdcCa <- exp(0.1252-1.96*0.0336)
higher_rdcCa <- exp(0.1252+1.96*0.0336)
lower_rdcCa
higher_rdcCa

rPostOPCa <- fit.mult.impute(HypoP ~ PostOPCa,
                             data = mi_test, xtrans = mi_test, fitter = lrm)
rPostOPCa
exp(coef(rPostOPCa))
lower_rPostOPCa <- exp(1.0245-1.96*0.3984)
higher_rPostOPCa <- exp(1.0245+1.96*0.3984)
lower_rPostOPCa
higher_rPostOPCa

rBSKinPA <- fit.mult.impute(HypoP ~ BSKinPA,
                            data = mi_test, xtrans = mi_test, fitter = lrm)
rBSKinPA
exp(coef(rBSKinPA))
lower_rBSKinPA <- exp(1.3523-1.96*0.3924)
higher_rBSKinPA <- exp(1.3523+1.96*0.3924)
lower_rBSKinPA
higher_rBSKinPA


#Multivariable model-------------
a <- fit.mult.impute(HypoP ~ dPTH + dCorrCa24u + BSKnietgezien,
                     data = mi_test, xtrans = mi_test, fitter = lrm)
a

#Calculate Shrinkage factor and optimism-adjusted performance-----------------------------------
# Fit the logistic regression model on the original dataset
fit <- fit.mult.impute(HypoP ~ dPTH + dCorrCa24u + BSKnietgezien,
                       data = mi_test, xtrans = mi_test, fitter = lrm, x = TRUE, y = TRUE)
# Calculate the shrinkage factor using the validate function with 100 bootstrap replicates
set.seed(15)
results <- validate(fit, method = "boot", B = 1000)

# Inspect the results
summary(results)
print(results)
shrinkage_factor <- 0.9369
shrunken_ORs <- exp(shrinkage_factor * fit$coefficients) #Intercept ook?

#performance model unadjusted
lp<-NULL
for (i in 1:5){
  lp.j<-predict(a,newdata=complete(mi_test,i))
  lp<-cbind(lp,lp.j)
}
ValidationMetrics::val.prob.mi(lp.mi=lp, y=as.numeric(mi_test$data$HypoP)-1)

#performance model adjusted for optimism with shrinkage
lp_s <- lp*0.9369
ValidationMetrics::val.prob.mi(lp.mi=lp_s, g=5, y=as.numeric(mi_test$data$HypoP)-1)

#Plotje met 3 ROC-curves (unadjusted for optimism)-----------------------------
a <- fit.mult.impute(HypoP ~ dPTH + dCorrCa24u + BSKnietgezien,
                     data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities <- as.vector(predict(a, type = "fitted"))
outcomes <- mi_test$data$HypoP
roc <- roc(outcomes, probabilities)

rdPTH <- fit.mult.impute(HypoP ~ dPTH,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities_rdPTH <- as.vector(predict(rdPTH, type = "fitted"))
outcomes_rdPTH <- mi_test$data$HypoP
roc_rdPTH <- roc(outcomes_rdPTH, probabilities_rdPTH)

rdcCa <- fit.mult.impute(HypoP ~ dCorrCa24u,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities_rdcCa <- as.vector(predict(rdcCa, type = "fitted"))
outcomes_rdcCa <- mi_test$data$HypoP
roc_rdcCa <- roc(outcomes_rdcCa, probabilities_rdcCa)

auc_final <- auc(outcomes, probabilities)
auc_rdPTH <- auc(outcomes_rdPTH, probabilities_rdPTH)
auc_rdcCa <- auc(outcomes_rdcCa, probabilities_rdcCa)
df_roc <- rbind(data.frame(FPR = 1 - roc$specificities, TPR = roc$sensitivities, Model = "Final Model"),
                data.frame(FPR = 1 - roc_rdPTH$specificities, TPR = roc_rdPTH$sensitivities, Model = "rdPTH Model"),
                data.frame(FPR = 1 - roc_rdcCa$specificities, TPR = roc_rdcCa$sensitivities, Model = "rdcCa Model"))

# Create ROC plot with AUC annotations
ggplot(data = df_roc, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "ROC Curve",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Model") +
  annotate("text", x = 0.13, y = 0.9, label = paste0("Final Model\nC-statistic: ", round(auc_final, 3)), size = 4, color= "red") +
  annotate("text", x = 0.4, y = 0.8, label = paste0("rdPTH Model\nC-statistic: ", round(auc_rdPTH, 3)), size = 4, color= "royalblue2") +
  annotate("text", x = 0.35, y = 0.5, label = paste0("rdcCa Model\nC-statistic: ", round(auc_rdcCa, 3)), size = 4, color= "forestgreen")

#Plotje met 4 ROC-curves (unadjusted for optimism)------------------------------
a <- fit.mult.impute(HypoP ~ dPTH + dCorrCa24u + BSKnietgezien,
                     data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities <- as.vector(predict(a, type = "fitted"))
outcomes <- mi_test$data$HypoP
roc <- roc(outcomes, probabilities)

rdPTH <- fit.mult.impute(HypoP ~ dPTH,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities_rdPTH <- as.vector(predict(rdPTH, type = "fitted"))
outcomes_rdPTH <- mi_test$data$HypoP
roc_rdPTH <- roc(outcomes_rdPTH, probabilities_rdPTH)

rdcCa <- fit.mult.impute(HypoP ~ dCorrCa24u,
                         data = mi_test, xtrans = mi_test, fitter = lrm)
probabilities_rdcCa <- as.vector(predict(rdcCa, type = "fitted"))
outcomes_rdcCa <- mi_test$data$HypoP
roc_rdcCa <- roc(outcomes_rdcCa, probabilities_rdcCa)

rBSKnietgezien <- fit.mult.impute(HypoP ~ BSKnietgezien,
                                  data = mi_test, xtrans = mi_test, fitter = lrm, x=TRUE, y=TRUE)
probabilities_rBSKnietgezien <- as.vector(predict(rBSKnietgezien, type = "fitted"))
outcomes_rBSKnietgezien <- mi_test$data$HypoP
roc_rBSKnietgezien <- roc(outcomes_rBSKnietgezien, probabilities_rBSKnietgezien)

auc_final <- auc(outcomes, probabilities)
auc_rdPTH <- auc(outcomes_rdPTH, probabilities_rdPTH)
auc_rdcCa <- auc(outcomes_rdcCa, probabilities_rdcCa)
auc_rBSKnietgezien <- auc(outcomes_rBSKnietgezien, probabilities_rBSKnietgezien)

df_roc <- rbind(data.frame(FPR = 1 - roc$specificities, TPR = roc$sensitivities, Model = "Combined model"),
                data.frame(FPR = 1 - roc_rdPTH$specificities, TPR = roc_rdPTH$sensitivities, Model = "dPTH"),
                data.frame(FPR = 1 - roc_rdcCa$specificities, TPR = roc_rdcCa$sensitivities, Model = "dcCa"),
                data.frame(FPR = 1 - roc_rBSKnietgezien$specificities, TPR = roc_rBSKnietgezien$sensitivities, Model = "Parathyroid identification"))


ggplot(data = df_roc, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "ROC Curve",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Model") +
  annotate("text", x = 0.13, y = 0.9, label = paste0("Combined model\nC-statistic: ", round(auc_final, 3)), size = 4, color= "red") +
  annotate("text", x = 0.37, y = 0.8, label = paste0("dPTH\nC-statistic: ", round(auc_rdPTH, 3)), size = 4, color= "royalblue2") +
  annotate("text", x = 0.4, y = 0.62, label = paste0("dcCa\nC-statistic: ", round(auc_rdcCa, 3)), size = 4, color= "forestgreen") +
  annotate("text", x = 0.22, y = 0.2, label = paste0("Parathyroid identification\nC-statistic: ", round(auc_rBSKnietgezien, 3)), size = 4, color= "orange") +
  scale_color_manual(values = c("red", "forestgreen", "royalblue2", "orange"))

#Internal-external cross-validation-------

   #Hier mag je aan de slag :)
   #Laat me weten als je additionele informatie nodig hebt, bijna altijd bereikbaar. Bellen vind ik top!


