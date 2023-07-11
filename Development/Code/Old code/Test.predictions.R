# load data
rm(list = ls(all.names = TRUE))
library(rms)
load(file=paste0("Z:/Project Predict Hypoparathyroidism/Development/Data/test.Rdata"))
m <- 5
dd<-rms::datadist(single.imputation)
options(datadist='dd')
options(digits=8)

# SINGLE IMPUTATION WITHOUT SHRINKAGE
test.model <- rms::lrm(HypoP ~ dPTH + BSKgezien + CorrCa24u,
                       data=single.imputation, x=TRUE, y=TRUE)
patient.nr <- 4
test.patient <- cbind(single.imputation[, c("BaselinePTH", "PTH24u", "dPTH",
                                            "Ca24u", "Albu24u", "CorrCa24u",
                                            "BSKgezien")],
                      predict(test.model,
                              data=single.imputation[, c("dPTH",
                                                         "CorrCa24u",
                                                         "BSKgezien")],
                              type="lp"),
                      predict(test.model,
                              data=single.imputation[, c("dPTH",
                                                         "CorrCa24u",
                                                         "BSKgezien")],
                              type="fitted"))[patient.nr, ]
test.coef <- summary(test.model,
                     dPTH=c(1, 2),
                     CorrCa24u=c(1, 2),
                     BSKgezien="No")[which(summary(test.model)[, "Type"]==1),
                                     c("Effect")]
# lp = 1.5003 + 0.0819*dPTH - 1.5727*BSKgezien=Yes - 4.2278*CorrCa24u
# -0.76812838
lp <- test.model$coefficients["Intercept"]+
  test.coef["dPTH"]*test.patient["dPTH"]+
  test.coef["CorrCa24u"]*test.patient["CorrCa24u"]+
  test.coef["BSKgezien - Yes:No"]*as.numeric(test.patient["BSKgezien"]=="Yes")
# p = exp(lp)/(1+exp(lp))
# 31.688411
p <- exp(as.numeric(lp))/(1+exp(as.numeric(lp)))*100

# MULTIPLE IMPUTATION AND WITHOUT SHRINKAGE
test.model.mi <- Hmisc::fit.mult.impute(HypoP ~ dPTH + BSKgezien + CorrCa24u,
                                        lrm,
                                        xtrans=imputed.data,
                                        data=work.data,
                                        n.impute=m, pr=FALSE, fit.reps=TRUE,
                                        fitargs=list(y=TRUE, x=TRUE, se.fit=TRUE, maxit=1000))
patient.nr <- 3
test.coef.mi <- summary(test.model.mi,
                     dPTH=c(1, 2),
                     CorrCa24u=c(1, 2),
                     BSKgezien="No")[which(summary(test.model.mi)[, "Type"]==1),
                                     c("Effect")]
# lp = 0.1755 + 0.0929*dPTH - 1.5208*BSKgezien=Yes - 4.0640*CorrCa24u
lp.mi <- test.model.mi$coefficients["Intercept"]+
  test.coef.mi["dPTH"]*test.patient["dPTH"]+
  test.coef.mi["CorrCa24u"]*test.patient["CorrCa24u"]+
  test.coef.mi["BSKgezien - Yes:No"]*as.numeric(test.patient["BSKgezien"]=="Yes")
p.mi <- exp(as.numeric(lp.mi))/(1+exp(as.numeric(lp.mi)))*100

# MULTIPLE IMPUTATION AND SHRINKAGE
shrinkage.factor <- 0.8740229

lp.shrunk <- shrinkage.factor*predict(test.model.mi, type="lp")
fit.extra <- lrm.fit(y=single.imputation$HypoP,
                     offset=lp.shrunk) # TODO: for each imputation?
intercept.extra <- coef(fit.extra)[1]
intercept.shrunk <- shrinkage.factor*test.model.mi$coefficients["Intercept"]+intercept.extra

test.coef.mi.shr <- shrinkage.factor*test.coef.mi

lp.mi.shr <- intercept.shrunk+
  test.coef.mi.shr["dPTH"]*test.patient["dPTH"]+
  test.coef.mi.shr["CorrCa24u"]*test.patient["CorrCa24u"]+
  test.coef.mi.shr["BSKgezien - Yes:No"]*as.numeric(test.patient["BSKgezien"]=="Yes")
p.mi.shr <- exp(as.numeric(lp.mi.shr))/(1+exp(as.numeric(lp.mi.shr)))*100

print(lp)        # -0.76812838
print(lp.mi)     # -0.6786823
print(lp.mi.shr) # -0.68060039

print(p)         # 31.688411
print(p.mi)      # 33.655546
print(p.mi.shr)  # 33.612731
