# new intercept after shrinkage
lp.full <- predict(full.model)
lp.shrunk.full <- shrinkage.factor*lp.full
full.extra <- lrm.fit(y=mice::complete(imputed.data, m.choose)$HypoP,
                      offset=lp.shrunk.full) # TODO: for each imputation?
full.intercept.extra <- coef(full.extra)[1]
full.intercept.shrunk <- shrinkage.factor*full.model$coefficients["Intercept"]+full.intercept.extra

# check
mean(mice::complete(imputed.data, m.choose)$HypoP=="Yes") # event probability
mean(plogis(lp.shrunk.full+full.intercept.extra))  # same

# save results
an.full <- stats::anova(full.model)
coef.full.model <- shrinkage.factor*summary(full.model,
                                            dPTH=c(0, 1),
                                            CorrCa24u=c(0, 0.1),
                                            BSKgezien="No")[which(summary(full.model)[, "Type"]==1), c("Effect")]
Int.full.CI.lower <- full.model$coefficients["Intercept"]+qnorm(0.025)*sqrt(diag(full.model$var)[1])
Int.full.CI.upper <- full.model$coefficients["Intercept"]+qnorm(0.975)*sqrt(diag(full.model$var)[1])
CI.full.model <- exp(rbind(c(shrinkage.factor*Int.full.CI.lower+full.intercept.extra,
                             shrinkage.factor*Int.full.CI.upper+full.intercept.extra),
                           shrinkage.factor*cbind(coef.full.model+qnorm(0.025)*sqrt(diag(full.model$var)[-1]),
                                                  coef.full.model+qnorm(0.975)*sqrt(diag(full.model$var)[-1]))))
OR.full.model <- exp(c(full.intercept.shrunk, coef.full.model))

# new intercept after shrinkage
lp.final <- predict(final.model)
lp.shrunk.final <- shrinkage.factor*lp.final
final.extra <- lrm.fit(y=mice::complete(imputed.data, m.choose)$HypoP,
                      offset=lp.shrunk.final)
final.intercept.extra <- coef(final.extra)[1]
final.intercept.shrunk <- shrinkage.factor*final.model$coefficients["Intercept"]+final.intercept.extra

# check
mean(mice::complete(imputed.data, m.choose)$HypoP=="Yes") # event probability
mean(plogis(lp.shrunk.final+final.intercept.extra))  # same

# save results
an.final <- stats::anova(final.model)
coef.final.model <- shrinkage.factor*summary(final.model,
                                            dPTH=c(0, 1),
                                            CorrCa24u=c(0, 0.1),
                                            BSKgezien="No")[which(summary(final.model)[, "Type"]==1), c("Effect")]
Int.CI.lower <- final.model$coefficients["Intercept"]+qnorm(0.025)*sqrt(diag(final.model$var)[1])
Int.CI.upper <- final.model$coefficients["Intercept"]+qnorm(0.975)*sqrt(diag(final.model$var)[1])
CI.final.model <- exp(rbind(c(shrinkage.factor*Int.CI.lower+final.intercept.extra,
                              shrinkage.factor*Int.CI.upper+final.intercept.extra),
                            shrinkage.factor*cbind(coef.final.model+qnorm(0.025)*sqrt(diag(final.model$var)[-1]),
                                                   coef.final.model+qnorm(0.975)*sqrt(diag(final.model$var)[-1]))))
OR.final.model <- exp(c(final.intercept.shrunk, coef.final.model))

# new intercept after shrinkage
lp.simple <- predict(simple.model)
lp.shrunk.simple <- shrinkage.factor*lp.simple
simple.extra <- lrm.fit(y=mice::complete(imputed.data, m.choose)$HypoP,
                      offset=lp.shrunk.simple)
simple.intercept.extra <- coef(simple.extra)[1]
simple.intercept.shrunk <- shrinkage.factor*simple.model$coefficients["Intercept"]+simple.intercept.extra

# check
mean(mice::complete(imputed.data, m.choose)$HypoP=="Yes") # event probability
mean(plogis(lp.shrunk.simple+simple.intercept.extra))  # same

# save results
an.simple <- stats::anova(simple.model)
coef.simple.model <- shrinkage.factor*summary(simple.model,
                                              dPTH=c(0, 1))[which(summary(simple.model)[, "Type"]==1), c("Effect")]
Int.simple.CI.lower <- simple.model$coefficients["Intercept"]+qnorm(0.025)*sqrt(diag(simple.model$var)[1])
Int.simple.CI.upper <- simple.model$coefficients["Intercept"]+qnorm(0.975)*sqrt(diag(simple.model$var)[1])
CI.simple.model <- exp(rbind(c(shrinkage.factor*Int.simple.CI.lower+simple.intercept.extra,
                               shrinkage.factor*Int.simple.CI.upper+simple.intercept.extra),
                             shrinkage.factor*cbind(coef.simple.model+qnorm(0.025)*sqrt(diag(simple.model$var)[-1]),
                                                    coef.simple.model+qnorm(0.975)*sqrt(diag(simple.model$var)[-1]))))
OR.simple.model <- exp(c(simple.intercept.shrunk, coef.simple.model))
