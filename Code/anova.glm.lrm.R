# LRT test with glm
drop1(glm(form.full.model, data=mice::complete(imputed.data, m), family="binomial"), test="Chisq")
# Wald-test = coef / SE^2
an.full <- anova(lrm(form.full.model, data=mice::complete(imputed.data, m)))
