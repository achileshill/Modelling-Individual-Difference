library(readxl)
library(ltm)
library(psych)
library(lavaan)
df <- readxl::read_xlsx("DarkTriad_Data File.xlsx")
narcissism <- df[,1:10]
machiavellianism <- df[,11:20]
psychopathy <- df[,21:30]

# Parallel Analysis
pcor_nar <- polychoric(narcissism)$rho
fa.parallel(pcor_nar, n.obs=117,fa='pc', main='Narcissism Parallel Analysis Scree Plots')

pcor_mach <- polychoric(machiavellianism)$rho
fa.parallel(pcor_mach, n.obs=117, fa='pc', main='Machiavellianism Parallel Analysis Scree Plots')

pcor_psyc <- polychoric(psychopathy)$rho
fa.parallel(pcor_psyc, n.obs=117, fa='pc', main='Psychopathy Parallel Analysis Scree Plots')

# Test of Equality of Discrimination
fit_nar <- grm(narcissism)
fit_nar_constrain <- grm(narcissism, constrained = T)
anova(fit_nar_constrain, fit_nar)

fit_mach <- grm(machiavellianism)
fit_mach_constrain <- grm(machiavellianism, constrained=T)
anova(fit_mach_constrain, fit_mach)

df_integer <- data.frame(lapply(df, function(x) factor(as.integer(x))))
psy1 <- df_integer[,21:30]
fit_psyc <- grm(psy1)
fit_psyc_constrain <-  grm(psy1, constrained = T)
anova(fit_psyc_constrain, fit_psyc)

# Test Information Function
plot(fit_nar, type='IIC', items=0, main='TIF Narcissism')
plot(fit_mach, type='IIC', items=0, main='TIF Machiavellianism')
plot(fit_psyc, type='IIC', items=0, main='TIF Psychopathy')
plot(fit_nar, type='IIC', items=c(1,2,6), main='Item 2 & 6 Poor Information Function')
plot(fit_mach, type='IIC', items=c(1,3,6,8,9,10), main='Item 13,16,18,19 & 20 Poor Information Function')
# Item-rest Correlation 
item_rest_corr <- function(data) {
  sapply(1:ncol(data), function(i) {
    cor(data[, i], rowSums(data[, -i], na.rm = TRUE), use = "pairwise.complete.obs")
  })
}

item_rest_corr(narcissism)
item_rest_corr(machiavellianism)
item_rest_corr(psychopathy)

# CTT Reliability estimate
# Reverse contraindicative items
nar_dup <- narcissism
neg_items_nar <- c(2,4,10)
nar_dup[, neg_items_nar] <- 6 - nar_dup[, neg_items_nar]
psych_dup <- psychopathy
neg_items_psyc <- c(2,4,6)
psych_dup[, neg_items_psyc] <- 6 - psych_dup[, neg_items_psyc]
neg_items_mach <- c(6,8,10)
machi_dup <- machiavellianism
machi_dup[, neg_items_mach] <- 6 - machi_dup[, neg_items_mach]

# Alpha and glb
alpha(nar_dup)$total
glb.algebraic(cov(nar_dup))[1]
alpha(machi_dup)$total
glb.algebraic(cov(machi_dup))[1]
alpha(psych_dup)$total
glb.algebraic(cov(psych_dup))[1]

# Re-estimate Reliability after removal
new_nar <- nar_dup[,c(1,3,4,5,7,8,9)]
new_mach <- machi_dup[, c(1,2,4,5,7)]
alpha(new_nar)$total
alpha(new_mach)$total
glb.algebraic(cov(new_nar))[1]
glb.algebraic(cov(new_mach))[1]


data <- df[,c(1,3,4,5,7,8,9,10,11,12,14,15,17,21,22,23,24,25,26,27,28,29,30)]

# CFA
# Comparison matrix
M <- matrix(nrow=4, ncol=10)
colnames(M) <- c('Model','Description','Chi-square','df','#par','RMSEA','CFI','TLI','BIC','AIC')
model <- '
Narcissism =~Q1+Q3+Q4+Q5+Q7+Q8+Q9+Q10
Machiavellianism =~Q11+Q12+Q14+Q15+Q17
Psychopathy =~ Q21+Q22+Q23+Q24+Q25+Q26+Q27+Q28+Q29+Q30
'
fit <- cfa(model, data=data, sample.nobs = 117)
M[1,c(1,2)] <- c('fit','Three-factor model')
M[1,c(3,4,5,6,7,8,9)] <- round(as.numeric(fitMeasures(fit,c('chisq','df','npar','rmsea','cfi','tli','bic'))),4)
M[1,10] <- round(as.numeric(-2*fitMeasures(fit,'logl')+2*fitMeasures(fit,'npar')),4)

max(modindices(fit)[,4])
modindices(fit)[which.max(modindices(fit)[,4]),c(1:4)]

# MI for Narcissism and Q28 is moderate, refine model to check if it fit better

model_fix1 <- '
Narcissism =~Q1+Q3+Q4+Q5+Q7+Q8+Q9+Q10+Q28
Machiavellianism =~Q11+Q12+Q14+Q15+Q17
Psychopathy =~ Q21+Q22+Q23+Q24+Q25+Q26+Q27+Q28+Q29+Q30
'
fit_1 <- cfa(model_fix1, data=data, sample.nobs = 117)
M[2,c(1,2)] <- c('fit_revised1','Three-factor model with cross-loading Narcissism and Q28')
M[2,c(3,4,5,6,7,8,9)] <- round(as.numeric(fitMeasures(fit_1,c('chisq','df','npar','rmsea','cfi','tli','bic'))),4)
M[2,10] <- round(as.numeric(-2*fitMeasures(fit_1,'logl')+2*fitMeasures(fit_1,'npar')),4)
max(modindices(fit_1)[,4])
modindices(fit_1)[which.max(modindices(fit_1)[,4]),c(1:4)]

# Moderate MI of residual covariance between Q5 and Q15

model_fix2 <- '
Narcissism =~Q1+Q3+Q4+Q5+Q7+Q8+Q9+Q10+Q28
Machiavellianism =~Q11+Q12+Q14+Q15+Q17
Psychopathy =~ Q21+Q22+Q23+Q24+Q25+Q26+Q27+Q28+Q29+Q30
Q5 ~~ Q15
'
fit_2 <- cfa(model_fix2, data=data, sample.nobs = 117)
M[3,c(1,2)] <- c('fit_revised2','Three-factor model with cross-loading Narcissism and Q28 and residual covariance between Q5 & Q15')
M[3,c(3,4,5,6,7,8,9)] <- round(as.numeric(fitMeasures(fit_2,c('chisq','df','npar','rmsea','cfi','tli','bic'))),4)
M[3,10] <- round(as.numeric(-2*fitMeasures(fit_2,'logl')+2*fitMeasures(fit_2,'npar')),4)

max(modindices(fit_2)[,4])
modindices(fit_2)[which.max(modindices(fit_2)[,4]),c(1:4)]

model_fix3 <- '
Narcissism =~Q1+Q3+Q4+Q5+Q7+Q8+Q9+Q10+Q28
Machiavellianism =~Q11+Q12+Q14+Q15+Q17
Psychopathy =~ Q21+Q22+Q23+Q24+Q25+Q26+Q27+Q28+Q29+Q30
Q5 ~~ Q15
Q4 ~~ Q10
'
fit_3 <- cfa(model_fix3, data=data, sample.nobs = 117)
M[4,c(1,2)] <- c('fit_revised3','Three-factor model with cross-loading Narcissism and Q28 and residual covariance between Q5 & Q15 and Q4 & Q10')
M[4,c(3,4,5,6,7,8,9)] <-round(as.numeric(fitMeasures(fit_3,c('chisq','df','npar','rmsea','cfi','tli','bic'))),4)
M[4,10] <- round(as.numeric(-2*fitMeasures(fit_3,'logl')+2*fitMeasures(fit_3,'npar')),4)

max(modindices(fit_3)[,4])
modindices(fit_3)[which.max(modindices(fit_3)[,4]),c(1:4)]

M

# EFA
M_EFA <- matrix(nrow=5, ncol=6)
cor_data <- cor(data)
fit_EFA_1f <- fa(r=cor_data,nfactor=1,n.obs=117,fm="ml",rotate="promax")
fit_EFA_2f <- fa(r=cor_data,nfactor=2,n.obs=117,fm="ml",rotate="promax")
fit_EFA_3f <- fa(r=cor_data,nfactor=3,n.obs=117,fm="ml",rotate="promax")
fit_EFA_4f <- fa(r=cor_data,nfactor=4,n.obs=117,fm="ml",rotate="promax")
fit_EFA_5f <- fa(r=cor_data,nfactor=5,n.obs=117,fm="ml",rotate="promax")
M_EFA[c(1,2,3,4,5),1] <- c(1,2,3,4,5)
colnames(M_EFA) <- c('nfactor','stat','df','pval','BIC','RMSEA')
M_EFA[1,c(2,3,4,5,6)] <- c(fit_EFA_1f$STATISTIC, fit_EFA_1f$dof, fit_EFA_1f$PVAL,fit_EFA_1f$BIC, fit_EFA_1f$RMSEA[1])
M_EFA[2,c(2,3,4,5,6)] <- c(fit_EFA_2f$STATISTIC, fit_EFA_2f$dof, fit_EFA_2f$PVAL,fit_EFA_2f$BIC, fit_EFA_2f$RMSEA[1])
M_EFA[3,c(2,3,4,5,6)] <- c(fit_EFA_3f$STATISTIC, fit_EFA_3f$dof, fit_EFA_3f$PVAL,fit_EFA_3f$BIC, fit_EFA_3f$RMSEA[1])
M_EFA[4,c(2,3,4,5,6)] <- c(fit_EFA_4f$STATISTIC, fit_EFA_4f$dof, fit_EFA_4f$PVAL,fit_EFA_4f$BIC, fit_EFA_4f$RMSEA[1])
M_EFA[5,c(2,3,4,5,6)] <- c(fit_EFA_5f$STATISTIC, fit_EFA_5f$dof, fit_EFA_5f$PVAL,fit_EFA_5f$BIC, fit_EFA_5f$RMSEA[1])
M_EFA

# Assess factor structure
fit_EFA_3f
fit_EFA_4f
fit_EFA_5f

