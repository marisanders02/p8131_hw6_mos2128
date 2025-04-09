##########################################################################################################
## Ex1: linear mixed effects models of the orthodontic example
# Note: we use the R pacakge nlme here. An alternative is the lme4 package.
# Note: The lmer function in lme4 is more suitable for modeling multiple non-nested random effects
##########################################################################################################
library(nlme)
library(ggplot2)
library (lattice)
Orthodont<- read.table ("orthodontic.dat",header=TRUE)
head(Orthodont)

Orth.new <- groupedData (distance ~ age | child,data = as.data.frame (Orthodont)) # not necessary for lme
head(Orth.new)
OrthFem <- subset(Orth.new, male==0)
OrthFem[1:5,]
plot(OrthFem) # default: ordered by max resp of each child
ggplot(Orth.new, aes(age, distance, group = child)) + 
  geom_line() +
  geom_point() +
  theme_classic()

######################################
# fit a random intercept model
LMM1 <- lme (distance ~ age, random = ~1 | child,  data = OrthFem, method='REML')  ## random intercept: ~1
summary (LMM1) # pay attention to: random effects, fixed effects, 
#
VarCorr(LMM1) # covariance estimates for random effects and variance for residuals
LMM1$sigma # std for residuals
vcov(LMM1) # covariance for fixed effects estimates (inv fisher info)
#
fixed.effects(LMM1) # fixed effects coeff: beta_0, beta_1
random.effects(LMM1) # ordered random effects, BLUP (in this case, just b_i), random effects for different subjects, row index are the subject index. 10th subject has the smallest random effects. 
fitted(LMM1) # fixed+random for each subj in each visit, beta_0 + b_i + beta_1*age
OrthFem$distance-fitted(LMM1) # residuals = measured values - fitted values. 4 residues per subject (we only subset male=0 subjects)
LMM1$residuals
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for subj1
points(OrthFem$age[1:4],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[1:4], pch=4,col=1) # global estimates for these four time point. beta_0 + beta_1 * age
lines(OrthFem$age[1:4],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[1:4], pch=4,col=1) 
lines(OrthFem$age[1:4],fitted(LMM1)[1:4], lty=5,type='b',pch=4,col=26) # subject-specific profiles
random.effects(LMM1) # check random effects b_i for the first subj: -1.229
#
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) # original data for subj2
points(OrthFem$age[5:8],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM1)[5:8], lty=5,type='b',pch=4,col=20)
random.effects(LMM1) #check random effects b_i for the second subj: 0.340


# check equivalence to marginal model with compound symmetry correlation structure
summary(gls(distance~age, OrthFem, correlation=corCompSymm(form = ~ 1 |child), method="REML"))
# check rho==sigma_b^2/(sigma_b^2+sigma^2)
## VarCorr(LMM1)
## sigma_b^2/(sigma_b^2+sigma^2): 4.2785689/(4.2785689+0.6084517) = 0.8754964 (==rho)


###########################################
# compare models (likelihood ratio test)
LMM.1 <- lme (distance ~ age, random = ~ 1 | child,  data = OrthFem, method='ML')  # do NOT use REML for likelihood ratio
LMM.2 <- lme (distance ~ 1, random = ~ 1 | child,  data = OrthFem, method='ML')
anova(LMM.2,LMM.1) 


#############################################
# fit a random intercept and slope model
LMM2 <- lme (distance ~ age, random = ~ 1+ age | child, data = OrthFem) # random effect include intercept and for age
summary (LMM2) 
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for child1
points(OrthFem$age[1:4],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[1:4], pch=4,col=1) # global
lines(OrthFem$age[1:4],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[1:4], pch=4,col=1) 
lines(OrthFem$age[1:4],fitted(LMM2)[1:4], lty=5,type='b',pch=4,col=26) 
random.effects(LMM2) # check BLUP 
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) 
points(OrthFem$age[5:8],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM2)[5:8], lty=5,type='b',pch=4,col=20)



########################################
# Ex2: TURBT
########################################
library(nlme)
library(ggplot2)
# load data
TURBT.data=read.table("TURBT.csv",header=TRUE,sep=',')
dim(TURBT.data) # 285*4
head(TURBT.data)

# process data
start='12/1/2017'
days = as.Date(as.character(TURBT.data$Date.Of.Procedure), format="%m/%d/%Y")-
  as.Date(as.character(start), format="%m/%d/%Y")
days = as.numeric(days) # days since the intervention
doc=as.character(TURBT.data$Performing.Provider)
unique(doc) # 18 physicians
score=TURBT.data$Count_OfElements
after=as.numeric(days>0)

# Q1: compare pre- and post-intervention scores
boxplot(score~after) # not quite right, need to account for subject-level change 
# fit GLS with compound symmetry cov
comsym <- gls(score~after, correlation=corCompSymm(form = ~ 1|doc),  method="REML")
summary(comsym) # interpret beta
# fit LMM with random intercept
LMM1 <- lme (score ~ after, random = ~1 | doc, method='REML') 
summary (LMM1) # Note: we can obtain the corr in gls with the two standard deviation estimates in LMM, the sd in the random effect sigma_b is 0.729
#
fixed.effects(LMM1) # fixed effects coeff 
random.effects(LMM1) # ordered random effects, BLUP (in this case, just b_i)


# Q2: trend of post-intervention scores
data1=data.frame(days,score,doc)
data2=subset(data1,days>0)
# speghetti plot
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  ggtitle("Score change after intervention") +
  theme_classic()
# fit LMM with random intercept
LMM2 = lme(score ~ days,random = ~ 1 | doc,  data = data2)
summary(LMM2)
# fit LMM with random intercept and slope
LMM3 <- lme(score ~ days, random = ~ 1+ days | doc, data = data2)
summary (LMM3) 
vcov(LMM3) # fisher information inverse
random.effects(LMM3) 
# translate this into plain language:  
# who has overall better performance? 
subj1=rownames(random.effects(LMM3))[which.max(random.effects(LMM3)[,1])]
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc==subj1,],aes(x=days, y=score,color=subj1))+
  ggtitle("Score change after intervention")  +
  theme_classic()
# who improve over time
subj2=rownames(random.effects(LMM3))[fixed.effects(LMM3)[2]+random.effects(LMM3)[,2]>0]
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc==subj2,],aes(x=days, y=score,color=subj2))+
  ggtitle("Score change after intervention")  +
  theme_classic()

# Caution: goodness of fit check is also important, but not covered here
# some outlying curves are not well represented in this example:
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc=='John Naitoh MD',],aes(x=days, y=score,color="John Naitoh MD"))+
  ggtitle("Score change after intervention")  +
  theme_classic()
# apparently this guy has increasing trend, but not reflected in the estimated parameters (due to model assumption that the random effect is Gaussian but this one seems to be an outlier) 


######################################
# Ex3: GCase
######################################
# load data
GCase.data=read.table("GCase.csv",header=TRUE,sep=',')
dim(GCase.data) # 1562*6
colnames(GCase.data)
table(GCase.data$DIAGNOSIS)

# process data
GCase.data$visit=as.numeric(as.factor(GCase.data$CLINICAL_EVENT))-1 # BL=0, V04=1, V06=2, V08=3
GCase.data$mutation=as.numeric(GCase.data$GBA.Mutation!='N') # 1=mutation, 0=no


# Q1: GCase ~ time + sex+PD+mutation
# fit a random intercept model, group by patient (PATNO)
LMM1 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation, random = ~1 | PATNO,  data = GCase.data, method='REML') 
summary (LMM1) 
# pay attention to visit fixed effect 0.21. overtime the average increase in GCase activity is 0.21 per unit change in time, taking all measure and all subject into account.  On average, male has higher CGase activity than female and the difference is 0.697. If a patient with mutation, then this patient has lower GCase activity. 


# Q2: PD rate of change =? Control rate of change
# fit subgroup analysis, separate LMM for PD and control
LMM1.PD <- lme (TESTVALUE ~ visit+GENDER+mutation, random = ~1 | PATNO,  data = GCase.data, subset=DIAGNOSIS=="PD", method='REML') 
summary (LMM1.PD)# rate of change: 0.24
LMM1.Control <- lme (TESTVALUE ~ visit+GENDER+mutation, random = ~1 | PATNO,  data = GCase.data, subset=DIAGNOSIS=="Control", method='REML') 
summary (LMM1.Control)# rate of change: 0.13 (not significant)

# fit one analysis, with interaction
LMM2.1 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation+DIAGNOSIS*visit, random = ~1 | PATNO,  data = GCase.data, method='ML') # do NOT use REML for likelihood ratio b/c we need LRT to compare models. If we use "REML" we do not have likelihood
LMM2.2 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation, random = ~1 | PATNO,  data = GCase.data, method='ML')
anova(LMM2.2,LMM2.1)  # LRT of interaction (not significant)
LMM2<- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation+DIAGNOSIS*visit, random = ~1 | PATNO,  data = GCase.data, method='REML') 
summary (LMM2) # Wald test of interaction ( not significant)
## the diff in rate of change is 0.1127. PD has higher change in GCase than ctrl but it is not significant
# 0.13: the rate of change for ctrl samples. 

# discrepency? Why? 
# -- subgroup analysis is different from unified analysis (where effect sizes in other variables 
# are enforced to be the same)
# -- (relevant to this case) the time effect for PD group is still significant in the unified 
# model
beta_PD=fixed.effects(LMM2)[2]+fixed.effects(LMM2)[6] # beta_visit+beta_visit*diagPD
beta_PD_std=sqrt(vcov(LMM2)[2,2]+vcov(LMM2)[6,6]+2*vcov(LMM2)[2,6]) # std(beta_visit+beta_visit*diagPD)
1-pnorm(beta_PD/beta_PD_std) # approx wald p-value

