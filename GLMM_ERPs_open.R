# This script fits General Linear Mixed Models to trial-level ERPs data
# DOI of publication:
# Authors: Eliana Nicolaisen-Sobesky; Dr. Álvaro Cabana

#------------------------------------------------------------------
#Load libraries 
#------------------------------------------------------------------

library(R.matlab)
library(lmerTest)
library(reshape2)
library(car)
library(R.utils)
library(emmeans)

#------------------------------------------------------------------
# LOAD ERPs DATA
#------------------------------------------------------------------

#Load ERPs data (mat matrix)
# mat is a matrix of x rows and y columns
# x has dimension number_of_trials*number_of_subjects
# y (columns) include: trial_number, task-condition (fair, mid-value, unfair), subject_ID, group (control, MD/SA) and voltage.
# The column voltage is the ERP voltage value for each trial for each subject (corresponding to the average over time window and electrodes of interest).

file_MFN = './path/to/ERPs/data/MFN'
mat_MFN  = read.csv(file_MFN)

file_LPP = './path/to/ERPs/data/LPP'
mat_LPP  = read.csv(file_LPP)

#------------------------------------------------------------------
# GENERAL LINEAR MODELS
#------------------------------------------------------------------

#FIT AND TEST MODELS
# Fit MFN models
mFull_MFN = lmer(voltage~group*condition*trial+(1|subject_ID),data=mat_MFN)
m2_MFN    = lmer(voltage~group*condition+trial+(1|subject_ID),data=mat_MFN)
# Test MFN models
AIC(mFull_MFN)
AIC(m2_MFN)

# Fit P3/LPP models
mFull_LPP = lmer(voltage~group*condition*trial+(1|subject_ID),data=mat_LPP)
m2_LPP    = lmer(voltage~group*condition+trial+(1|subject_ID),data=mat_LPP)
# Test P3/LPP models
AIC(mFull_LPP)
AIC(m2_LPP)

summary(mFull_MFN)
summary(m2_LPP)

#------------------------------------------------------------------
# STATISTICAL TESTS
#------------------------------------------------------------------

Anova(mFull_MFN,type=3)
Anova(m2_LPP,type=3)

#------------------------------------------------------------------
# POST-HOCS
#------------------------------------------------------------------

# POST-HOCs MFN

emm_options(pbkrtest.limit=15000)
mFull.comp1             = emmeans(mFull_MFN,"condition")
mFull.comp2             = emmeans(mFull_MFN,"group", by="condition")
comp1_MFN.cond          = pairs(mFull.comp1)
comp2_MFN.cond_by_group = pairs(mFull.comp2)
summary(comp1_MFN.cond,by=NULL, adjust='bonferroni')
summary(comp2_MFN.cond_by_group,by=NULL, adjust='bonferroni')

#POST-HOC SLOPES MFN

# Test vs zero 
mFull_MFN_slopes = lmer(voltage~trial:condition:group+group*condition+(1|subject_ID),data=mat_MFN)
summary(mFull_MFN_slopes)

# Test slopes between groups
mFull_MFN_groups = lmer(voltage~trial:condition+trial:condition:group+group*condition+(1|subject_ID),data=mat_MFN)
summary(mFull_MFN_groups)

# Test slops by condition within groups
mFull_MFN_cond_by_groups = lmer(voltage~group+group:condition+group:trial:condition+trial:group+(1|subject_ID),data=mat_MFN)
summary(mFull_MFN_cond_by_groups)

# POST-HOCs LPP

emm_options(pbkrtest.limit=15000)
m2.comp1       = emmeans(m2_LPP,"condition")
comp1_LPP.cond = pairs(m2.comp1)
summary(comp1_LPP.cond,by=NULL, adjust='bonferroni')
