# Seroy and Grunbaum 2020

# All statistical analyses described in manuscript

# --------------------------------------------------------------------------------------------
# FIGURE 1
# DOES DENSITY + SETTLEMENT TIME CONTROL COLONY SIZE IN THE FIELD?
# --------------------------------------------------------------------------------------------

# load in data
fig1_data = read.csv(file = 'data/field_data.csv', header = T)
head(fig1_data)

# check that plates are factors
levels(fig1_data$plate)

# preliminary stats
mean(fig1_data$field_area)
sd(fig1_data$field_area)
min(fig1_data$field_area)
max(fig1_data$field_area)
length(fig1_data$field_area)


# Create linear models for data

# what do the data look like?
hist(fig1_data$field_area) # skewed
hist(log(fig1_data$field_area)) # looks better log-transformed

# Linear model with log-transformed areas in the field
# using density as continuous, not factor
fig1_mod0 = lm(log(field_area) ~ 1, data = fig1_data)
fig1_mod1 = lm(log(field_area) ~ density, data = fig1_data)
fig1_mod2 = lm(log(field_area) ~ settle_time, data = fig1_data)
fig1_mod3 = lm(log(field_area) ~ density + settle_time, data = fig1_data)
fig1_mod4 = lm(log(field_area) ~ density * settle_time, data = fig1_data)


# select the best model using AIC
AIC(fig1_mod0) # 662.4108
AIC(fig1_mod1) # 645.7763
AIC(fig1_mod2) # 540.0819
AIC(fig1_mod3) # 520.7702 - best model
AIC(fig1_mod4) # 522.7639

# check model assumptions
plot(fig1_mod3)

# model summary
summary(fig1_mod3) 
anova(fig1_mod3)



# --------------------------------------------------------------------------------------------
# FIGURE 2
# DOES THE MODEL QUANTITATIVLY PREDICT FIELD CONDTIONS?
# --------------------------------------------------------------------------------------------

# uses same data set at Fig 1 (fig1_data)

# create linear models to describe the relationship btwn model and field areas

# what do the data look like?
hist(log(fig1_data$mod_perc))
hist(log(fig1_data$field_perc))

# log transformed, with and without forced intercept at 0
fig2_mod1 = lm(log(field_perc)~log(mod_perc)-1, fig1_data)
fig2_mod2 = lm(log(field_perc)~log(mod_perc), fig1_data)

AIC(fig2_mod1) # 433.7741 - even though this model is not lowest, going with this one because no intercept makes more biological sense
AIC(fig2_mod2) # 426.9539

# check model assumptions
plot(fig2_mod1)

# model summary
summary(fig2_mod1)



# --------------------------------------------------------------------------------------------
# FIGURE 4
# DOES DENSITY + pH AFFECT COST OF DEFENSE?
# --------------------------------------------------------------------------------------------

# load in data
fig4_data = read.csv(file = 'data/Figure4_data.csv', header=T)
head(fig4_data)
dim(fig4_data)


# calculate normalized and  relative metrics
fig4_data$norm_COD = fig4_data$COD/5000         # calculate NORMALIZED COD for each run
fig4_data$norm_def = fig4_data$mean_pred/5000   # calculate NORMALIZED def. col size for each run
fig4_data$rel_COD = fig4_data$COD/fig4_data$mean_cont       # calculate RELATIVE COD for each run
fig4_data$rel_def = fig4_data$mean_pred/fig4_data$mean_cont # calculate RELATIVE def. col size for each run

head(fig4_data)


# make pH a factor
fig4_data$pH = as.factor(fig4_data$pH)

# reorder pH values to put in analysis so that everything is compared to pH 7.9
fig4_data$fact_pH = as.character(fig4_data$pH)
fig4_data$fact_pH[fig4_data$fact_pH == 7.9] <- "a7.9"
fig4_data$fact_pH[fig4_data$fact_pH == 7.6] <- "b7.6"
fig4_data$fact_pH[fig4_data$fact_pH == 7.3] <- "c7.3"

head(fig4_data)


# preliminary stats

# general stats - normalized COD
mean(fig4_data$norm_COD)
sd(fig4_data$norm_COD)
length(fig4_data$norm_COD)

# general stats - normalized def. col size
mean(fig4_data$norm_def)
sd(fig4_data$norm_def)
length(fig4_data$norm_def)

# general stats - rel COD
mean(fig4_data$rel_COD)
sd(fig4_data$rel_COD)
length(fig4_data$rel_COD)

# general stats - normalized def. col size
mean(fig4_data$rel_def)
sd(fig4_data$rel_def)
length(fig4_data$rel_def)


# Linear models for Fig 4a: Normalized COD vs. density + pH
# -------------------------------------------------------------------

# what do the data look like? ...transform?
hist(fig4_data$norm_COD)

# GLM with gamma distribution and log link
fig4_mod0a = glm(norm_COD ~ 1, family = Gamma(link='log'), data = fig4_data)             
fig4_mod1a = glm(norm_COD ~ pop_dens, family = Gamma(link='log'), data = fig4_data)      
fig4_mod2a = glm(norm_COD ~ fact_pH, family = Gamma(link='log'), data = fig4_data)            
fig4_mod3a = glm(norm_COD ~ pop_dens + fact_pH, family = Gamma(link='log'), data = fig4_data) 
fig4_mod4a = glm(norm_COD ~ pop_dens * fact_pH, family = Gamma(link='log'), data = fig4_data) 

# use AIC to select best model
AIC(fig4_mod0a) # -5677.473
AIC(fig4_mod1a) # -6263.011
AIC(fig4_mod2a) # -5772.982
AIC(fig4_mod3a) # -6464.192 - best model, most parsimonious of best two
AIC(fig4_mod4a) # -6462.685

# check model assumptions
plot(fig4_mod3a)

# model summary
summary(fig4_mod3a)



# Linear models for Fig 4b: Normalized def colony size vs. density + pH
# -------------------------------------------------------------------

# what does the data look like? ... transform?
hist(fig4_data$norm_def)

# GLM with gamma distribution and log link 
fig4_mod0b = glm(norm_def ~ 1, family = Gamma(link='log'), data = fig4_data)             
fig4_mod1b = glm(norm_def ~ pop_dens, family = Gamma(link='log'), data = fig4_data)      
fig4_mod2b = glm(norm_def ~ fact_pH, family = Gamma(link='log'), data = fig4_data)            
fig4_mod3b = glm(norm_def ~ pop_dens + fact_pH, family = Gamma(link='log'), data = fig4_data) 
fig4_mod4b = glm(norm_def ~ pop_dens * fact_pH, family = Gamma(link='log'), data = fig4_data) 

# use AIC to select best model
AIC(fig4_mod0b) # -6501.198
AIC(fig4_mod1b) # -7618.766
AIC(fig4_mod2b) # -6567.982  
AIC(fig4_mod3b) # -7878.141 best model, most parsimonious of best two
AIC(fig4_mod4b) # -7875.437

# check model assumptions
plot(fig4_mod3b) 

# model summary
summary(fig4_mod3b)


# Linear models for Fig 4c: Relative COD vs. density + pH
# -------------------------------------------------------------------

# what do the data look like? looks pretty normal, no transform
hist(fig4_data$rel_COD)

# Regular linear model
fig4_mod0c = lm(rel_COD ~ 1, data = fig4_data)             
fig4_mod1c = lm(rel_COD ~ pop_dens, data = fig4_data)      
fig4_mod2c = lm(rel_COD ~ fact_pH, data = fig4_data)            
fig4_mod3c = lm(rel_COD ~ pop_dens + fact_pH, data = fig4_data) 
fig4_mod4c = lm(rel_COD ~ pop_dens * fact_pH, data = fig4_data) 

# use AIC to select best model
AIC(fig4_mod0c) # -1029.101
AIC(fig4_mod1c) # -1057.987
AIC(fig4_mod2c) # -1276.75
AIC(fig4_mod3c) # -1315.831 best model, most parsimonious of best two
AIC(fig4_mod4c) # -1314.105

# check model assumptions
plot(fig4_mod3c) 

# model summary
summary(fig4_mod3c)


# Linear models for Fig 4d: Relative def colony size vs. density + pH
# -------------------------------------------------------------------

# what does the data look like? looks pretty normal, no transform
hist(fig4_data$rel_def)

# Regular linear model
fig4_mod0d = lm(rel_def ~ 1, data = fig4_data)             
fig4_mod1d = lm(rel_def ~ pop_dens, data = fig4_data)      
fig4_mod2d = lm(rel_def ~ fact_pH, data = fig4_data)            
fig4_mod3d = lm(rel_def ~ pop_dens + fact_pH, data = fig4_data) 
fig4_mod4d = lm(rel_def ~ pop_dens * fact_pH, data = fig4_data) 

# use AIC to select best model
AIC(fig4_mod0d) # -1029.101
AIC(fig4_mod1d) # -1057.987
AIC(fig4_mod2d) # -1276.75  
AIC(fig4_mod3d) # -1315.831 best model, most parsimonious of best two
AIC(fig4_mod4d) # -1314.105

# check model assumptions
plot(fig4_mod3d)

# model summary
summary(fig4_mod3d)



# --------------------------------------------------------------------------------------------
# FIGURE 5
# DOES DENSITY, SETTLEMENT TIME + pH AFFECT COST OF DEFENSE?
# --------------------------------------------------------------------------------------------

# load in data 
fig5_data = read.csv(file = "data/Figure5_data.csv", header = T) # all data
head(fig5_data)
dim(fig5_data)


# calculate normalized and relative metrics
fig5_data$norm_COD = fig5_data$COD_late/5000 # add normalized COD column
fig5_data$norm_def = fig5_data$mean_pred_late/5000 # add normalized defended colony size column
fig5_data$rel_COD = fig5_data$COD/fig5_data$mean_cont # add relative COD column
fig5_data$rel_def = fig5_data$mean_pred/fig5_data$mean_cont # add relative defended colony size column


# preliminary statistics
# mean, min, max and sd for normalized COD
mean(fig5_data$norm_COD)
sd(fig5_data$norm_COD)
length(fig5_data$norm_COD)

# mean, min, max and sd for normalized defended colony size
mean(fig5_data$norm_def)
sd(fig5_data$norm_def)
length(fig5_data$norm_def)

# mean, min, max and sd for relative COD
mean(fig5_data$rel_COD)
sd(fig5_data$rel_COD)
length(fig5_data$rel_COD)

# mean, min, max and sd for relative defended colony size
mean(fig5_data$rel_def)
sd(fig5_data$rel_def)
length(fig5_data$rel_def)

# make pH a factor
fig5_data$pH = as.factor(fig5_data$pH)

# reorder pH values to put in analysis so that everything is compared to pH 7.9
fig5_data$fact_pH = as.character(fig5_data$pH)
fig5_data$fact_pH[fig5_data$fact_pH == 7.9] <- "a7.9"
fig5_data$fact_pH[fig5_data$fact_pH == 7.6] <- "b7.6"
fig5_data$fact_pH[fig5_data$fact_pH == 7.3] <- "c7.3"

head(fig5_data)
tail(fig5_data)


# Create linear models for Fig 5a: NORMALIZED COD vs. density, settlement time + pH
# -------------------------------------------------------------------

# model of the mean
fig5_mod0a = lm(norm_COD ~ 1, fig5_data) 

# single predictor - no interaction
fig5_mod1a = lm(norm_COD ~ fact_pH, fig5_data)
fig5_mod2a = lm(norm_COD ~ pop_dens, fig5_data)
fig5_mod3a = lm(norm_COD ~ settle_time, fig5_data)

# two predictors - no interaction
fig5_mod4a = lm(norm_COD ~ pop_dens + fact_pH, fig5_data)
fig5_mod5a = lm(norm_COD ~ settle_time + fact_pH, fig5_data)
fig5_mod6a = lm(norm_COD ~ pop_dens + settle_time, fig5_data)

# three predictors - no interaction
fig5_mod7a = lm(norm_COD ~ pop_dens + settle_time + fact_pH, fig5_data)

# interaction models - two way interactions
fig5_mod8a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens, fig5_data)
fig5_mod9a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time, fig5_data)
fig5_mod10a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + pop_dens:settle_time, fig5_data)

fig5_mod11a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time, fig5_data)
fig5_mod12a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + pop_dens:settle_time, fig5_data)
fig5_mod13a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time + pop_dens:settle_time, fig5_data)
fig5_mod14a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time, fig5_data)

# with a three-way interaction
fig5_mod15a = lm(norm_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time + fact_pH:pop_dens:settle_time , fig5_data)

# Use AIC to determine best model
AIC(fig5_mod0a) # -10490.07
AIC(fig5_mod1a) # -10601.08
AIC(fig5_mod2a) # -10947.93
AIC(fig5_mod3a) # -10765.59
AIC(fig5_mod4a) # -11093.84
AIC(fig5_mod5a) # -10896.5
AIC(fig5_mod6a) # -11312.98
AIC(fig5_mod7a) # -11494.6 

AIC(fig5_mod8a) # -11513.67
AIC(fig5_mod9a) # -11504.54
AIC(fig5_mod10a) # -11548.13
AIC(fig5_mod11a) # -11523.79
AIC(fig5_mod12a) # -11567.93
AIC(fig5_mod13a) # -11558.51
AIC(fig5_mod14a) # -11578.5 **best model**

AIC(fig5_mod15a) # -11576.53

# check model assumptions
plot(fig5_mod14a)

# model summary
summary(fig5_mod14a)
anova(fig5_mod14a)



# Create linear models for Fig 5b: NORMALIZED DEF COL SIZE vs. density, settlement time + pH
# ------------------------------------------------------------------- 

# model of the mean
fig5_mod0b = lm(norm_def ~ 1, fig5_data) 

# single predictor - no interaction
fig5_mod1b = lm(norm_def ~ fact_pH, fig5_data)
fig5_mod2b = lm(norm_def ~ pop_dens, fig5_data)
fig5_mod3b = lm(norm_def ~ settle_time, fig5_data)

# two predictors - no interaction
fig5_mod4b = lm(norm_def ~ pop_dens + fact_pH, fig5_data)
fig5_mod5b = lm(norm_def ~ settle_time + fact_pH, fig5_data)
fig5_mod6b = lm(norm_def ~ pop_dens + settle_time, fig5_data)

# three predictors - no interaction
fig5_mod7b = lm(norm_def ~ pop_dens + settle_time + fact_pH, fig5_data)

# interaction models
fig5_mod8b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens, fig5_data)
fig5_mod9b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time, fig5_data)
fig5_mod10b = lm(norm_def ~ pop_dens + settle_time + fact_pH + pop_dens:settle_time, fig5_data)

fig5_mod11b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time, fig5_data)
fig5_mod12b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + pop_dens:settle_time, fig5_data)
fig5_mod13b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time + pop_dens:settle_time, fig5_data)
fig5_mod14b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time, fig5_data)

fig5_mod15b = lm(norm_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time + pH:pop_dens:settle_time , fig5_data)



AIC(fig5_mod0b) # -12108.64
AIC(fig5_mod1b) # -12188.16
AIC(fig5_mod2b) # -13158.89
AIC(fig5_mod3b) # -12417.37
AIC(fig5_mod4b) # -13307.62
AIC(fig5_mod5b) # -12513.08
AIC(fig5_mod6b) # -13759.43
AIC(fig5_mod7b) # -13972.69 

AIC(fig5_mod8b) # -13990.77
AIC(fig5_mod9b) # -13982
AIC(fig5_mod10b) # -14021.76
AIC(fig5_mod11b) # -14000.25
AIC(fig5_mod12b) # -14040.49 
AIC(fig5_mod13b) # -14031.46
AIC(fig5_mod14b) # -14050.35 - best model

AIC(fig5_mod15b) # -14051.63 - lowest but within 2 of a more parsimonious model

# check model assumptions
plot(fig5_mod14b)

# model summary
summary(fig5_mod14b)
anova(fig5_mod14b)


# Create linear models for Fig 5c: RELATIVE COD vs. density, settlement time + pH
# ------------------------------------------------------------------- 

# model of the mean
fig5_mod0c = lm(rel_COD ~ 1, fig5_data) 

# single predictor - no interaction
fig5_mod1c = lm(rel_COD ~ fact_pH, fig5_data)
fig5_mod2c = lm(rel_COD ~ pop_dens, fig5_data)
fig5_mod3c = lm(rel_COD ~ settle_time, fig5_data)

# two predictors - no interaction
fig5_mod4c = lm(rel_COD ~ pop_dens + fact_pH, fig5_data)
fig5_mod5c = lm(rel_COD ~ settle_time + fact_pH, fig5_data)
fig5_mod6c = lm(rel_COD ~ pop_dens + settle_time, fig5_data)

# three predictors - no interaction
fig5_mod7c = lm(rel_COD ~ pop_dens + settle_time + fact_pH, fig5_data)

# interaction models - two way interactions
fig5_mod8c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens, fig5_data)
fig5_mod9c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time, fig5_data)
fig5_mod10c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + pop_dens:settle_time, fig5_data)

fig5_mod11c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time, fig5_data)
fig5_mod12c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + pop_dens:settle_time, fig5_data)
fig5_mod13c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time + pop_dens:settle_time, fig5_data)
fig5_mod14c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time, fig5_data)

# with a three-way interaction
fig5_mod15c = lm(rel_COD ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time + pH:pop_dens:settle_time , fig5_data)

AIC(fig5_mod0c) # -349.3991
AIC(fig5_mod1c) # -575.5623
AIC(fig5_mod2c) # -356.6021
AIC(fig5_mod3c) # -419.152
AIC(fig5_mod4c) # -584.0243
AIC(fig5_mod5c) # -655.3281
AIC(fig5_mod6c) # -426.7303
AIC(fig5_mod7c) # -664.2777 ** best model **

AIC(fig5_mod8c) # -661.7561
AIC(fig5_mod9c) # -660.3495
AIC(fig5_mod10c) # -662.5852 
AIC(fig5_mod11c) # -657.8279
AIC(fig5_mod12c) # -660.0638
AIC(fig5_mod13c) # -658.657
AIC(fig5_mod14c) # -656.1356 

AIC(fig5_mod15c) # -654.0526

# check model assumptions
plot(fig5_mod7c) 


# model summary
summary(fig5_mod7c)
anova(fig5_mod7c)


# Create linear models for Fig 5d: RELATIVE DEF COL SIZES vs. density, settlement time + pH
# ------------------------------------------------------------------- 

# model of the mean
fig5_mod0d = lm(rel_def ~ 1, fig5_data) 

# single predictor - no interaction
fig5_mod1d = lm(rel_def ~ fact_pH, fig5_data)
fig5_mod2d = lm(rel_def ~ pop_dens, fig5_data)
fig5_mod3d = lm(rel_def ~ settle_time, fig5_data)

# two predictors - no interaction
fig5_mod4d = lm(rel_def ~ pop_dens + fact_pH, fig5_data)
fig5_mod5d = lm(rel_def ~ settle_time + fact_pH, fig5_data)
fig5_mod6d = lm(rel_def ~ pop_dens + settle_time, fig5_data)

# three predictors - no interaction
fig5_mod7d = lm(rel_def ~ pop_dens + settle_time + fact_pH, fig5_data)

# interaction models - two way interactions
fig5_mod8d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens, fig5_data)
fig5_mod9d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time, fig5_data)
fig5_mod10d = lm(rel_def ~ pop_dens + settle_time + fact_pH + pop_dens:settle_time, fig5_data)

fig5_mod11d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time, fig5_data)
fig5_mod12d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + pop_dens:settle_time, fig5_data)
fig5_mod13d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:settle_time + pop_dens:settle_time, fig5_data)
fig5_mod14d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time, fig5_data)

# with a three-way interaction
fig5_mod15d = lm(rel_def ~ pop_dens + settle_time + fact_pH + fact_pH:pop_dens + fact_pH:settle_time + pop_dens:settle_time + pH:pop_dens:settle_time , fig5_data)


AIC(fig5_mod0d) # -349.3991
AIC(fig5_mod1d) # -575.5623
AIC(fig5_mod2d) # -356.6021
AIC(fig5_mod3d) # -419.152
AIC(fig5_mod4d) # -584.0243
AIC(fig5_mod5d) # -655.3281
AIC(fig5_mod6d) # -426.7303
AIC(fig5_mod7d) # -664.2777 ** best model **

AIC(fig5_mod8d) # -661.7561
AIC(fig5_mod9d) # -660.3495
AIC(fig5_mod10d) # -662.5852
AIC(fig5_mod11d) # -657.8279
AIC(fig5_mod12d) # -660.0638 
AIC(fig5_mod13d) # -658.657
AIC(fig5_mod14d) # -656.1356

AIC(fig5_mod15d) # -654.0526

# check model assumptions
plot(fig5_mod7d)

# model summary 
summary(fig5_mod7d)
anova(fig5_mod7d)

