mtop = glm(formula = trt.effect ~ lag + pre.richness + prcnt_LZ +
             offset(preEWM) + YSI, 
           family = binomial(link = "logit"), data = pre)

summary(mtop)

mtop2 = glm(formula = nab.effect ~ lag + pre.richness + prcnt_LZ +
              + YSI, 
            family = binomial(link = "logit"), data = pre)

summary(mtop2)

# top models for each response variable
# NOTE: swapping p.area for prcnt_LZ. We can use p.area (% permitted, treated)
# but the swap would mean we lose a lot of data 
# NOTE: remove the active_ing from the models. 
# NOTE: try to return the models to the form that Brooks had these in once I fix my station numbers. 

# here are two different version 

# 8.11 notes from Jake
# don't worry about YSI because we might not have the same data to inform this.
# consider looking at specific species and how they are affected by certain variables. 
# consider looking at APM and working it into model
# could consider water quality parameters, but since George didn't use them, it may not be a priority. 
# calculate the proportion of lake treated and look at LZ and lake percentages. 


#NOTE: avoid highly correlated variables in the same model, but in the case of lag here, 
# the lag:YSI is the relationship between the variables which is different than the variables themselves. 

# NOTE: add some variable that is related to consecutive treatment years

# treatment effect models -----
mtop = glm(formula = trt.effect ~ prcnt_LZ + lag + YSI +  
             lag:YSI + offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

mtop = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + lag +
             offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

mtop = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + YSI +
             offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

mtop = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + YSI + lag +
             offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

mtop = glm(formula = trt.effect ~ pre.richness + prcnt_LZ +
             offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

mtop = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + depth_ft +
             offset(preEWM), 
           family = binomial(link = "logit"), data = pre)

summary(mtop)

# native abundance model -----
mtop2 = glm(
  formula = nab.effect ~ preEWM * prcnt_LZ + offset(pre.nab), 
  family = binomial(link = "logit"), 
  data = pre
)

# community changes model -----
mtop3 = glm(formula = com.effect ~ lag * prcnt_LZ + active_ing, family = binomial(link = "logit"), data = pre)