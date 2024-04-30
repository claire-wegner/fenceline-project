# Install and load the following packages
library (tidyverse)
library (nlme)
library (lmtest)

# Import the IQ_full dataset
IQ_full <- read.csv("data/IQ_full.csv")

# FUNCTIONS ----
# The following functions need to be run prior to running the analyses and plots

--------------------------------------------------------------------------------
# This function will calculate the difference in behaviour from baseline
# Specify either steps or lying
adjust2BL <- function(data1, var1="steps") {
  if (var1=="steps") {
    data2 <- data1 %>%
      group_by(id) %>%
      mutate(steps_diff = daily_steps-daily_steps[obs_day == "BL"])
    print("New variable steps_diff created")
  }
  if (var1 == "lying") {
    data2 <- data1 %>%
      group_by(id) %>%
      mutate(lying_diff = daily_L_m-daily_L_m[obs_day == "BL"])
    print("New variable lying_diff created")
  }
  return(data2)
}

--------------------------------------------------------------------------------
# This function will plot mean values as a line
geom_hpline <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHpline,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}
  
#' @rdname geom_hpline
#' @format NULL
#' @usage NULL
#' @export
GeomHpline <- ggproto("GeomHpline", GeomSegment,
                      required_aes = c("x", "y"),
                      non_missing_aes = c("size", "colour", "linetype", "width"),
                      default_aes = aes(
                        width = 0.5, colour = "black", size = 2, linetype = 1,
                        alpha = NA
                      ),
                      
                      draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                                            lineend = "butt", linejoin = "round", na.rm = FALSE) {
                        data <- mutate(data, x = x - width/2, xend = x + width, yend = y)
                        ggproto_parent(GeomSegment, self)$draw_panel(
                          data, panel_params, coord, arrow = arrow, arrow.fill = arrow.fill,
                          lineend = lineend, linejoin = linejoin, na.rm = na.rm
                        )
                      }
)

--------------------------------------------------------------------------------
# Function to compute standard errors
# Extract the estimates for SE to use in the creation of error bars
get.SE <- function(x, mod, age=4) {
  #mod is the fitted model in nlme
  if (age==4) y_calf_indicator_vector <- age4.vector(x)
  if (age==6) y_calf_indicator_vector <- age6.vector(x)
  est <- crossprod(y_calf_indicator_vector, summary(mod)$coef$fixed)
  var1 <- t(y_calf_indicator_vector)%*%vcov(mod)%*%y_calf_indicator_vector
  SE <- sqrt(var1)
  return(c(est, est-sqrt(var1), est+sqrt(var1), SE))
}

--------------------------------------------------------------------------------

# ANALYSES ----

# Correct variables
IQ_full$trt <- as.factor(IQ_full$trt)
IQ_full$id <- as.factor(IQ_full$id)

# Create new variables, obs_day2 and parity2, which are the numeric versions of obs_day and parity
IQ_full$obs_day2 <- as.numeric(IQ_full$obs_day)
IQ_full$parity2 <- as.numeric(IQ_full$parity)

# Ensure that the original obs_day is ordered correctly
IQ_full$obs_day <- factor(IQ_full$obs_day, 
                          levels = c("BL", "1", "2", "3", "4", "5", "6", "7", 
                                     "8", "9", "10", "11"))

## 1) Feed-seeking behaviour ---------------------------------------------------

# A mixed linear model with a quadratic relationship between obs_day2 and F_min_hourly
lme.feed2 <- lme(F_min_hourly ~ trt + obs_day2 + I(obs_day2^2)+ trt:obs_day2 + trt:I(obs_day2^2) + percent_near,
                 random = ~ 1 | id,
                 correlation = corAR1(form = ~ 1 | id),
                 data = IQ_full,
                 na.action=na.omit)

# The same model but with no interaction term
lme.feed3 <- lme(F_min_hourly ~ trt + obs_day2 + I(obs_day2^2) + percent_near,
                 random = ~ 1 | id,
                 correlation = corAR1(form = ~ 1 | id),
                 data = IQ_full,
                 na.action=na.omit)

# Use a likelihood ratio test (LRT) to see if the model with the interaction term is a better fit than that without
lrtest(lme.feed2, lme.feed3)
# No, so we can remove the interaction
# The best model is therefore lme.feed2

# Check residuals for final model
qqnorm(residuals(lme.feed3))
qqline(residuals(lme.feed3))
plot(lme.feed3)

# Output for the final model
summary(lme.feed3)
anova(lme.feed3)

# Intra-class correlation coefficient (ICC)
# This is calculated by first extracting the variance of random effects and residual variance from the model output,
# then dividing the variance of random effects by the sum of the two variances.
calf_feeding_ICC <- ((1.438605^2)/((1.438605^2)+(3.953969^2)))
print(calf_feeding_ICC)

# We can also see what the model looks like when we remove the twins
feed.twins.remov <- lme(F_min_hourly ~ trt + obs_day2 + I(obs_day2^2) + percent_near,
                        random = ~ 1 | id,
                        correlation = corAR1(form = ~ 1 | id),
                        data = IQ_full %>% dplyr::filter(single_twins == 'single'),
                        na.action=na.omit)

summary(feed.twins.remov)
anova(feed.twins.remov)

## 2) Step count ---------------------------------------------------------------

IQ_full <- adjust2BL(IQ_full, var1="steps")# Calculate step difference for each day compared to baseline
IQ_full <- IQ_full %>% #Make positive by adding an integer based on the smallest numbers for cows and calves
  mutate(steps_diff2 = ifelse(cow_calf == 'cow', steps_diff + 1506, steps_diff + 1777)) %>%
  mutate(steps_diff_fourthroot = (steps_diff2^0.25)) #Box-cox transformation suggested a fourth-root transformation was best

# COWS

# A mixed linear model with a quadratic relationship between obs_day2 and steps_diff_fourthroot
quad.cow.steps <- lme(steps_diff_fourthroot ~ trt + trt*obs_day2 + trt*I(obs_day2^2) + percent_near + parity2,
                      random = list(id = ~ 1),
                      correlation = corAR1(form = ~1 | id),
                      data = IQ_full %>% dplyr::filter(cow_calf == 'cow'),
                      na.action=na.omit)

# The same model but with no interaction term
poly.cow.2 <- lme(steps_diff_fourthroot ~ trt + obs_day2 + I(obs_day2^2) + percent_near + parity2,
                  random = list(id = ~ 1),
                  correlation = corAR1(form = ~1 | id),
                  data = IQ_full %>% dplyr::filter(cow_calf == 'cow'),
                  na.action=na.omit)

# Use a LRT to see if interaction effect should be included in the model
lrtest(quad.cow.steps, poly.cow.2) 
# Since the test is not significant, we know that the model containing interactions does not explain more of the variation
# We will thus use the model that does not contain an interaction between trt and obs_day2 terms

# Check the residuals of the final model
qqnorm(residuals(poly.cow.2))
qqline(residuals(poly.cow.2))
plot(poly.cow.2)

# Results from the final model
summary(poly.cow.2)
anova(poly.cow.2)

# Intra-class correlation coefficient (ICC)
cow_steps_ICC <- ((0.3121076^2)/((0.3121076^2)+(0.9767383^2)))
print(cow_steps_ICC)

# CALVES

# A mixed linear model with a quadratic relationship between obs_day2 and steps_diff_fourthroot
poly.calf <- lme(steps_diff_fourthroot ~ trt  + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near,
                 random = ~ 1 | id,
                 correlation = corAR1(form = ~ 1 | id),
                 data = IQ_full %>% dplyr::filter(cow_calf == 'calf'),
                 na.action=na.omit)

# The same model but with no interaction term
poly.calf2 <- lme(steps_diff_fourthroot ~ trt  + obs_day2 + I(obs_day2^2) + percent_near,
                  random = ~ 1 | id,
                  correlation = corAR1(form = ~ 1 | id),
                  data = IQ_full %>% dplyr::filter(cow_calf == 'calf'),
                  na.action=na.omit) 

# Use a LRT to see if interaction effect should be included in the model
lrtest(poly.calf, poly.calf2)
# Since the test is not significant, we know the model that without interaction terms explains more of the variance

# Check the residuals
qqnorm(residuals(poly.calf2))
qqline(residuals(poly.calf2))
plot(poly.calf2)

# Results from the final model
summary(poly.calf2)
anova(poly.calf2) 

# Intra-class correlation coefficient (ICC)
calf_steps_ICC <- ((0.3903269^2)/((0.3903269^2)+(0.9246458^2)))
print(calf_steps_ICC)

# Test the model without the twins
steps.twins.remov <- lme(steps_diff_fourthroot ~ trt  + obs_day2 + I(obs_day2^2) + percent_near,
                         random = ~ 1 | id,
                         correlation = corAR1(form = ~ 1 | id),
                         data = IQ_full %>% dplyr::filter(cow_calf == 'calf' & single_twins == 'single'),
                         na.action=na.omit) 

summary(steps.twins.remov)
anova(steps.twins.remov)

## 3) Lying time ---------------------------------------------------------------

IQ_full <- adjust2BL(IQ_full, var1="lying") # Calculate difference in lying time compared to baseline

# COWS

# A mixed linear model with a quadratic relationship between obs_day2 and lying_diff
lme.cow.LT2 <- lme(lying_diff ~ trt + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near + parity2,
                   random = ~ 1 | id,
                   correlation = corAR1(form = ~ 1 | id),
                   data = IQ_full %>% dplyr::filter(cow_calf == 'cow'),
                   na.action=na.omit)

# The same model but with no interaction term
lme.cow.LT4 <- lme((lying_diff) ~ trt + obs_day2 + I(obs_day2^2) + percent_near + parity2,
                   random = ~ 1 | id,
                   correlation = corAR1(form = ~ 1 | id),
                   data = IQ_full %>% dplyr::filter(cow_calf == 'cow'),
                   na.action=na.omit)

# Use a LRT to see if interaction effect should be included in the model
lrtest(lme.cow.LT2, lme.cow.LT4) 
# Test is significant, so interaction effect needs to stay

# Check the residuals for the final model
qqnorm(residuals(lme.cow.LT2))
qqline(residuals(lme.cow.LT2))
plot(lme.cow.LT2)

# Results from the final model
summary(lme.cow.LT2)
anova(lme.cow.LT2)

cow_lying_ICC <- ((32.57274^2)/((32.57274^2)+(110.0537^2)))
print(cow_lying_ICC)

# CALVES

# A mixed linear model with a quadratic relationship between obs_day2 and lying_diff
lme.calf.LT2 <- lme(lying_diff ~ trt + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near,
                    random = ~ 1 | id,
                    correlation = corAR1(form = ~ 1 | id),
                    data = IQ_full %>% dplyr::filter(cow_calf == 'calf'),
                    na.action=na.omit)

# The same model but with no interaction term
lme.calf.LT3 <- lme(lying_diff ~ trt + obs_day2 + I(obs_day2^2) + percent_near,
                    random = ~ 1 | id,
                    correlation = corAR1(form = ~ 1 | id),
                    data = IQ_full %>% dplyr::filter(cow_calf == 'calf'),
                    na.action=na.omit)

# Use a LRT to see if interaction effect should be included in the model
lrtest(lme.calf.LT2, lme.calf.LT3)
# Test is significant, so interaction effect needs to stay

# Check the residuals for the final model
qqnorm(residuals(lme.calf.LT2))
qqline(residuals(lme.calf.LT2))
plot(lme.calf.LT2)

# Results from the final model
summary(lme.calf.LT2)
anova(lme.calf.LT2)

# Test the model without the twins
lying.twins.remov <- lme(lying_diff ~ trt + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near,
                         random = ~ 1 | id,
                         correlation = corAR1(form = ~ 1 | id),
                         data = IQ_full %>% dplyr::filter(cow_calf == 'calf' & single_twins == 'single'),
                         na.action=na.omit)

summary(lying.twins.remov)
anova(lying.twins.remov)

# DESCRIPTIVE STATISTICS -------------------------------------------------------

# For all descriptive statistics, view stats per treatment by changing "obs_day" for "obs_day + trt"

## 1) Vocalizations ------------------------------------------------------------

# Mean and SD for vocalizations (based on raw data)
mean_voc_day <- aggregate(voc_perc_yes ~ obs_day + trt, data = IQ_full, FUN = mean)
print(mean_voc_day)

sd_voc_day <- aggregate(voc_perc_yes ~ obs_day + trt, data = IQ_full, FUN = sd)
print(sd_voc_day)

## 2) Feed-seeking -------------------------------------------------------------

# Mean and SD time spent on feeding behaviour (based on raw data)
mean_feed_day <- aggregate(F_min_hourly ~ obs_day, data = IQ_full, FUN = mean)
print(mean_feed_day)

sd_feed_day <- aggregate(F_min_hourly ~ obs_day, data = IQ_full, FUN = sd)
print(sd_feed_day)

## 3) Step count ---------------------------------------------------------------

# COWS
# Daily step count (based on raw data)
mean_steps_day_cow <- aggregate(daily_steps ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = mean)
print(mean_steps_day_cow)

sd_steps_day_cow <- aggregate(daily_steps ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = sd)
print(sd_steps_day_cow)

# Difference in daily step count (based on raw data)
mean_diff_steps_cow <- aggregate(steps_diff ~ obs_day2, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = mean)
print(mean_diff_steps_cow)

sd_diff_steps_cow <- aggregate(steps_diff ~ obs_day2, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = sd)
print(sd_diff_steps_cow)

# CALVES
# Daily step count (based on raw data)
mean_steps_day_calf <- aggregate(daily_steps ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = mean)
print(mean_steps_day_calf)

sd_steps_day_calf <- aggregate(daily_steps ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = sd)
print(sd_steps_day_calf)

# Difference in daily step count (based on raw data)
mean_diff_steps_calf <- aggregate(steps_diff ~ obs_day2, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = mean)
print(mean_diff_steps_calf)

sd_diff_steps_calf <- aggregate(steps_diff ~ obs_day2, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = sd)
print(sd_diff_steps_calf)

## 4) Lying time ---------------------------------------------------------------

# COWS
# Daily lying time in h/d (based on raw data)
mean_LT_day_cow <- aggregate((daily_L_m/60) ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = mean)
print(mean_LT_day_cow)

sd_LT_day_cow <- aggregate((daily_L_m/60)  ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = sd)
print(sd_LT_day_cow)

# CALVES
# Daily lying time in h/d (based on raw data)
mean_LT_day <- aggregate((daily_L_m/60) ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = mean)
print(mean_LT_day)

sd_LT_day <- aggregate((daily_L_m/60) ~ obs_day, data = IQ_full %>% dplyr::filter(cow_calf == 'calf'), FUN = sd)
print(sd_LT_day)

## 5) Parity -------------------------------------------------------------------

mean_parity <- aggregate(parity2 ~ obs_day + trt, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = mean)
print(mean_parity)

sd_parity <- aggregate(parity2 ~ obs_day + trt, data = IQ_full %>% dplyr::filter(cow_calf == 'cow'), FUN = sd)
print(sd_parity)

# PLOTS ------------------------------------------------------------------------
## 1) Lying time  ----

### a) Calves (Fib 3A) ---------------------------------------------------------
# Run a model with lying time converted to h/d
# This is purely for the purpose of plotting lying time in h/d
# Models using h/d and min/d were previously computed and deemed identical, save for units
LT.calf.hours <- lme((lying_diff/60) ~ trt + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near,
                     random = ~ 1 | id,
                     correlation = corAR1(form = ~ 1 | id),
                     data = IQ_full %>% dplyr::filter(cow_calf == 'calf'),
                     na.action=na.omit)

# Use results from summary to plot the polynomial
summary(LT.calf.hours)$coef$fixed

# First, we use the coefficients to plot data for the 4_mo group, since it is the reference in the model
# 4 mo 
x <- 1:11
Coef4_calf_LT <- summary(LT.calf.hours)$coef$fixed[c(1,3,4)] #Look at summary of fixed effects to see which belong only to 4_mo group
y4_calf_LT <- Coef4_calf_LT[1]+Coef4_calf_LT[2]*x+Coef4_calf_LT[3]*x^2
plot(x, y4_calf_LT)

# Then, we make the 6_mo plot by adding the effects of the 4_mo group to the fixed effects that only pertain to the 6_mo group
# (i.e., interactions of 6_mo and obs_day)
# 6 mo
x <- 1:11
Coef6_calf_LT <- summary(LT.calf.hours)$coef$fixed[c(2,6,7)] + Coef4_calf_LT
y6_calf_LT <- Coef6_calf_LT[1]+Coef6_calf_LT[2]*x+Coef6_calf_LT[3]*x^2
plot(x, y6_calf_LT)

# Create functions for 4_mo and 6_mo treatments, to create regression lines based on these functions that appear continuous
func_4_calfLT <- function(x) {
  return(Coef4_calf_LT[1]+Coef4_calf_LT[2]*x+Coef4_calf_LT[3]*x^2)
}
func_6_calfLT <- function(x) {
  return(Coef6_calf_LT[1]+Coef6_calf_LT[2]*x+Coef6_calf_LT[3]*x^2)
}

# Extract the estimates for SE to use in the creation of error bars
#The following two functions produce a vector to extract the estimates
age4.vector <- function(x) {
  matrix(c(1,0,x,x^2,0,0,0), 7, 1)
}
age6.vector <- function(x) {
  age4.vector(x) + matrix(c(0,1,0,0,0,x,x^2), 7, 1)
}

### This function extracts the estimate and the SE
get.SE <- function(x, mod, age=4) {
  #mod is the fitted model in nlme
  if (age==4) y_calf_indicator_vector <- age4.vector(x)
  if (age==6) y_calf_indicator_vector <- age6.vector(x)
  est <- crossprod(y_calf_indicator_vector, summary(mod)$coef$fixed)
  var1 <- t(y_calf_indicator_vector)%*%vcov(mod)%*%y_calf_indicator_vector
  SE <- sqrt(var1)
  return(c(est, est-sqrt(var1), est+sqrt(var1), SE))
}

x.time = 1:11
age=4
M4 <- matrix(0, length(x.time), 4)
colnames(M4) = c("Est", "Est-SE", "Est+SE", "SE")
for(i.x in x.time) M4[i.x,1:4]=get.SE(i.x, LT.calf.hours, age)

age=6
M6 <- matrix(0, length(x.time), 4)
colnames(M6) = c("Est", "Est-SE", "Est+SE", "SE")
for(i.x in x.time) M6[i.x,1:4]=get.SE(i.x, LT.calf.hours, age)

SE <- c(M4[,4], M6[,4])

# Combine the 3 vectors into a new dataframe so we can plot it
df_calf_LT <- data.frame(x,y4_calf_LT,y6_calf_LT) %>%
  rename("y4_calf" = "y4_calf_LT", "y6_calf" = "y6_calf_LT") %>%
  pivot_longer(cols= c("y4_calf", "y6_calf"),
               names_to = "Treatment",
               values_to = "Lying_diff") %>%
  arrange(Treatment , x) %>%
  cbind(SE) %>%
  mutate(ymin=Lying_diff-SE, ymax = SE+Lying_diff)

# Plot the dataframe, overlaying the functions as regression lines that perfectly fit the discrete data points
calf_LT <- df_calf_LT %>%
  ggplot(aes(x = x, y = Lying_diff, color= Treatment)) +
  geom_point(shape = 16, size = 3) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1.5,  color = "light gray") +
  geom_function(fun = func_4_calfLT, color = "#ff835a", lwd = 1.3) +
  geom_function(fun = func_6_calfLT, color = "#085a73", lwd = 1.3) +
  theme_classic() +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text = element_text(color = "black"),
        text = element_text(size = 24, color = "black"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(),
        legend.position = "none") + #First number is position on x-axis, second is height 
  labs(#title = "Lying time",
    y = "Difference in daily lying time (h/d)",
    x = "Time (days)",
    tag = "A") +
  scale_color_manual(values = c("y4_calf" = "#ff835a", "y6_calf" = "#085a73"),
                     labels = c("4MO", "6MO")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(-8, -6, -4, -2, 0)) +
  scale_x_continuous(breaks = x, label = x) +
  geom_errorbar(aes(ymin=ymin,ymax=ymax), width =0.2, size=0.75, position = position_dodge(width = 0.3))

# View the final calf plot
calf_LT

### b) Cows (Fig 3B) -----------------------------------------------------------
# Run a model with lying time converted to h/d
LT.cow.hours <- lme((lying_diff/60) ~ trt + obs_day2 + I(obs_day2^2) + trt:obs_day2 + trt:I(obs_day2^2) + percent_near + parity2,
                    random = ~ 1 | id,
                    correlation = corAR1(form = ~ 1 | id),
                    data = IQ_full %>% dplyr::filter(cow_calf == 'cow'),
                    na.action=na.omit)

# Use results from summary to plot the polynomial
summary(LT.cow.hours)$coef$fixed

# First, we use the coefficients to plot data for the 4_mo group, since it is the reference in the model
# 4 mo 
x <- 1:11
Coef4_cow_LT <- summary(LT.cow.hours)$coef$fixed[c(1,3,4)] #Look at summary of fixed effects to see which belong only to 4_mo group
y4_cow_LT <- Coef4_cow_LT[1]+Coef4_cow_LT[2]*x+Coef4_cow_LT[3]*x^2
plot(x, y4_cow_LT)

# Then, we make the 6_mo plot by adding the effects of the 4_mo group to the fixed effects that only pertain to the 6_mo group
# (i.e., interactions of 6_mo and obs_day)
# 6 mo
x <- 1:11
Coef6_cow_LT <- summary(LT.cow.hours)$coef$fixed[c(2,7,8)] + Coef4_cow_LT
y6_cow_LT <- Coef6_cow_LT[1]+Coef6_cow_LT[2]*x+Coef6_cow_LT[3]*x^2
plot(x, y6_cow_LT)

# Extract the estimates for SE to use in the creation of error bars
#The following two functions produce a vector to extract the estimates
age4.vector <- function(x) {
  matrix(c(1,0,x,x^2,0,0,0,0), 8, 1)
}
age6.vector <- function(x) {
  age4.vector(x) + matrix(c(0,1,0,0,0,0,x,x^2), 8, 1)
}

# This is the main code calling the functions above
# Run your code first to fit the model poly.calf first
x.time = 1:11
age=4
M4 <- matrix(0, length(x.time), 4)
colnames(M4) = c("Est", "Est-SE", "Est+SE", "SE")
for(i.x in x.time) M4[i.x,1:4]=get.SE(i.x, LT.cow.hours, age)

age=6
M6 <- matrix(0, length(x.time), 4)
colnames(M6) = c("Est", "Est-SE", "Est+SE", "SE")
for(i.x in x.time) M6[i.x,1:4]=get.SE(i.x, LT.cow.hours, age)

SE <- c(M4[,4], M6[,4])

# Comine the 3 vectors into a new dataframe so we can plot it
df_cow_LT <- data.frame(x,y4_cow_LT,y6_cow_LT) %>%
  pivot_longer(cols= c("y4_cow_LT", "y6_cow_LT"),
               names_to = "Treatment",
               values_to = "Lying_diff") %>%
  arrange(Treatment , x) %>% 
  cbind(SE) %>% 
  mutate(ymin=Lying_diff-SE, ymax = SE+Lying_diff)

# Create functions for 4_mo and 6_mo treatments, to create regression lines based on these functions that appear continuous
func_4_cowLT <- function(x) {
  return(Coef4_cow_LT[1]+Coef4_cow_LT[2]*x+Coef4_cow_LT[3]*x^2)
}
func_6_cowLT <- function(x) {
  return(Coef6_cow_LT[1]+Coef6_cow_LT[2]*x+Coef6_cow_LT[3]*x^2)
}

# Plot the dataframe, overlaying the functions as regression lines that perfectly fit the discrete data points
cow_LT <- df_cow_LT %>%
  ggplot(aes(x = x, y = Lying_diff, color= Treatment)) +
  geom_hline(yintercept = 0, linetype = "longdash", size = 1.5,  color = "light gray") +
  geom_point(shape = 16, size = 4) +
  geom_function(fun = func_4_cowLT, color = "#ff835a", lwd = 1.3) +
  geom_function(fun = func_6_cowLT, color = "#085a73", lwd = 1.3) +
  theme_classic() +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        text = element_text(size = 24),
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(),
        legend.position = c(0.85, 0.2)) + #First number is position on x-axis, second is height 
  labs(#title = "Lying time",
    y = "Difference in daily lying time (h/d)",
    x = "Time (days)",
    tag = "B") +
  scale_color_manual(values = c("y4_cow_LT" = "#ff835a", "y6_cow_LT" = "#085a73"),
                     labels = c("4MO", "6MO")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1)) +
  scale_x_continuous(breaks = x, label = x) +
  geom_errorbar(aes(ymin=ymin,ymax=ymax), width =0.2, size=0.75, position = position_dodge(width = 0.3)) 

# View the final cow plot
cow_LT

## 2) Step count ----

# We ultimately need two datasets for both cows and calves: one with individual 
# cow- or calf-day datapoints, and another with only median values and IQR for each day

### a) Calves (Fig 4) ----------------------------------------------------------

# Isolate the calf data and remove the baseline values
# This dataset will contain the individual calf-day differences in step count
IQ_full_calf <- IQ_full %>%
  filter(cow_calf == 'calf') %>%
  filter(obs_day != "BL")

# Put back treatment-days with no data in IQ_full_calf (EXPANDED DATA SET 1)

# One value per observation day and treatment
trt_d <- expand.grid(unique(IQ_full_calf$obs_day), 
                     unique(IQ_full_calf$trt))
trt_d <- trt_d %>% rename(obs_day = Var1, trt = Var2)

# One value per calf and treatment
trt_id <- unique(dplyr::select(IQ_full_calf, id, trt))
d_full <- full_join(trt_d, trt_id)

# Join to get one value per calf-day and correct treatment
pl_IQ_full_calf <- left_join(d_full, IQ_full_calf, 
                                 by = c("id", "obs_day", "trt"))

# Create a second calf-only dataset that contains only the median values and IQR
IQ_median_calf <- IQ_full %>%
  filter(cow_calf == 'calf') %>%
  group_by(obs_day, trt) %>%
  summarise(median_steps_diff = median(steps_diff, na.rm = TRUE),
            IQR_steps_diff = IQR(steps_diff, na.rm = TRUE),
            Q1 = quantile(steps_diff, 0.25, na.rm = TRUE),
            Q3 = quantile(steps_diff, 0.75, na.rm = TRUE)) %>%
  filter(obs_day != "BL")

# Put back treatment-days with no data in IQ_median_calf (EXPANDED DATA SET 2)
# pl_ added to name to indicate only use only for plotting
pl_IQ_median_calf <- left_join(trt_d, IQ_median_calf, 
                               by = c("obs_day", "trt"))

# Plot the data
calf_raw_steps_median <- ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5,  
             color = "light gray") +
  geom_jitter(data=pl_IQ_full_calf, # EXPANDED DATA SET 1 - TO INDUCE INDIVIDUAL NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE POINTS OF THE OTHER GROUP TO KEEP CORRECT SIDE
              aes(x = obs_day, y=steps_diff, colour = trt), 
              size = 2, alpha = 1, pch = 1, # pch = force point format to what you want
              position = position_jitterdodge(0.6)) + 
  geom_errorbar(data=pl_IQ_median_calf, # EXPANDED DATA SET 2 - TO INDUCE SUMMARY NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE ERRORBARS OF THE OTHER GROUP TO KEEP CORRECT SIDE
                aes(x=obs_day, ymin=Q1, ymax=Q3, colour = trt), 
                width = 0.2, size = 0.8, 
                position = position_dodge(width = 0.8))+
  geom_hpline(data=pl_IQ_median_calf, # EXPANDED DATA SET 2 - TO INDUCE SUMMARY NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE MEDIAN BAR OF THE OTHER GROUP TO KEEP CORRECT SIDE
              aes(x=obs_day, y=median_steps_diff, colour = trt), 
              position = position_dodge(width = 0.8), width = 0.4, size = 2.5) +
  theme_classic() +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        text = element_text(size = 24, color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(),
        legend.position = c(0.85,0.9)) +
  labs(color = "Treatment",
       y = "Difference in daily steps (steps/d)",
       x = "Time (days)") +
  scale_color_manual(values = c("4_mo" = "#ff835a", "6_mo" = "#085a73"),
                     labels = c("4MO", "6MO")) +
  scale_y_continuous(breaks = c(-2000, 0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000))

# View the plot
calf_raw_steps_median

### a) Cows (Fig 5) ------------------------------------------------------------

# Isolate the cow data and remove the baseline values
# This dataset will contain the individual cow-day differences in step count
IQ_full_cow <- IQ_full %>%
  filter(cow_calf == 'cow') %>%
  filter(obs_day != "BL")

# Put back treatment-days with no data in IQ_full_cow (EXPANDED DATA SET 1)

# One value per observation day and treatment
trt_d <- expand.grid(unique(IQ_full_cow$obs_day), 
                     unique(IQ_full_cow$trt))
trt_d <- trt_d %>% rename(obs_day = Var1, trt = Var2)

# One value per cow and treatment
trt_id <- unique(dplyr::select(IQ_full_cow, id, trt))
d_full <- full_join(trt_d, trt_id)

# Join to get one value per cow-day and correct treatment
pl_IQ_full_cow <- left_join(d_full, IQ_full_cow, 
                                by = c("id", "obs_day", "trt"))

# Create a second cow dataset containing median and IQR values for each treatment and day
IQ_median_cow <- IQ_full %>%
  filter(cow_calf == 'cow') %>%
  group_by(obs_day, trt) %>%
  summarise(median_steps_diff = median(steps_diff, na.rm = TRUE),
            IQR_steps_diff = IQR(steps_diff, na.rm = TRUE),
            Q1 = quantile(steps_diff, 0.25, na.rm = TRUE),
            Q3 = quantile(steps_diff, 0.75, na.rm = TRUE)) %>%
  filter(obs_day != "BL")

# Put back treatment-days with no data in IQ_median_cow (EXPANDED DATA SET 2) 
# pl_ added to name to indicate only use only for plotting
pl_IQ_median_cow <- left_join(trt_d, IQ_median_cow, 
                              by = c("obs_day", "trt"))

# Plot the data
cow_raw_steps_median <- ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5,  
             color = "light gray") +
  geom_jitter(data=pl_IQ_full_cow, # EXPANDED DATA SET 1 - TO INDUCE INDIVIDUAL NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE POINTS OF THE OTHER GROUP TO KEEP CORRECT SIDE
              aes(x = obs_day, y=steps_diff, colour = trt), 
              size = 2, alpha = 1, pch = 1, # pch = force point format to what you want
              position = position_jitterdodge(0.6)) + 
  geom_errorbar(data=pl_IQ_median_cow, # EXPANDED DATA SET 2 - TO INDUCE SUMMARY NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE ERRORBARS OF THE OTHER GROUP TO KEEP CORRECT SIDE
                aes(x=obs_day, ymin=Q1, ymax=Q3, colour = trt), 
                width = 0.2, size = 0.8, 
                position = position_dodge(width = 0.8))+
  geom_hpline(data=pl_IQ_median_cow, # EXPANDED DATA SET 2 - TO INDUCE SUMMARY NA VALUES FOR DAYS WITH NO VALUES FOR ONE OF THE GROUPS, TO FORCE THE MEDIAN BAR OF THE OTHER GROUP TO KEEP CORRECT SIDE
              aes(x=obs_day, y=median_steps_diff, colour = trt), 
              position = position_dodge(width = 0.8), width = 0.4, size = 2.5) +
  theme_classic() +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        text = element_text(size = 24, color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(),
        legend.position = c(0.8,0.9)) +
  labs(color = "Treatment",
       y = "Difference in daily steps (steps/d)",
       x = "Time (days)") +
  scale_color_manual(values = c("4_mo" = "#ff835a", "6_mo" = "#085a73"),
                     labels = c("4MO", "6MO")) +
  scale_y_continuous(breaks = c(- 2000, 0, 2000, 4000, 6000, 8000))

# View the plot
cow_raw_steps_median