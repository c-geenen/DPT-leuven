library(tidyverse)
library(readxl)
library(lubridate)
library(investr)
library(DescTools)
library(MASS)

# Read data
cases <- read_excel("Input/Supplementary Data.xlsx", sheet="cases", na="NA", guess_max=10^5)
cases_age <- read_excel("Input/Supplementary Data.xlsx", sheet="cases_age", na="NA", guess_max=10^5)

# Initialise results tibble
cascade <- tibble(
  step = character(),
  successes = integer(),
  total = integer(),
  x90_est = numeric(),
  mean_est = numeric(),
  mean_log = numeric(),
  mean_low = numeric(),
  mean_high = numeric(),
  sdlog_est = numeric(),
  sdlog_se = numeric(),
  sdlog_low = numeric(),
  sdlog_high = numeric(),
  Asym_est = numeric(),
  Asym_se = numeric(),
  Asym_low = numeric(),
  Asym_high = numeric()
)

#############
# Describe included cases
############

paste0("n included cases: ", nrow(cases))

# Demographics
table(cases$Sex, useNA="always")
print(cases_age, width=Inf)
cases %>%
  summarise(
    prop_female = sum(Sex=="Vrouw", na.rm=T)/sum(!is.na(Sex), na.rm=T),
    missing_sex = sum(is.na(Sex))/n()
  )

################
# Delays
################

# PCR sample to result
delays <- cases %>%
  #filter(!is.na(PCR_sample_time),!is.na(PCR_result_time)) %>%
  reframe(
    ci = c("mid","low","high"),
    n=n(),
    sample_to_result_q10q50q90 = quantile(sample_to_result_time, c(0.5,0.1,0.9), na.rm=T),
    sample_to_result_median  = MedianCI(sample_to_result_time, conf.level=0.95, sides="two.sided", method="exact", na.rm=T),
    sample_to_result_mean  = MeanCI(sample_to_result_time, conf.level=0.95, sides="two.sided", method="classic", na.rm=T),
    sample_to_result_na = sum(is.na(sample_to_result_time)),
    result_to_interview_q10q50q90 = quantile(result_to_interview_time, c(0.5,0.1,0.9), na.rm=T),
    result_to_interview_median  = MedianCI(result_to_interview_time, conf.level=0.95, sides="two.sided", method="exact", na.rm=T),
    result_to_interview_mean  = MeanCI(result_to_interview_time, conf.level=0.95, sides="two.sided", method="classic", na.rm=T),
    result_to_interview_na = sum(is.na(result_to_interview_time))
  )
print(delays, width=Inf)

# for a fair comparison with DPT, where we assume a lognormal distribution:
# fit time from PCR result to interview to lognormal distribution
result_to_interview <- cases %>%
  filter(
    !is.na(result_to_interview_time),
    result_to_interview_time > 0
  )
(fit <- fitdistr(result_to_interview$result_to_interview_time, "lognormal"))
# check fit
ggplot(result_to_interview, aes(result_to_interview_time)) +
  geom_histogram(aes(y=after_stat(density)), binwidth=0.5, alpha=0.5) +
  geom_function(fun=dlnorm, args=list(fit$estimate["meanlog"],fit$estimate["sdlog"]), colour="red", linewidth=1, n=1000) +
  theme_bw() +
  labs(x="PCR result to case interview",y="Density")
# determine mean and its 95% CI of lognormal distribution on a log scale
meanlog_ci <- as.numeric(c(
  meanlog_low <- fit$estimate["meanlog"] - fit$sd["meanlog"]*1.96,
  fit$estimate["meanlog"],
  meanlog_high <- fit$estimate["meanlog"] + fit$sd["meanlog"]*1.96
))
# convert the mean to a linear scale
mean_ci <- exp( meanlog_ci + (fit$estimate["sdlog"]^2)/2)
mean_ci

##################
# Step 1: case is active user
#######################

## Show number of cases who are active users (by sex)
cases %>%
  group_by(Active_user) %>%
  summarise(
    n= n(),
    n_female=sum(Sex=="Vrouw", na.rm=T),
    n_male=sum(Sex=="Man", na.rm=T),
    prop_female=n_female/(n_male+n_female),
    n_na_sex=sum(is.na(Sex))
  )

## Show number of cases who are active users (with age)
print(cases_age, width=Inf)

# app use when data is missing
cases %>% filter(Missing_age) %>% summarise(app_use=sum(Active_user, na.rm=T)/sum(!is.na(Active_user)))
cases %>% filter(is.na(Sex)) %>% summarise(app_use=sum(Active_user, na.rm=T)/sum(!is.na(Active_user)))

# Add totals to summary
cascade <- cascade %>% add_row(
  step = "Case app use",
  successes = sum(cases$Active_user, na.rm=T), # n active users
  total = sum(!is.na(cases$Active_user)   # total included cases
))

########################
# Step 2: case uploads identifier
########################

# Filter only cases who are app users and were interviewed
cases_users_traced <- cases %>%
  ### exclude non-active-users
  filter(Active_user==TRUE) %>%
  ### Exclude cases not interviewed
  filter( Traced == TRUE ) %>%
  ### Exclude cases interviewed before test result was available
  filter( 
    sample_to_interview_days >= 0,
    result_to_interview_time >= 0 | is.na(result_to_interview_time)
  ) %>%
  ### Exclude cases with missing data on identifier upload
  filter( !is.na(Confirmed_linked) )
paste0( "cases surveyed on identifier upload status: ", nrow(cases_users_traced) )

## Show number of interviewed app users who uploaded their identifier (by sex)
cases_users_traced %>%
  group_by(Confirmed_linked) %>%
  summarise(
    #mean_age=mean(Age, na.rm=T),
    n= n(),
    n_female=sum(Sex=="Vrouw", na.rm=T),
    n_male=sum(Sex=="Man", na.rm=T),
    n_na_sex=sum(is.na(Sex)),
    #n_na_age=sum(is.na(Age))
  )

# Confidence interval for traced users who have uploaded the identifier (at any time after the PCR result)
binom.test(
  sum(cases_users_traced$Confirmed_linked, na.rm=T), # n active users
  sum(!is.na(cases_users_traced$Confirmed_linked)), # total cases with data
  100/100, alternative =  "two.sided", conf.level = 0.95)

##################################
# Step 2: identifier upload by case
# Plot step 2, by time from PCR test result to case interview
###################################

plot_data_step2 <- cases_users_traced %>%
  mutate(
    Confirmed_linked = as.integer(Confirmed_linked)
  )

# show NAs for timing, these cases are excluded from the graph
plot_data_step2 %>% group_by(is.na(result_to_interview_time)) %>%
  summarise(n_alert=sum(Confirmed_linked), n=n())
plot_data_step2_noNA <- plot_data_step2 %>%
  filter(!is.na(result_to_interview_time))

# Fit data to a cumulative distribution function multiplied with a maximum (the final proportion reached after time infinity)
step2_model <- nls(Confirmed_linked ~ Asym * plnorm(result_to_interview_time, meanlog, sdlog), plot_data_step2_noNA, start = list(meanlog = 1, sdlog = 1, Asym=0.5))
# Print fitted parameters
(coef <- coef(step2_model))
# determine confidence interval for parameters:
# Asym: the maximum proportion of cases that will eventually trigger AEN
# mean: the mean time to AEN triggering
param_raw <- summary(step2_model)$parameters %>%
  as_tibble() %>%
  mutate(parameter=rownames(summary(step2_model)$parameters)) %>%
  mutate(
    CI_low = Estimate - `Std. Error`*1.96,
    CI_high = Estimate + `Std. Error`*1.96
  ) %>%
  dplyr::select(parameter, Estimate,`Std. Error`,CI_low,CI_high)
param <- param_raw %>%
  add_row(
    parameter="mean",
    Estimate=exp(param_raw$Estimate[1] + (param_raw$Estimate[2]^2)/2),
    `Std. Error`=NA,
    CI_low=exp(param_raw$CI_low[1] + (param_raw$Estimate[2]^2)/2),
    CI_high=exp(param_raw$CI_high[1] + (param_raw$Estimate[2]^2)/2)
  )
param
# use the inverse formula to find the time x at which a certain proportion of identifier upload y is reached
time_to_proportion <- function(y) return( qlnorm(y/coef["Asym"], coef["meanlog"], coef["sdlog"]) )
# determine time to 10%, 50% and 90% of maximum, convert to hours
time_to_proportion(coef["Asym"]* c(0.1, 0.5, 0.90) ) * 24
# determine time to 50% and 90% of maximum
(percentiles <- time_to_proportion(coef["Asym"]* c(0.1, 0.5, 0.9) ) )

# plot model
new_data_step2 <- data.frame(result_to_interview_time=seq(0, 10, by = 0.01))
new_data_step2_pred <- as_tibble(predFit(step2_model, newdata = new_data_step2, interval = "confidence", level= 0.95)) %>% 
  mutate(result_to_interview_time = new_data_step2$result_to_interview_time) %>%
  rename(Confirmed_linked=fit)
plot1 <- ggplot(new_data_step2_pred, aes(result_to_interview_time, Confirmed_linked, ymin=lwr, ymax=upr)) +
  geom_line(aes(colour="Cumulative log-normal"), linewidth=1.5) +
  geom_ribbon(alpha=0.2, aes(fill="Cumulative log-normal")) +
  geom_point(
    data=plot_data_step2_noNA, aes(result_to_interview_time, Confirmed_linked), inherit.aes=F,
    alpha = 0.5, position = position_jitter(w=0, h=0.04)) +
  geom_smooth(
    data=plot_data_step2_noNA, aes(result_to_interview_time, Confirmed_linked, colour="LOESS", fill="LOESS"), inherit.aes=F,
    method = "loess", se = T, linewidth=1.5, alpha = 0.2, level=0.95) +
  geom_segment(yend=coef["Asym"]*0.9, xend=time_to_proportion(coef["Asym"]* 0.9), x=-1, y=coef["Asym"]*0.9, linewidth=0.5 ) +
  geom_segment(yend=coef["Asym"]*0.9, xend=time_to_proportion(coef["Asym"]* 0.9), x=time_to_proportion(coef["Asym"]* 0.9), y=-1, linewidth=0.5 ) +
  scale_colour_manual(name = "Fitted curve: ", values = c("blue","red"))  + 
  scale_fill_manual(name = "Fitted curve: ", values = c("blue","red")) +
  labs(y = "Probability of triggering AEN", x = "Days from PCR result to AEN trigger query") +
  coord_cartesian(ylim=c(0,1),xlim=c(0,7)) +
  theme_bw() +
  theme(
    legend.position = "top"
  )
plot1

ggsave("output/cascade_delays_1.pdf", width=4, height=4)

##########################################
# Step 2: Only include cases queried at least 13h after PCR result
############################################

# Confidence interval for traced users who have uploaded the identifier,
# given interview was at least 9.66 hours (the 90th percentile) after result was available
at_least_10h <- cases_users_traced %>%
  filter(
    # if exact PCR result and tracing time are available, require 13 hour interval
    (result_to_interview_time > time_to_proportion(coef["Asym"]*0.90) & !(is.na(result_to_interview_time))) | 
    # OR if exact times are not both available, make sure Tracing date is at least 2 days after PCR sample date
    (sample_to_interview_days >=  2 & is.na(result_to_interview_time))
  )

# add totals and estimates to summary
cascade <- cascade %>% add_row(
  step = "Case identifier upload",
  successes = sum(at_least_10h$Confirmed_linked, na.rm=T),
  total = sum(!is.na(at_least_10h$Confirmed_linked)),   # total included cases
  x90_est = percentiles[3],
  mean_est = param$Estimate[4],
  mean_log = param$Estimate[1],
  mean_low = param$CI_low[4],
  mean_high = param$CI_high[4],
  sdlog_est = param$Estimate[2],
  sdlog_se = param$`Std. Error`[2],
  sdlog_low = param$CI_low[2],
  sdlog_high = param$CI_high[2],
  Asym_est = param$Estimate[3],
  Asym_se = param$`Std. Error`[3],
  Asym_low = param$CI_low[3],
  Asym_high = param$CI_high[3]
)

################################################
### Read list of included case-contact-test trios ####
#################################################

# Load combinations of cases who uploaded their identifier, their manually traced contacts and
# the tests of these contacts at our test centre, within 14 days before or after case PCR sampling.
# (this dataframe also includes an observation for each case-contact pair where the contact did
# not book a test at our test centre.)
case_contact_test_trios <- read_excel("Input/Supplementary Data.xlsx", sheet="case_contact_test_trios", na="NA", guess_max=10^5) %>%
  
  ### Only include relevant preregistration forms
  ### Remove lines from the filter function to determine numbers of excluded pairs
  filter(
    # Contact must not already have a recent positive test
    !already_pos,
    # Contact must have a test within two weeks before or after the case PCR test
    !is.na(PCR_result_to_prereg),
    ### Pre-registration must be after PCR test result of case (assumed 21 PM same day if unknown)
    PCR_result_to_prereg >= 0,
    ## Exclude contacts already traced in the previous 7 days to a different case
    !already_exposed
  )

## Number of case_contact pairs
case_contact_test_trios %>%
  distinct(case_contact_pair_id) %>%
  nrow()

####################
# Delay from case PCR result to MCT contact preregistration
####################

delay_prereg <- case_contact_test_trios %>%
  reframe(
    ci = c("mid","low","high"),
    n=n(),
    PCR_result_to_prereg_q10q50q90 = quantile(PCR_result_to_prereg, c(0.1,0.5,0.9), na.rm=T),
    PCR_result_to_prereg_median  = MedianCI(PCR_result_to_prereg, conf.level=0.95, sides="two.sided", method="exact", na.rm=T),
    PCR_result_to_prereg_mean  = MeanCI(PCR_result_to_prereg, conf.level=0.95, sides="two.sided", method="classic", na.rm=T),
    PCR_result_to_prereg_na = sum(is.na(PCR_result_to_prereg)),
  )
print(delay_prereg, width=Inf)

##########################################
### Step 3: app use amongst contacts ###
##########################################

# Show how many contacts used the app, with confidence interval

# Summarise to list of contacts
contact_coronalert_use <- case_contact_test_trios %>%
  distinct(case_contact_pair_id, contact_has_coronalert)

table(contact_coronalert_use$contact_has_coronalert, useNA = "always")
cascade <- cascade %>% add_row(
  step = "Contact app use",
  successes = sum(contact_coronalert_use$contact_has_coronalert, na.rm=T),
  total = sum(!is.na(contact_coronalert_use$contact_has_coronalert))
)

##################################
# Plot step 4, by time from case PCR test result to contact test booking
###################################

# This plot shows one point per test, so allows multiple points per contact person.
plot_data_step4 <- case_contact_test_trios %>%
  filter(contact_has_coronalert) %>%
  mutate(
    contact_got_alert = as.integer(prereg_app_alert_received),
  ) %>%
  # maximum one observation per person per day
  group_by(round(PCR_result_to_prereg,0),case_contact_pair_id) %>%
  arrange(desc(PCR_result_to_prereg)) %>%
  slice(1) %>%
  ungroup

# Fit data to a cumulative distribution function multiplied with a maximum (the final proportion reached after time infinity)
step4_model <- nls(contact_got_alert ~ Asym * plnorm(PCR_result_to_prereg, meanlog, sdlog), plot_data_step4, start = list(meanlog = 1, sdlog = 1, Asym=0.5))
# Print fitted parameters
summary(step4_model)
(coef <- coef(step4_model))
# determine confidence interval for parameters
param_raw <- summary(step4_model)$parameters %>%
  as_tibble() %>%
  mutate(parameter=rownames(summary(step4_model)$parameters)) %>%
  mutate(
    CI_low = Estimate - `Std. Error`*1.96,
    CI_high = Estimate + `Std. Error`*1.96
  ) %>%
  dplyr::select(parameter, Estimate,`Std. Error`,CI_low,CI_high)
param <- param_raw %>%
  add_row(
    parameter="mean",
    Estimate=exp(param_raw$Estimate[1] + (param_raw$Estimate[2]^2)/2),
    `Std. Error`=NA,
    CI_low=exp(param_raw$CI_low[1] + (param_raw$Estimate[2]^2)/2),
    CI_high=exp(param_raw$CI_high[1] + (param_raw$Estimate[2]^2)/2)
  )
param
# use the inverse formula to find the time x at which a certain proportion of identifier upload y is reached
time_to_proportion <- function(y) return( qlnorm(y/coef["Asym"], coef["meanlog"], coef["sdlog"]) )
# determine time to 50% and 95% of maximum, with maximum given by a, convert to hours
(percentiles <- time_to_proportion(coef["Asym"]* c(0.1, 0.5, 0.9) ) )
percentiles * 24

# plot model
new_data <- data.frame(PCR_result_to_prereg=seq(0, 10, by = 0.01))
new_data_pred <- as_tibble(predFit(step4_model, newdata = new_data, interval = "confidence", level= 0.95)) %>% 
  mutate(PCR_result_to_prereg = new_data$PCR_result_to_prereg) %>%
  rename(contact_got_alert=fit)
plot2 <- ggplot(new_data_pred, aes(PCR_result_to_prereg, contact_got_alert, ymin=lwr, ymax=upr)) +
  geom_line(aes(colour="Cumulative log-normal"), linewidth=1.5) +
  geom_ribbon(alpha=0.2, aes(fill="Cumulative log-normal")) +
  geom_point(
    data=plot_data_step4, aes(PCR_result_to_prereg, contact_got_alert), inherit.aes=F,
    alpha = 0.5, position = position_jitter(w=0, h=0.04)) +
  geom_smooth(
    data=plot_data_step4, aes(PCR_result_to_prereg, contact_got_alert, colour="LOESS", fill="LOESS"), inherit.aes=F,
    method = "loess", se = T, linewidth=1.5, alpha = 0.2, level=0.95) +
  geom_segment(yend=coef["Asym"]*0.9, xend=time_to_proportion(coef["Asym"]* 0.9), x=-1, y=coef["Asym"]*0.9, linewidth=0.5 ) +
  geom_segment(yend=coef["Asym"]*0.9, xend=time_to_proportion(coef["Asym"]* 0.9), x=time_to_proportion(coef["Asym"]* 0.9), y=-1, linewidth=0.5 ) +
  scale_colour_manual(name = "Fitted curve: ", values = c("blue","red"))  + 
  scale_fill_manual(name = "Fitted curve: ", values = c("blue","red")) +
  labs(y = "Probability of AEN receipt", x = "Days from PCR result to AEN receipt query") +
  coord_cartesian(ylim=c(0,1),xlim=c(0,7)) +
  theme_bw() +
  theme(
    legend.position = "top"
  )
plot2

ggsave("output/cascade_delays_2.pdf", width=4, height=4)

################################################
### Step 4: receipt of notification ###
################################################

# Summarise the test observations of each contact into one observation per case-contact pair
app_using_contacts <- case_contact_test_trios %>%
  
  filter(contact_has_coronalert) %>%
  
  # Based on the plot in the previous section, only include preregistration forms
  # at least 53 hours after the case PCR test result.
  filter(PCR_result_to_prereg>=percentiles[3]) %>%

  ### Each combination of source_ID and PCR_date is a unique infection
  ### (some index cases had multiple infections, but if it was less than
  ### 60 days before, it has already been removed from the input file.)
  group_by(case_contact_pair_id) %>%
  
  ### Sort preregistration forms by date and time
  arrange(PCR_result_to_prereg) %>%
  
  ### For each unique contact of a unique infection of a case, summarise to one row.
  summarise(
    contact_got_alert=any(prereg_app_alert_received),
  ) %>%
  ungroup

# Show how many reveived AEN
table(app_using_contacts$contact_got_alert, useNA = "always")

cascade <- cascade %>% add_row(
  step = "Contact AEN receipt",
  successes = sum(app_using_contacts$contact_got_alert, na.rm=T),
  total = sum(!is.na(app_using_contacts$contact_got_alert)),
  x90_est = percentiles[3],
  mean_est = param$Estimate[4],
  mean_log = param$Estimate[1],
  mean_low = param$CI_low[4],
  mean_high = param$CI_high[4],
  sdlog_est = param$Estimate[2],
  sdlog_se = param$`Std. Error`[2],
  sdlog_low = param$CI_low[2],
  sdlog_high = param$CI_high[2],
  Asym_est = param$Estimate[3],
  Asym_se = param$`Std. Error`[3],
  Asym_low = param$CI_low[3],
  Asym_high = param$CI_high[3]
)

###########################
# Calculate probabilities and save results to file
#########################

# Function for determining exact binomial confidence interval
binom_CI <- function(x, n, p, l) {
  if (n>0) {
    CI <- binom.test(x, n, p/100, alternative =  "two.sided", conf.level = 0.95)$conf.int[l]
    return(CI)
  } else {
    return(c(0,1)[l])
  }
}

cascade_summary <- cascade %>%
  rowwise %>%
  mutate(
    p_est = successes/(total),
    p_CI_low = binom_CI(successes, total, 0, 1),
    p_CI_high = binom_CI(successes, total, 0, 2)
  ) %>%
  ungroup

write_csv(cascade_summary,"output/cascade_summary.csv")

