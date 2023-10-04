library(tidyverse)
library(epitools)
library(lubridate)
library(patchwork)
library(DescTools)

# Function for determining exact binomial confidence interval
binom_CI <- function(x, n, p, l) {
  if (n>0) {
    CI <- binom.test(x, n, p/100, alternative =  "two.sided", conf.level = 0.95)$conf.int[l]
    return(CI)
  } else {
    return(c(0,1)[l])
  }
}

# Function for determining confidence interval for ratio of two proportions
risk_ratio_CI <- function(pos_test, neg_test, pos_control, neg_control, l) {
  if ((neg_test+pos_test>0) & (neg_control+pos_control>0)) {
    CI <- riskratio(c(neg_control, pos_control, neg_test, pos_test), conf.level=0.95, method="small")
    output <- c( CI$measure[2], CI$measure[4], CI$measure[6], CI$p.value[6])
    return(output[l])
  } else {
    return(c(0,0,0,0)[l])
  }
}


##########################
### Load input files and check basic stats #####
#########################

observations <- read_csv("Input/tests.csv") %>%
  mutate(
    app_status = if_else(
      has_coronalert,
      if_else(alert_received,"alert","no alert"),
      "non-user"
    )
  ) %>%
  mutate(
    infected=if_else(
      outcome=="pos",
      1,
      if_else(
        outcome=="neg",
        0,
        NA
      )
    ),
    indication = case_match(indication, "close_contact" ~ "MCT", "coronalert" ~ "AEN", "symptoms" ~ "Symptoms",  .default = NA)
  ) %>%
  
  # Only include study period (change to FALSE to generate data for extended study periiod)
  filter(study_period==T) %>%   # only the period used in the analysis of the notification cascade
  # Leave out missing test indications
  filter(
    !is.na(has_coronalert),
    !is.na(indication)
  ) %>%
  # Exclude bookings with a positive test in the last 60 days.
  filter( Already_pos == FALSE ) %>%
  # Exclude bookings with a previous test in the last 14 days.
  filter( Already_tested == FALSE ) %>%
  # Exclude bookings with a previous booking in the last 14 days.
  filter( Already_booked == FALSE ) %>%
  # Exclude conflicting info
  filter( !(indication=="AEN" & app_status!="alert") ) %>%
  ungroup

# app use over time
app_use <- observations %>% group_by(week) %>%
  summarise(
    n = n(),
    n_users = sum(!app_status=="non-user"),
    user_rate = n_users/n,
  )
app_use
# app use overall
binom.test(sum(app_use$n_users), sum(app_use$n), 0/100, alternative =  "two.sided", conf.level = 0.95)

# app use and alert based on indication
observations %>% group_by(indication) %>%
  summarise(
    n = n(),
    n_users = sum(!app_status=="non-user"),
    user_rate = n_users/n,
    n_alert = sum(app_status=="alert"),
    alert_rate = n_alert/n_users,
    alert_rate_total = n_alert/n,
    pos_rate = sum(outcome=="pos")/sum(outcome!="l2fu")
  )

# alert based on indication and outcome
observations %>% group_by(indication,outcome) %>%
  summarise(
    n = n(),
    n_users = sum(!app_status=="non-user"),
    user_rate = n_users/n,
    n_alert = sum(app_status=="alert"),
    alert_rate_users = n_alert/n_users,
    alert_rate_total = n_alert/n
  )

########################
# Delay symptoms to sample
#######################

delay <- observations %>%
  reframe(
    ci = c("mid","low","high"),
    n=n(),
    symptoms_to_sample_q10q50q90 = quantile(symptoms_to_sample, c(0.1,0.5,0.9), na.rm=T),
    symptoms_to_sample_median  = MedianCI(symptoms_to_sample, conf.level=0.95, sides="two.sided", method="exact", na.rm=T),
    symptoms_to_sample_mean  = MeanCI(symptoms_to_sample, conf.level=0.95, sides="two.sided", method="classic", na.rm=T),
    symptoms_to_sample_na = sum(is.na(symptoms_to_sample)),
  )
delay

###########################
# Generate totals for main results table
##########################

summary <- observations %>%
  # mutate(app_status="any") %>%    # use this line to get positivity rates for each group
  group_by(indication, app_status) %>%
  summarise(
    n_neg = sum(outcome=="neg"),
    n_pos = sum(outcome=="pos"),
    n_l2fu = sum(outcome=="l2fu")
  ) %>%
  ungroup
summary

analysis <- summary %>%
  bind_rows(
    summary %>%
      filter(
        app_status == "alert",
        indication != "Symptoms"
      ) %>%
      summarise(
        indication = "AEN with or without MCT",
        n_neg = sum(n_neg),
        n_pos = sum(n_pos),
        n_l2fu = sum(n_l2fu)
      ),
    summary %>%
      filter(
        indication == "MCT"
      ) %>%
      summarise(
        indication = "MCT with or without AEN",
        n_neg = sum(n_neg),
        n_pos = sum(n_pos),
        n_l2fu = sum(n_l2fu)
      )
  ) %>%
  
  mutate(n=n_pos+n_neg) %>%				# n is sum of negatives and positives
  mutate(pos_rate = n_pos/n) %>%				# add column for positivity rate
  
  # use non-users with symptoms as reference for symptomatic group. use non-users manually traced as reference for all others
  mutate(ref = case_when(
    indication=="Symptoms" ~ "symp",
    indication=="AEN" | indication=="MCT" ~ "non-user",
    .default = "overlapping"
  )) %>%
  group_by(ref) %>%
  mutate(
    n_neg_ref = max(n_neg), # reference are non-app users in the same group. there are always more non-app users.
    n_pos_ref = max(n_pos)
  ) %>%
  
  # calculate SAR confidence intervals, relative risks and their confidence intervals, and p-values.
  rowwise() %>%
  mutate(CI_low = binom_CI(n_pos, n, 0, 1),		# add columns for upper and lower bounds of confidence interval
         CI_high = binom_CI(n_pos, n, 0, 2)) %>%
   mutate(
        RR = risk_ratio_CI(n_pos, n_neg, n_pos_ref, n_neg_ref, 1),
        CI_RR_low = risk_ratio_CI(n_pos, n_neg, n_pos_ref, n_neg_ref, 2),		# add columns for upper and lower bounds of confidence interval
        CI_RR_high = risk_ratio_CI(n_pos, n_neg, n_pos_ref, n_neg_ref, 3),
        p = risk_ratio_CI(n_pos, n_neg, n_pos_ref, n_neg_ref, 4)
  ) %>%
  ungroup %>%
  
  mutate(app_status=fct_relevel(as.factor(app_status),c("non-user","alert","no alert"))) %>%
  mutate(indication=fct_relevel(as.factor(indication),c("MCT","AEN","Symptoms"))) %>%
  arrange(indication, app_status) %>%
  
  dplyr::select(-n_neg_ref, -n_pos_ref)

print(analysis, n=50)

write_csv(analysis, "output/positivity_rates.csv")

###########################
# Changes of positive predictive value (risk ratio) over time
##########################

indication_clert <- observations %>%
  drop_na(indication)
indication_clert

plot_pos_rate <- ggplot(indication_clert,aes(week,infected, colour=indication)) +
  #geom_jitter(width=1,height=0.02,size=0.1, alpha=0.3) +
  geom_smooth(aes(fill=indication),span=0.3,alpha=0.2) +
  labs(y="Positivity rate",x="Date") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.5)) +
  theme(
    legend.position = "none"
  )

plot_n_test_bookings <- ggplot(indication_clert,aes(week, fill=indication)) +
  geom_bar() +
  labs(y="Number of tests",x="Date", fill="Test indication") +
  theme_bw()

plot_pos_rate | plot_n_test_bookings

###########################
# Same results table, but showing AEN rate by test indication
##########################

summary <- observations %>%
  # mutate(app_status="any") %>%    # use this line to get positivity rates for each group
  # mutate(outcome = if_else(indication==NA,"any",outcome)) %>%
  group_by(indication, outcome) %>%
  summarise(
    n = n(),
    n_no_aen = sum(app_status=="no alert"),
    n_alert = sum(app_status=="alert"),
    n_non_user = sum(app_status=="non-user")
  )
summary

analysis <- summary %>%
  
  mutate(n=n_alert+n_no_aen) %>%				# n is sum of negatives and positives
  mutate(aen_rate = n_alert/n) %>%				# add column for positivity rate
  
  mutate(
    # use people with symptoms but testing negative as reference
    n_no_aen_ref = as.numeric(summary[summary$indication == "Symptoms" & summary$outcome == "neg","n_no_aen"]),
    n_aen_ref = as.numeric(summary[summary$indication == "Symptoms" & summary$outcome == "neg","n_alert"]),
  ) %>%
  
  # calculate SAR confidence intervals, relative risks and their confidence intervals, and p-values.
  rowwise() %>%
  mutate(CI_low = binom_CI(n_alert, n, 0, 1),		# add columns for upper and lower bounds of confidence interval
         CI_high = binom_CI(n_alert, n, 0, 2)) %>%
  mutate(
    RR = risk_ratio_CI(n_alert, n_no_aen, n_aen_ref, n_no_aen_ref, 1),
    CI_RR_low = risk_ratio_CI(n_alert, n_no_aen, n_aen_ref, n_no_aen_ref, 2),		# add columns for upper and lower bounds of confidence interval
    CI_RR_high = risk_ratio_CI(n_alert, n_no_aen, n_aen_ref, n_no_aen_ref, 3),
    p = risk_ratio_CI(n_alert, n_no_aen, n_aen_ref, n_no_aen_ref, 4)
  ) %>%
  ungroup %>%
  
  mutate(outcome=fct_relevel(as.factor(outcome),c("neg","pos","l2fu"))) %>%
  mutate(indication=fct_relevel(as.factor(indication),c("Symptoms","MCT","AEN"))) %>%
  arrange(indication, outcome) %>%
  
  select(-n_no_aen_ref, -n_aen_ref)

print(analysis, n=50)

# proportion of AEN receipt amongst persons traced with MCT (irrespective of app use)
binom.test(
  sum(summary[summary$indication=="MCT",]$n_alert), 
  sum(summary[summary$indication=="MCT",]$n),
  0/100, alternative =  "two.sided", conf.level = 0.95
)

