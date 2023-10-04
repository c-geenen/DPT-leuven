# devtools::install_github("HopkinsIDD/tti")
library(tti)
library(tidyverse)
library(patchwork)
library(ggtext)

cascade <- read_csv("output/cascade_summary.csv")

symptoms_to_test = 1.9 # mean time from symptoms to sampling for symptomatic screening (calculated in SAR_analysis.R)
test_to_result = 0.44 # mean time from sampling to positive result (calculated in cascade_analysis.R)
result_to_AEN = cascade$mean_est[4] # mean time from case PCR to digital notification (calculated cascade_analysis.R)
result_to_MCT = 2.264929 # mean time from case PCR result to manual contact notification, as determined using a fitted lognormal distribution (cascade_analysis.R)

DPT_succes_rate = 0.043
MCT_succes_rate = DPT_succes_rate / 60 * 616
combo_succes_rate = DPT_succes_rate / 60 * 648

################
# Determine mean time from case PCR to notifying contacts through either
# DPT or MCT in a combined strategy.
####################

# Probability of being notified by time from PCR test result
# (to get confidence interval, plug in upper and lower bounds in result_to_AEN and result_to_MCT above)
prob_MCT_notified <- function(t) MCT_succes_rate * plnorm(t, log(result_to_MCT)-(1.12706035^2)/2, sdlog=1.12706035 ) # values from cascade_analysius.R
prob_DPT_notified <- function(t) DPT_succes_rate * plnorm(t, log(result_to_AEN)-(cascade$sdlog_est[4]^2)/2, cascade$sdlog_est[4])
prob_combo_notified <- function(t) prob_MCT_notified(t) + prob_DPT_notified(t) - prob_MCT_notified(t)*prob_DPT_notified(t)
# visualise proportion notified
ggplot() +
  geom_function(fun=prob_MCT_notified, aes(colour="MCT")) +
  geom_function(fun=prob_DPT_notified, aes(colour="DPT")) +
  geom_function(fun=prob_combo_notified, aes(colour="Combo")) +
  xlim(0,7)
# determine probability density and visualise
combo_notify_delay <- tibble( t = seq(0,1000, by=0.0001) ) %>%
  mutate(
    prob_cum = prob_combo_notified(t),
    prob = replace_na(prob_cum - lag(prob_cum),0)
  )
ggplot(combo_notify_delay, aes(t,prob)) + geom_line() + xlim(c(0,7))
# determine mean
(result_to_AEN_or_MCT  = as.numeric( summarise(combo_notify_delay, mean_delay = sum(t*prob) / sum(prob)) ))
  
####################
# Run the model
###################
  
get_r_for_each_strategy <- function(
    R,
    rho
  ) {
  
  strategies <- c("base","isolation","MCT","DPT","combo")
  results <- tibble()
  
  for (strategy in strategies) {
    print(strategy)
    
    # Probability of being traced given (community or household) exposure
    omega_c <- case_when(
      strategy == "base" ~ 0,
      strategy == "isolation" ~ 0,
      strategy == "MCT" ~ MCT_succes_rate, # include CI
      strategy == "DPT" ~ DPT_succes_rate, # include CI
      strategy == "combo" ~ combo_succes_rate, # include CI
      TRUE ~ NA
    )
    omega_h <- omega_c
    omega_q <- omega_c
    
    # Time delay from index cases symptom onset to quarantine of community contacts.
    t_qcs <- case_when(
      strategy == "base" ~ 100,
      strategy == "isolation" ~ 100,
      strategy == "MCT" ~ symptoms_to_test + test_to_result + result_to_MCT,
      strategy == "DPT" ~ symptoms_to_test + test_to_result + result_to_AEN,
      strategy == "combo" ~ symptoms_to_test + test_to_result + result_to_AEN_or_MCT,
      TRUE ~ NA
    )
    t_qhs <- t_qcs
    t_q <- t_qcs
    t_qca <- 100 #t_qcs + additional_delay_asymp
    t_qha <- 100 #t_qhs + additional_delay_asymp
    
    # Time delay from symptom onset to isolation in  person detected through screening
    t_ds <- case_when(
      strategy == "base" ~ 100,
      strategy == "isolation" ~ symptoms_to_test + test_to_result,
      strategy == "MCT" ~ symptoms_to_test + test_to_result,
      strategy == "DPT" ~ symptoms_to_test + test_to_result,
      strategy == "combo" ~ symptoms_to_test + test_to_result,
      TRUE ~ NA
    )
    t_da = 100 #t_ds + additional_delay_asymp
    
    df <- get_r_effective_df(
      alpha = 0.2, # Probability of asymptomatic infection
      R = R,
      kappa = 0.35, # Relative transmissibility of asymptomatic individual
      eta = 0, # All contacts are considered community contacts
      nu = 1, # No differentiation between household/community contacts
      t_ds = t_ds, # Time delay from symptom onset to isolation in detected symptomatic person
      t_da = t_da, # Time delay from symptom onset to isolation in detected asymptomatic person
      t_qcs = t_qcs,
      t_qca = t_qca,
      t_qhs = t_qhs,
      t_qha = t_qha,
      t_q = t_q,
      omega_c = omega_c,
      omega_h = omega_h,
      omega_q = omega_q,
      quarantine_days = 7,
      isolation_days = 10,
      rho_s = rho,
      rho_a = 0,
      t_incubation = 4.2,
      offset = -2.31,
      shape = 1.65,
      rate = 0.5,
      stoch = F,
      theta = 0.1,
      n_inf = 20,
      n_iter = 100000
    ) %>%
      dplyr::select(r_effective,prop_identified) %>%
      mutate(
        strategy=strategy,
        R=R,
        rho=rho
      )
    
    results <- bind_rows(results,df)
  }
  
  return(results)
}



########################
# Vary rho
######################

iterations <- tibble()

set.seed(1)

for (rho in seq(0.1,0.9,by=0.1)) {
  print(rho)
  
  for (R in c(1.5)) {
    results <- get_r_for_each_strategy(R=R,rho=rho)
    
    iterations <- bind_rows(iterations,results)
  }
}

# if using computation-intensive stochastic model, use these lines to save and load the results
#write_csv(iterations, "output/iterations_r_effective_200k.csv")
#iterations <- read_csv("output/iterations_r_effective_200k.csv")

summarise_iterations <- iterations %>%
  group_by(strategy, rho, R) %>%
  summarise(
    r_effective = mean(r_effective),
    prop_identified = mean(prop_identified)
  ) %>%
  ungroup %>%
  mutate(
    strategy = recode_factor(
      strategy,
      "base" = "No intervention",
      "isolation" = "Case isolation only",
      "DPT" = "Case isolation and DPT",
      "MCT" = "Case isolation and MCT",
      "combo" = "Case isolation, DPT and MCT"
    )
  )

plot1 <- ggplot(summarise_iterations,aes(rho,r_effective,colour=strategy,linetype=strategy)) +
  #facet_wrap(~paste0("R = ",R)) +
  #geom_point(alpha=0.7) +
  geom_line(linewidth=1) +
  coord_cartesian(ylim=c(1,1.5),xlim=c(0,1)) +
  theme_bw() +
  geom_hline(yintercept=1,colour="grey") +
  labs(
    colour = "Tracing strategy",
    linetype = "Tracing strategy",
    y = "R<sub>eff</sub>",
    x = "Proportion isolated through symptomatic screening"
  ) +
  theme(
    axis.title.y = element_markdown(),
    axis.title.x = element_blank()
  )
plot1

relative_to_isolation <- summarise_iterations %>%
  group_by(rho, R) %>%
  arrange(strategy) %>%
  mutate(
    r_diff_case_isolation = (r_effective[2]-R)/R,
    r_diff = (r_effective - R) / R,
    r_diff_multiplier = r_diff / r_diff_case_isolation
  ) %>%
  ungroup

plot2 <- ggplot(relative_to_isolation,aes(rho,r_diff_multiplier,colour=strategy,linetype=strategy)) +
  #facet_wrap(~paste0("R = ",R)) +
  #geom_point(alpha=0.7) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim=c(0,2), xlim=c(0,1)) +
  theme_bw() +
  geom_hline(yintercept=0,colour="grey") +
  labs(
    colour = "Tracing strategy",
    linetype = "Tracing strategy",
    y = "Reduction of R<sub>eff</sub> relative to isolation only",
    x = "Proportion identified through symptomatic screening"
  ) +
  theme(
    axis.title.y = element_markdown()
  )
plot2

(plot <- plot1 / plot2 + plot_layout(guides = 'collect'))

ggsave("output/r_effective_varying_rho.png",plot=plot,width=6,height=6)

mean_relative_to_isolation <- relative_to_isolation %>%
  group_by(strategy) %>%
  summarise(
    r_diff_multiplier_mean = mean(r_diff_multiplier),
    r_diff_multiplier_q10 = quantile(r_diff_multiplier,0.1),
    r_diff_multiplier_q90 = quantile(r_diff_multiplier,0.9)
  )
mean_relative_to_isolation
# effect of DPT only compared to MCT
(mean_relative_to_isolation$r_diff_multiplier_mean[3]-1) / (mean_relative_to_isolation$r_diff_multiplier_mean[4] - 1)
# additional effect of DPT in a combined strategy
(mean_relative_to_isolation$r_diff_multiplier_mean[5]-1) / (mean_relative_to_isolation$r_diff_multiplier_mean[4] - 1)

#################
# Vary R (not used)
################@

iterations <- tibble()

for (R in seq(0.1,3,by=0.1)) {
  print(R)
  
  for (rho in c(0.1,0.3,0.5)) {
    results <- get_r_for_each_strategy(R=R,rho=rho)
    
    iterations <- bind_rows(iterations,results)
  }
}

summarise_iterations <- iterations %>%
  group_by(strategy, R, rho) %>%
  summarise(
    r_effective = mean(r_effective),
    prop_identified = mean(prop_identified)
  )

ggplot(summarise_iterations,aes(R,r_effective,colour=strategy)) +
  facet_wrap(~paste0("rho = ",rho), ncol=2) +
  geom_point(alpha=0.7) +
  geom_line() +
  coord_equal() +
  theme_bw() +
  theme(legend.position=c(0.8,0.2))



########################
# Vary rho and R (not used)
######################

iterations <- tibble()

for (rho in seq(0,1,by=0.05)) {
  print(rho)
  
  for (R in seq(1,2.5,by=0.1)) {
    results <- get_r_for_each_strategy(R=R,rho=rho)
    
    iterations <- bind_rows(iterations,results)
  }
}

summarise_iterations <- iterations %>%
  group_by(strategy, rho, R) %>%
  summarise(
    r_effective = mean(r_effective),
    prop_identified = mean(prop_identified)
  )

ggplot(summarise_iterations,aes(rho,R,fill=r_effective)) +
  geom_raster() +
  scale_fill_gradient2(low = "green", high = "red",mid="white", midpoint=1) +
  scale_y_continuous() +
  facet_wrap(~strategy) +
  theme_bw()
