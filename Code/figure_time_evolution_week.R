library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)

loess_span = 0.2
colours = c("darkorange3","#00BA38", "blue2")

# Function for determining exact binomial confidence interval
binom_CI <- function(x, n, p, l) {
  if (n>0) {
    CI <- binom.test(x, n, p/100, alternative =  "two.sided", conf.level = 0.95)$conf.int[l]
    return(CI)
  } else {
    return(c(0,1)[l])
  }
}

###########################
# read data file
#######################

# Read data file, exclude conflicting info
observations <- read_excel("Input/Supplementary Data.xlsx", sheet="tests", na="NA", guess_max=10^5) %>%
  # Exclusions
  filter(
    !is.na(has_coronalert),
    !is.na(indication)
  ) %>%
  filter( Already_pos == FALSE ) %>%
  filter( Already_tested == FALSE ) %>%
  filter( Already_booked == FALSE ) %>%
  filter( !(indication=="coronalert" & !alert_received) ) %>%
  
  # add timepoint column
  mutate(
    indication = case_match(
      indication,
      "symptoms" ~ "Symptomatic screening",
      "close_contact" ~ "Manually traced",
      "coronalert" ~ "Digitally traced",
    ),
    indication = fct_relevel(indication, c("Symptomatic screening","Manually traced","Digitally traced")),
    timepoint=as_date(week) # can be changed to "week"
  ) %>%
  filter(!is.na(indication))

plot_test_indication_absolute <- ggplot(filter(observations,outcome=="pos"), aes(x=timepoint,fill=indication)) +
  geom_bar(width=7, alpha=0.7) +
  labs(y="Cases confirmed at test centre", x = "Date", fill="Test indication") +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  theme(
    legend.position = c(0.2,0.7)
  ) +
  scale_fill_manual(values=colours)
plot_test_indication_absolute

plot_test_indication_relative <- ggplot(filter(observations,outcome=="pos"), aes(x=timepoint,fill=indication)) +
  geom_bar(position="fill", width=7, alpha=0.7) +
  labs(y="Proportion of detected cases", x = "Date") +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  theme(
    legend.position = "none"
  ) +
  scale_fill_manual(values=colours)
plot_test_indication_relative

####################
# Positivity rates
##################

pos_rates <- observations %>%
  group_by(timepoint, indication) %>%
  summarise(
    n_pos = sum(outcome=="pos"),
    n_not_l2fu = sum(outcome!="l2fu"),
    pos_rate = n_pos/n_not_l2fu,
  ) %>%
  rowwise %>%
  mutate(
    # add confidence intervals
    pos_rate_lower = binom_CI(n_pos, n_not_l2fu, 0, 1),
    pos_rate_upper = binom_CI(n_pos, n_not_l2fu, 0, 2)
  )

plot_pos_rate_both <- ggplot(filter(pos_rates,indication!="Symptomatic screening"),aes(x=timepoint,y=pos_rate,ymin=pos_rate_lower,ymax=pos_rate_upper, fill=indication,colour=indication)) +
  #geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  #facet_wrap(~indication) +
  theme_bw() +
  theme(legend.position = c(0.4,0.7)) +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  geom_vline(xintercept=ymd(c("2021-04-26")), colour=colours[3]) +
  labs(y = "Infection risk", x = "Date", colour="Test indication", fill="Test indication") +
  coord_cartesian(ylim=c(0,0.65)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95) +
  scale_colour_manual(values=c(colours[2],colours[3])) +
  scale_fill_manual(values=c(colours[2],colours[3]))
plot_pos_rate_both

####################
# App use
####################

app_use <- observations %>%
  filter(!indication=="Digitally traced") %>%
  group_by(timepoint) %>%
  summarise(
    # app use amongst all included persons
    n_users = sum(has_coronalert),
    n_total = n(),
    app_use_est = n_users/n_total,
    # proportion receiving AEN amongst HRC
    n_hrc_aen = sum(indication=="Manually traced" & alert_received),
    n_hrc = sum(indication=="Manually traced"),
    aen_rate_hrc_est = n_hrc_aen/n_hrc
  ) %>%
  rowwise %>%
  mutate(
    # add confidence intervals
    app_use_lower = binom_CI(n_users, n_total, 0, 1),
    app_use_upper = binom_CI(n_users, n_total, 0, 2),
    aen_rate_hrc_lower = binom_CI(n_hrc_aen, n_hrc, 0, 1),
    aen_rate_hrc_upper = binom_CI(n_hrc_aen, n_hrc, 0, 2),
  )

plot_app_use <- ggplot(app_use,aes(x=timepoint,y=app_use_est,ymin=app_use_lower,ymax=app_use_upper)) +
  geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "App uptake", x = "Date") +
  coord_cartesian(ylim=c(0,0.7)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[3], colour=colours[3])
plot_app_use 

plot_aen_rate <- ggplot(app_use,aes(x=timepoint,y=aen_rate_hrc_est,ymin=aen_rate_hrc_lower,ymax=aen_rate_hrc_upper)) +
  geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  geom_vline(xintercept=ymd(c("2021-04-26")), colour=colours[3]) +
  labs(y = "Proportion of MCT contacts with AEN", x = "Date") +
  coord_cartesian(ylim=c(0,0.25)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[3], colour=colours[3])
plot_aen_rate 

#####################
# Incidence
####################

cases <- read_excel("Input/Supplementary Data.xlsx", sheet="extended_cases", na="NA", guess_max=10^5) %>%
  mutate(
    timepoint=as_date(week), # can be changed to "week"
  )

local <- cases %>%
  group_by(timepoint) %>%
  summarise(
    n_cases = n()#/50000*100000
  ) %>%
  mutate(location="Leuven students")

belgium <- read_csv("https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv") %>% # or "Input/COVID19BE_CASES_AGESEX.csv"
  # add timepoint column
  mutate(
    month = floor_date(DATE, unit = "months"),
    week = floor_date(DATE, unit = "weeks", week_start = getOption("lubridate.week.start", 1)
    ),
    timepoint=week # can be changed to "week"
  ) %>%
  filter(
    DATE>=ymd("20210201"),
    DATE<=ymd("20220331")
  ) %>%
  group_by(timepoint) %>%
  summarise(
    n_cases = sum(CASES)/500#/11584008*100000
  ) %>%
  mutate(location="Belgium")

incidence <- bind_rows(
  local,
  belgium
)

plot_incidence <- ggplot(incidence, aes(timepoint,n_cases,fill=location, colour=location)) +
  #geom_col(position=position_dodge()) +
  geom_point() +
  geom_ribbon(aes(ymax=n_cases),ymin=0,alpha=0.4) +
  #scale_y_log10() +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "Confirmed cases (Leuven)", x = "Date", colour="Population", fill="Population") +
  scale_y_continuous(sec.axis = sec_axis(~.*500, name = "Confirmed cases (Belgium)", labels = scales::label_number(suffix = "k", scale = 1e-3))) +
  theme(
    legend.position = c(0.2,0.7)
  )
plot_incidence

#######################
# AEN vs manual tracing capacity
######################

by_indication <- observations %>%
  group_by(timepoint) %>%
  summarise(
    n_aen = sum(indication=="Digitally traced"),
    n_hrc = sum(indication=="Manually traced"),
    #n_symp = sum(indication=="symptoms")
  ) %>%
  pivot_longer(
    !timepoint,
    names_prefix = "n_",
    names_to="indication",
    values_to="n"
  )

p2 <- ggplot(by_indication, aes(x=timepoint,y=n, colour=indication)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "top") +
  scale_y_log10() +
  labs(y = "Number of tests", x = "Date")
p2

proportion_digital <- by_indication %>%
  filter(indication != "symp") %>%
  pivot_wider(names_from = "indication",values_from="n") %>%
  rowwise %>%
  mutate(
    p_digital = aen/(hrc),
    p_digital_lower = binom_CI(aen, hrc, 0, 1),
    p_digital_upper = binom_CI(aen, hrc, 0, 2)
  )

plot_proportion_digital <- ggplot(proportion_digital,aes(x=timepoint,y=p_digital,ymin=p_digital_lower,ymax=p_digital_upper)) +
  geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  geom_vline(xintercept=ymd(c("2021-04-26")), colour=colours[3]) +
  labs(y = "Ratio of DPT to MCT test indication", x = "Date") +
  coord_cartesian(ylim=c(0,.2)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[3], colour=colours[3])
plot_proportion_digital  

#################
# Contacts traced by team (per case)
# (measure of workload / yield)
###################

cases_summary <- cases %>%
  group_by(timepoint) %>%
  summarise(
    n_cases=n(),
    n_contacts_total=sum(n_contacts),
    n_contacts_per_case_median=median(n_contacts),
    n_contacts_per_case_low = quantile(n_contacts,probs=0.1),
    n_contacts_per_case_high = quantile(n_contacts,probs=0.9),
    n_at_least_1_contact = sum(n_contacts>0),
    prob_at_least_1_contact = n_at_least_1_contact / n_cases,
    prob_at_least_1_contact_low = binom_CI(n_at_least_1_contact, n_cases, 0, 1),
    prob_at_least_1_contact_high = binom_CI(n_at_least_1_contact, n_cases, 0, 2),
    sample_to_tracing_mean = mean(sample_to_tracing, na.rm=T),
    sample_to_tracing_median = median(sample_to_tracing, na.rm=T),
    sample_to_tracing_low = quantile(sample_to_tracing, na.rm=T, probs=0.1),
    sample_to_tracing_high = quantile(sample_to_tracing, na.rm=T, probs=0.9)
  )

# total contacts traced
plot_total_contacts <- ggplot(cases_summary, aes(timepoint,n_contacts_total)) +
  geom_col(width=7, fill=colours[2], colour=colours[2]) +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "MCT: total traced contacts", x = "Date")
plot_total_contacts

# tracing delay per case
plot_tracing_delay <- ggplot(filter(cases,!is.na(sample_to_tracing)), aes(timepoint,sample_to_tracing)) +
  geom_jitter(width=2, height=0.4, size=0.1, alpha=0.4, colour="darkgrey") +
  #geom_smooth(method="loess",span=loess_span,colour="red",fill="red",alpha=0.4) +
  #geom_point(data=cases_summary,aes(x=timepoint,y=sample_to_tracing_median),colour="black") +
  #geom_line(data=cases_summary,aes(x=timepoint,y=sample_to_tracing_median),colour="black") +
  #geom_errorbar(data=cases_summary,aes(x=timepoint,y=sample_to_tracing_median,ymin=sample_to_tracing_low,ymax=sample_to_tracing_high),colour="black", alpha=0.8) +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "MCT: days from case sampling to interview", x = "Date") +
  coord_cartesian(ylim=c(-0.3,12)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[2], colour=colours[2])
plot_tracing_delay

# contacts traced per case
plot_contacts_per_case <- ggplot(cases, aes(timepoint,n_contacts)) +
  geom_jitter(width=2, height=0.4, size=0.1, alpha=0.4, colour="darkgrey") +
  #geom_point(data=cases_summary,aes(x=timepoint,y=n_contacts_per_case_median),colour="black") +
  #geom_line(data=cases_summary,aes(x=timepoint,y=n_contacts_per_case_median),colour="black") +
  #geom_errorbar(data=cases_summary,aes(x=timepoint,y=n_contacts_per_case_median,ymin=n_contacts_per_case_low,ymax=n_contacts_per_case_high),colour="black", alpha=0.8) +
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "MCT: contacts per case", x = "Date") +
  coord_cartesian(ylim=c(-0.3,20)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[2], colour=colours[2])
plot_contacts_per_case

# proportion of cases with at least 1 contact traced
plot_proportion_traced <- ggplot(cases_summary,aes(x=timepoint,y=prob_at_least_1_contact,ymin=prob_at_least_1_contact_low,ymax=prob_at_least_1_contact_high)) +
  geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "MCT: proportion with ≥1 traced contact", x = "Date")  +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[2], colour=colours[2]) +
  coord_cartesian(ylim=c(0,1))
plot_proportion_traced  


#########################
# Proportion of exposures to multiple cases (not used)
########################

all_contacts <- read_excel("Input/Supplementary Data.xlsx", sheet="extended_contacts", na="NA", guess_max=10^5) %>%
  mutate(
    timepoint=as_date(week) # can be changed to "week"
  )

proportion_multi <- all_contacts %>%
  group_by(timepoint) %>%
  summarise(
    n=n(),
    n_multi = sum(additional_exposure)
  ) %>%
  rowwise %>%
  mutate(
    p_multi = n_multi/(n),
    p_multi_lower = binom_CI(n_multi, n, 0, 1),
    p_multi_upper = binom_CI(n_multi, n, 0, 2)
  )

plot_proportion_multi <- ggplot(proportion_multi,aes(x=timepoint,y=p_multi,ymin=p_multi_lower,ymax=p_multi_upper)) +
  geom_errorbar(colour="grey") +
  geom_point() +
  #geom_line()+
  theme_bw() +
  geom_vline(xintercept=ymd(c("2022-01-09","2021-10-18"))) +
  labs(y = "Proportion with multiple exposures", x = "Date") +
  coord_cartesian(ylim=c(0,0.5)) +
  geom_smooth(method = "loess", span=loess_span, se = T, linewidth=1, alpha = 0.2, level=0.95, fill=colours[2], colour=colours[2])
plot_proportion_multi 

#######################
# Combine plots
#########################

#all
((plot_incidence / plot_test_indication_absolute / plot_app_use / plot_proportion_digital / plot_aen_rate) |
   (plot_total_contacts / plot_tracing_delay / plot_contacts_per_case / plot_proportion_traced / plot_proportion_multi))  + plot_annotation(tag_levels = 'a')

#ggsave("Output/figure_time_evolution.png",width=12,height=14)

# only graphs for main figure
(plot_test_indication_absolute + plot_pos_rate_both + plot_app_use + plot_proportion_traced + plot_proportion_digital + 
    plot_tracing_delay)  + plot_annotation(tag_levels = c("a"), tag_suffix = ".") + plot_layout(ncol=2)

ggsave("Output/figure_time_evolution_main.pdf",width=12,height=11)
 
# other plots for supplementary info
(plot_incidence + plot_total_contacts + plot_aen_rate  +
    plot_contacts_per_case )  + plot_annotation(tag_levels = 'a', tag_suffix = ".") + plot_layout(ncol=2)

ggsave("Output/figure_time_evolution_supplementary.png",width=12,height=7)

