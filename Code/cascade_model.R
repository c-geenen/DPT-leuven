library(tidyverse)
library(fitdistrplus)
library(gridExtra)
#library(metR)
library(ggtext)
library(readxl)
library(ggrepel)

set.seed(1) 

p <- tibble(p=numeric())   # Initialise results table
cascade <- read_csv("output/cascade_summary.csv")

# In our stochastic model, we use a beta distribution for the uncertain success probability of each step
# in the cascade, based the numbers of successes and failures from our results.
# alpha minus 1 is the number of successes and beta minus 1 the number of failures:
step1q <- function(x) { qbeta(x, cascade$successes[1]+1, cascade$total[1]-cascade$successes[1]+1) }
step2q <- function(x) { qbeta(x, cascade$successes[2]+1, cascade$total[2]-cascade$successes[2]+1 ) }
step3q <- function(x) { qbeta(x, cascade$successes[3]+1, cascade$total[3]-cascade$successes[3]+1) }
step4q <- function(x) { qbeta(x, cascade$successes[4]+1, cascade$total[4]-cascade$successes[4]+1 ) }

cascade_completion_sample <- function() {
  # use a random number between 0 and 1 to sample from the probability distribution of
  # the actual probability to progress through step 1
  p1 <- step1q( runif(1) )
  # Repeat for the other steps
  p2 <- step2q( runif(1) )
  p3 <- step3q( runif(1) )
  p4 <- step4q( runif(1) )
  
  # Multiply these conditional probabilities to determine probability to complete entire cascade
  return(p1*p2*p3*p4)
}

for (i in 1:100000) {
  if (i %% 10000==0) print(paste0(i/1000, "k iterations"))
  
  p <- p %>%
    add_row(p = cascade_completion_sample())
}

#write_csv(p, "Output/model_iterations.csv")
#p <- read_csv("Output/model_iterations.csv")

# Fit a beta distribution to the data
(fit_params <- fitdist(p$p,"beta"))

# Median and 95% confidence interval for the probability of completing the cascade.
qbeta(c(0.025,0.5,0.975), fit_params$estimate[1], fit_params$estimate[2])
# Mean (estimate) of the probability of completing the cascade
(estimate_completion <- fit_params$estimate[1] / (fit_params$estimate[1]+fit_params$estimate[2]))

# Plot the density functions
# alpha minus 1 is the number of successes and beta minus 1 the numbers of failures
step1d <- function(x) { dbeta(x, cascade$successes[1]+1, cascade$total[1]-cascade$successes[1]+1) }
step2d <- function(x) { dbeta(x, cascade$successes[2]+1, cascade$total[2]-cascade$successes[2]+1 ) }
step3d <- function(x) { dbeta(x, cascade$successes[3]+1, cascade$total[3]-cascade$successes[3]+1) }
step4d <- function(x) { dbeta(x, cascade$successes[4]+1, cascade$total[4]-cascade$successes[4]+1 ) }
fitted <- function(x) { dbeta(x, fit_params$estimate[1], fit_params$estimate[2]) }

plot1 <- ggplot(p, aes(p)) + 
  geom_histogram(binwidth=0.0001, aes(y=after_stat(density)), alpha=0.5) +
  geom_function(fun=fitted, aes(colour="Entire cascade"), n=1000) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,0.08)

plot2 <- ggplot() + 
  geom_function(fun=fitted, aes(colour="Entire cascade"), n=1000) +
  geom_function(fun=step1d, aes(colour="Step 1: app use by case"), n=1000) +
  geom_function(fun=step2d, aes(colour="Step 2: identifier upload"), n=1000) +
  geom_function(fun=step3d, aes(colour="Step 3: app use by contact"), n=1000) +
  geom_function(fun=step4d, aes(colour="Step 4: receipt of notification"), n=1000) +
  theme_bw() +
  xlim(0,0.7) +
  labs(
    colour = "Step of cascade",
    x = "Probability of success",
    y = "Density"
  )

plot <- grid.arrange(plot1, plot2)
plot2
ggsave("Output/cascade_completion.png", plot2, width=8, height=4)


########################
# Sensitivity of MCT
########################

# DPT only
(DPT_ci <- c(
  qbeta(c(0.025), fit_params$estimate[1], fit_params$estimate[2]),
  estimate_completion,
  qbeta(c(0.975), fit_params$estimate[1], fit_params$estimate[2])
))

# MCT only
(MCT_ci <- DPT_ci * 616/60)

# DPT+MCT
(combo_ci <- DPT_ci * 648/60)


###############################
# Different app uptake (not used)
###############################

app_uptake <- seq(0.01,1,by=0.01)
result <- tibble(app_uptake=numeric(),p=numeric(),upper=numeric(), lower=numeric())

for (a in app_uptake) {
  p_pop <- tibble(p=numeric())   # Initialise results table
  
  for (i in 1:1000) {
    if (i %% 10000==0) print(paste0("App uptake ",a,", ",i/1000, "k iterations"))
    
    p1 <- a
    p2 <- step2q( runif(1) )
    p3 <- a
    p4 <- step4q( runif(1) )
    
    p_pop <- p_pop %>%
      add_row(p = p1*p2*p3*p4)
  }
  
  fit <- fitdist(p_pop$p,"beta")
  
  result <- result %>%
    add_row(
      app_uptake = a,
      p = qbeta(0.5, fit$estimate[1], fit$estimate[2]),
      lower = qbeta(0.025, fit$estimate[1], fit$estimate[2]),
      upper = qbeta(0.975, fit$estimate[1], fit$estimate[2]) 
    )

}

ggplot(result, aes(app_uptake,p)) +
  geom_smooth() +
  geom_point() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.30))

ggsave("output/cascade_completion_by_app_uptake.png",width=4, height=4)


###############################
# Different app uptake and AEN triggering
###############################

our_app_uptake = mean(c(
  cascade$successes[1] / cascade$total[1], # case app uptake
  cascade$successes[3] / cascade$total[3]  # contact app uptake
))
our_AEN_trigger_prob = cascade$successes[2] / cascade$total[2]

app_uptake <- seq(0,1,by=0.005)
AEN_trigger_prob <- seq(0,1,by=0.005)
result_b <- tibble(app_uptake=numeric(),AEN_trigger_prob=numeric(),p=numeric())

set.seed(1) 

for (a in app_uptake) {
  print(paste0("App uptake ",a,", "))
  
  for (b in AEN_trigger_prob) {
  
    p_pop <- tibble(p=numeric())   # Initialise results table
    
    p1 <- a
    p2 <- b
    p3 <- a
    p4 <- cascade$successes[4] / cascade$total[4]
    
    result_b <- result_b %>%
      add_row(
        app_uptake = a,
        AEN_trigger_prob = b,
        p = p1 * p2 * p3 * p4
      )
    
  }
}

#write_csv(result_b, "Output/model_iterations_b.csv")
#result_b <- read_csv("Output/model_iterations_b.csv")

other_studies <- read_csv("input/selected_studies.csv") %>%
  mutate(
    `DPT system`=factor(replace_na(`DPT system`,"Other"), levels=c("Coronalert","NHS COVID-19 app","SwissCovid","Other")),
    use_label = (`Authors`=="Ballouz et al"),
    Authors = str_replace_all(Authors," - ","\n")
  )

ggplot(result_b, aes(app_uptake,AEN_trigger_prob)) +
  geom_raster(aes(fill=p)) +
  geom_contour(binwidth=0.01, colour="black", aes(z=p), alpha=0.1) +
  geom_contour(binwidth=0.1, colour="black", aes(z=p)) +
  scale_fill_gradient2(low = "white", high = "#ffc000") +
  metR::geom_text_contour(aes(z = p), rotate=F, binwidth=0.1, skip=0, stroke=0.2, nudge_x=-0.03, nudge_y=0.01, colour="black") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 11, hjust = 0.5)
  ) +
  coord_fixed(ylim=c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    fill = "Probability\n of cascade\n completion",
    y = "Probability of triggering notifications",
    x = "App uptake"
  ) +
  #geom_point(size=5, colour="white", shape=15, data=filter(other_studies,use_label), aes(x=`App uptake`, y=`Infected app users consenting to contact tracing`)) +
  geom_point(size=3, data=other_studies, aes(x=`App uptake`, y=`Infected app users consenting to contact tracing`,colour=`DPT system`, shape=`DPT system`)) + 
  geom_text_repel(point.padding=2, size=2.5, data=other_studies, aes(x=`App uptake`, y=`Infected app users consenting to contact tracing`,label=`Authors`, colour=`DPT system`), show.legend=F) +
  #geom_label_repel(point.padding=5, size=3, data=filter(other_studies,use_label), aes(x=`App uptake`, y=`Infected app users consenting to contact tracing`,label=`Authors`, colour=`DPT system`), show.legend=F) +
  scale_shape_manual(values=c(15,16,17,18))
  

ggsave("output/cascade_completion_by_app_uptake_and_trigger_rate.pdf",width=7, height=4.8)
