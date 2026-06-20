# Each self-contained code is separated by three lines as follows.
# Please Run 1 section at a time.

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# Scenario 1: Modelling A(H1N1) Child Vaccination Prioritization

rm(list=ls()) 

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Child Vaccination Prioritization
# 1.1 Defining parameters (transition rates in units of days^-1)
child_priority_params <- list(
  beta_cc = 4.4,       # the infection rate between children
  beta_aa = 2.2,       # the infection rate between adults
  beta_ca = 0.88,      # the infection rate from adults to children 
  beta_ac = 1.32,      # the infection rate from children to adults
  nu_e = 0.6,          # the vaccine effectiveness rate
  sigma = 0.625,       # the progression rate from exposed to infectious state, which is the same for children and adults
  gamma = 0.2,         # the recovery rate, which is the same for children and adults
  mu_c = 0.001,        # the death rate for children
  mu_a = 0.05,         # the death rate for adults
  Nc = 14400000,       # the UK child population (c)
  Na = 53900000,       # the UK adult population (a)
  N_total = 68300000   # the total UK population (c + a)
)

# 1.2 Defining initial conditions
isv_child_priority_vaccination <- c(
  Sc = (child_priority_params$Nc - 1) * 0.05,  # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 95% of the child population; this is because 20% of 68,300,000 is 13,660,000, which is 95% of 14,400,000
  Vc = (child_priority_params$Nc - 1) * 0.95,  # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 95% of the child population; this is because 20% of 68,300,000 is 13,660,000, whihc is 95% of 14,400,000
  Ec = 0,                                      # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                      # the pandemic starts when an infected person is infectious
  Rc = 0,                                      # we assume no one has recovered in initial state values
  Dc = 0,                                      # we assume no one has died in isv
  
  Sa = (child_priority_params$Na * 1),         # because we are prioritizing children, all 20% of the allocatable vaccine is to be given to the children                       # we make the same assumtions for adults as we did for children, except for Ia
  Va = 0,
  Ea = 0,
  Ia = 0,                                      # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ra = 0,
  Da = 0
)

# 1.3 Defining timesteps
child_priority_vaccination_time <- seq(from = 0, to = 105, by = 1)        # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 1.4 Defining model function
child_priority_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 1.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 1.4.2 Describing differential equations
    dSc <- -lambda_c * Sc                                                # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- -((1 - nu_e) * lambda_c * Vc)
    dEc <- (lambda_c * Sc) + ((1 - nu_e) * lambda_c * Vc) - (sigma * Ec)
    dIc <- (sigma * Ec) - (gamma * Ic) - (mu_c * Ic)
    dRc <- gamma * Ic
    dDc <- mu_c * Ic
    
    dSa <- -lambda_a * Sa
    dVa <- 0
    dEa <- (lambda_a * Sa) - (sigma * Ea)
    dIa <- (sigma * Ea) - (gamma * Ia) - (mu_a * Ia)
    dRa <- gamma * Ia
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dRc, dDc, dSa, dVa, dEa, dIa, dRa, dDa)))
  })
}

# 1.5 Solving differential equations using ODE
child_priority_vaccination_output <- as.data.frame(ode(                  # the 'ODE' function can be used to solve differential equations
  y = isv_child_priority_vaccination,
  times = child_priority_vaccination_time,
  func = child_priority_vaccination_model,
  parms = child_priority_params
  ))

child_priority_vaccination_output_long <- melt(as.data.frame(child_priority_vaccination_output), id = "time")

child_priority_vaccination_output_long$proportion <- child_priority_vaccination_output_long$value/sum(isv_child_priority_vaccination)


# 1.6 Plotting the output
ggplot(data = child_priority_vaccination_output_long,
       aes(x = time, y = proportion, colour = variable, group = variable)) +
  geom_line() +
  xlab('Time (days)') +
  ylab('Proportion of the Population') +
  labs(title = paste('Pandemic A(H1N1) Model - Child Priority Vaccination'), colour = 'Compartment')

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# Scenario 2: Modelling A(H1N1) Adult Vaccination Prioritization

rm(list=ls())  

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Adult Vaccination Prioritization
# 2.1 Defining parameters (transition rates in units of days^-1)
adult_priority_params <- list(
  beta_cc = 4.4,       # the infection rate between children
  beta_aa = 2.2,       # the infection rate between adults
  beta_ca = 0.88,      # the infection rate from adults to children 
  beta_ac = 1.32,      # the infection rate from children to adults
  nu_e = 0.6,          # the vaccine effectiveness rate
  sigma = 0.625,       # the progression rate from exposed to infectious state, which is the same for children and adults
  gamma = 0.2,         # the recovery rate, which is the same for children and adults
  mu_c = 0.001,        # the death rate for children
  mu_a = 0.05,         # the death rate for adults
  Nc = 14400000,       # the UK child population (c)
  Na = 53900000,       # the UK adult population (a)
  N_total = 68300000   # the total UK population (c + a)
)

# 2.2 Defining initial conditions
isv_adult_priority_vaccination <- c(
  Sc = (adult_priority_params$Nc - 1),         # the '-1' represents the one child in the 'Ic' compartment; because we are prioritizing adults, all 20% of the allocatable vaccine is to be given to the adults 
  Vc = 0,                                      # no children are vaccinated in this scenario
  Ec = 0,                                      # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                      # the pandemic starts when an infected person is infectious
  Rc = 0,                                      # we assume no one has recovered in initial state values
  Dc = 0,                                      # we assume no one has died in isv
  
  Sa = (adult_priority_params$Na * 0.75),      # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 25% of the adult population; this is because 20% of 68,300,000 is 13,660,000, which is 25% of 59,300,000                       # we make the same assumtions for adults as we did for children, except for Ia
  Va = (adult_priority_params$Na * 0.25),      # this same logic follows
  Ea = 0,
  Ia = 0,                                      # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ra = 0,
  Da = 0
)

# 2.3 Defining timesteps
adult_priority_vaccination_time <- seq(from = 0, to = 105, by = 1)       # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 2.4 Defining model function
adult_priority_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 2.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 2.4.2 Describing differential equations
    dSc <- -lambda_c * Sc                                                # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- 0
    dEc <- (lambda_c * Sc) - (sigma * Ec)
    dIc <- (sigma * Ec) - (gamma * Ic) - (mu_c * Ic)
    dRc <- gamma * Ic
    dDc <- mu_c * Ic
    
    dSa <- -lambda_a * Sa
    dVa <- -((1 - nu_e) * lambda_a * Va)
    dEa <- (lambda_a * Sa) + ((1 - nu_e) * lambda_a * Va) - (sigma * Ea)
    dIa <- (sigma * Ea) - (gamma * Ia) - (mu_a * Ia)
    dRa <- gamma * Ia
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dRc, dDc, dSa, dVa, dEa, dIa, dRa, dDa)))
  })
}

# 2.5 Solving differential equations using ODE
adult_priority_vaccination_output <- as.data.frame(ode(                 # the 'ODE' function can be used to solve differential equations
  y = isv_adult_priority_vaccination,
  times = adult_priority_vaccination_time,
  func = adult_priority_vaccination_model,
  parms = adult_priority_params
))

adult_priority_vaccination_output_long <- melt(as.data.frame(adult_priority_vaccination_output), id = "time")

adult_priority_vaccination_output_long$proportion <- adult_priority_vaccination_output_long$value/sum(isv_adult_priority_vaccination)


# 2.6 Plotting the output
ggplot(data = adult_priority_vaccination_output_long,
       aes(x = time, y = proportion, colour = variable, group = variable)) +
  geom_line() +
  xlab('Time (days)') +
  ylab('Proportion of the Population') +
  labs(title = paste('Pandemic A(H1N1) Model - Adult Priority Vaccination'), colour = 'Compartment')

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# END


