
# Each self-contained code is separated by three lines as follows.
# Please Run 1 section at a time.

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# 1: Modelling A(H1N1) No Vaccination

rm(list=ls()) 

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) No Vaccination
# 1.1 Defining parameters (transition rates in units of days^-1)
no_vaccination_params <- list(
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
isv_no_vaccination <- c(
  Sc = (no_vaccination_params$Nc),
  Vc = 0,                                      # no one is vaccinated in this scenario
  Ec = 0,                                      # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                      # the pandemic starts when an infected person is infectious
  Rc = 0,                                      # we assume no one has recovered in initial state values
  Dc = 0,                                      # we assume no one has died in isv
  
  Sa = (no_vaccination_params$Na),        
  Va = 0,
  Ea = 0,
  Ia = 0,                                      # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ra = 0,
  Da = 0
)

# 1.3 Defining timesteps   
no_vaccination_time <- seq(from = 0, to = 105, by = 1)                   # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 1.4 Defining model function
no_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 1.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 1.4.2 Describing differential equations
    dSc <- -lambda_c * Sc                                                # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- 0
    dEc <- (lambda_c * Sc) - (sigma * Ec)
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
no_vaccination_output <- as.data.frame(ode(                              # the 'ODE' function can be used to solve differential equations
  y = isv_no_vaccination,
  times = no_vaccination_time,
  func = no_vaccination_model,
  parms = no_vaccination_params
  ))

no_vaccination_output <- as.data.frame(no_vaccination_output)

# 1.6 Plotting the output
ggplot(no_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Ic, color = 'Infected Children')) +
  geom_line(aes(y = Ia, color = 'Infected Adults')) +
  geom_line(aes(y = Dc, color = 'Deaths (Children)')) +
  geom_line(aes(y = Da, color = 'Deaths (Adults)')) +
  labs(title = 'Pandemic A(H1N1) Model - No Population Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
### 2: Modelling A(H1N1) Full Vaccination

rm(list=ls()) 

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Full Vaccination
# 2.1 Defining parameters (transition rates in units of days^-1)
full_vaccination_params <- list(
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

# 2.2 Defining initial conditions (100% of children and adults vaccinated)
isv_full_vaccination <- c(
  Sc = 0,
  Vc = (full_vaccination_params$Nc),           # everyone is vaccinated in this scenario
  Ec = 0,                                      # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                      # the pandemic starts when an infected person is infectious
  Rc = 0,                                      # we assume no one has recovered in initial state values
  Dc = 0,                                      # we assume no one has died in isv
  
  Sa = 0,        
  Va = (full_vaccination_params$Na),
  Ea = 0,
  Ia = 0,                                      # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ra = 0,
  Da = 0
)

# 2.3 Defining timesteps   
full_vaccination_time <- seq(from = 0, to = 105, by = 1)                # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 2.4 Defining model function
full_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 2.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 2.4.2 Describing differential equations
    dSc <- 0                                                            # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- -((1 - nu_e) * lambda_c * Vc)
    dEc <- ((1 - nu_e) * lambda_c * Vc) - (sigma * Ec)
    dIc <- (sigma * Ec) - (gamma * Ic) - (mu_c * Ic)
    dRc <- gamma * Ic
    dDc <- mu_c * Ic
    
    dSa <- 0
    dVa <- -((1 - nu_e) * lambda_a * Va)
    dEa <- ((1 - nu_e) * lambda_a * Vc) - (sigma * Ea)
    dIa <- (sigma * Ea) - (gamma * Ia) - (mu_a * Ia)
    dRa <- gamma * Ia
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dRc, dDc, dSa, dVa, dEa, dIa, dRa, dDa)))
  })
}

# 2.5 Solving differential equations using ODE
full_vaccination_output <- as.data.frame(ode(                          # the 'ODE' function can be used to solve differential equations
  y = isv_full_vaccination,
  times = full_vaccination_time,
  func = full_vaccination_model,
  parms = full_vaccination_params
))

full_vaccination_output <- as.data.frame(full_vaccination_output)

# 2.6 Plotting the output
ggplot(full_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Ic, color = 'Infected Children')) +
  geom_line(aes(y = Ia, color = 'Infected Adults')) +
  geom_line(aes(y = Dc, color = 'Deaths (Children)')) +
  geom_line(aes(y = Da, color = 'Deaths (Adults)')) +
  labs(title = 'Pandemic A(H1N1) Model - Full Population Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
### 3: Modelling Child Vaccination Prioritization (Simplified Graphs)

rm(list=ls())  

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling Child Vaccination Prioritization (Simplified Graphs)
# 3.1 Defining parameters (transition rates in units of days^-1)
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

# 3.2 Defining initial conditions 
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

# 3.3 Defining timesteps
child_priority_vaccination_time <- seq(from = 0, to = 105, by = 1)       # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 3.4 Defining model function
child_priority_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 3.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 3.4.2 Describing differential equations
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

###-----------------------------------------------------------------------------

# 3.5 Solving differential equations using ODE
child_priority_vaccination_output <- ode(                               # the 'ODE' function can be used to solve differential equations
  y = isv_child_priority_vaccination,
  times = child_priority_vaccination_time,
  func = child_priority_vaccination_model,
  parms = child_priority_params
)

child_priority_vaccination_output <- as.data.frame(child_priority_vaccination_output)

# 3.6 Plotting the output
ggplot(child_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = 'Sc')) +
  geom_line(aes(y = Sa, color = 'Sa')) +
  geom_line(aes(y = Ic, color = 'Ic')) +
  geom_line(aes(y = Ia, color = 'Ia')) +
  geom_line(aes(y = Rc, color = 'Rc')) +
  geom_line(aes(y = Ra, color = 'Ra')) +
  geom_line(aes(y = Dc, color = 'Dc')) +
  geom_line(aes(y = Da, color = 'Da')) +
  labs(title = 'Pandemic A(H1N1) Model - Child Priority Vaccination', y = 'Population', x = 'Time (days)')
  
  ###-----------------------------------------------------------------------------

# Alternative 3.5-3.6 Plotting the output for children only
child_priority_vaccination_output <- ode(                               
  y = isv_child_priority_vaccination,
  times = child_priority_vaccination_time,
  func = child_priority_vaccination_model,
  parms = child_priority_params
)

child_priority_vaccination_output <- as.data.frame(child_priority_vaccination_output)

ggplot(child_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = 'Sc')) +
  geom_line(aes(y = Ic, color = 'Ic')) +
  geom_line(aes(y = Rc, color = 'Rc')) +
  geom_line(aes(y = Dc, color = 'Dc')) +
  labs(title = 'Pandemic A(H1N1) Model - Child Priority Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------

# Alternative 3.5-3.6 Plotting the output for adults only
child_priority_vaccination_output <- ode(                               
  y = isv_child_priority_vaccination,
  times = child_priority_vaccination_time,
  func = child_priority_vaccination_model,
  parms = child_priority_params
)

child_priority_vaccination_output <- as.data.frame(child_priority_vaccination_output)

ggplot(child_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sa, color = 'Sa')) +
  geom_line(aes(y = Ia, color = 'Ia')) +
  geom_line(aes(y = Ra, color = 'Ra')) +
  geom_line(aes(y = Da, color = 'Da')) +
  labs(title = 'Pandemic A(H1N1) Model - Child Priority Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# 4: Modelling Adult Vaccination Prioritization (Simplified Graphs)

rm(list=ls())  

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling Adult Vaccination Prioritization (Simplified Graphs)
# 4.1 Defining parameters (transition rates in units of days^-1)
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

# 4.2 Defining initial conditions 
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

# 4.3 Defining timesteps
adult_priority_vaccination_time <- seq(from = 0, to = 105, by = 1)        # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 4.4 Defining model function
adult_priority_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 4.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia))  # instead of the denominators using Nc or Na as denoted on the slides, I have used the sum of the following compartments because the denominator should reflect the current number of individuals in both groups that are currently available to be infected or contribute to infection
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic))
    
    # 4.4.2 Describing differential equations
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

###-----------------------------------------------------------------------------

# 4.5 Solving differential equations using ODE
adult_priority_vaccination_output <- ode(                               # the 'ODE' function can be used to solve differential equations
  y = isv_adult_priority_vaccination,
  times = adult_priority_vaccination_time,
  func = adult_priority_vaccination_model,
  parms = adult_priority_params
)

adult_priority_vaccination_output <- as.data.frame(adult_priority_vaccination_output)

# 4.6 Plotting the output
ggplot(adult_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = 'Sc')) +
  geom_line(aes(y = Sa, color = 'Sa')) +
  geom_line(aes(y = Ic, color = 'Ic')) +
  geom_line(aes(y = Ia, color = 'Ia')) +
  geom_line(aes(y = Rc, color = 'Rc')) +
  geom_line(aes(y = Ra, color = 'Ra')) +
  geom_line(aes(y = Dc, color = 'Dc')) +
  geom_line(aes(y = Da, color = 'Da')) +
  labs(title = 'Pandemic A(H1N1) Model - Adult Priority Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------

# Alternative 4.5-4.6 Plotting the output for children only
adult_priority_vaccination_output <- ode(                               
  y = isv_adult_priority_vaccination,
  times = adult_priority_vaccination_time,
  func = adult_priority_vaccination_model,
  parms = adult_priority_params
)

adult_priority_vaccination_output <- as.data.frame(adult_priority_vaccination_output)

ggplot(adult_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = 'Sc')) +
  geom_line(aes(y = Ic, color = 'Ic')) +
  geom_line(aes(y = Rc, color = 'Rc')) +
  geom_line(aes(y = Dc, color = 'Dc')) +
  labs(title = 'Pandemic A(H1N1) Model - Adult Priority Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------

# Alternative 4.5-4.6 Plotting the output for adults only
adult_priority_vaccination_output <- ode(                               
  y = isv_adult_priority_vaccination,
  times = adult_priority_vaccination_time,
  func = adult_priority_vaccination_model,
  parms = adult_priority_params
)

adult_priority_vaccination_output <- as.data.frame(adult_priority_vaccination_output)

ggplot(adult_priority_vaccination_output, aes(x = time)) +
  geom_line(aes(y = Sa, color = 'Sa')) +
  geom_line(aes(y = Ia, color = 'Ia')) +
  geom_line(aes(y = Ra, color = 'Ra')) +
  geom_line(aes(y = Da, color = 'Da')) +
  labs(title = 'Pandemic A(H1N1) Model - Adult Priority Vaccination', y = 'Population', x = 'Time (days)')

###-----------------------------------------------------------------------------
###----------------------------------------------------------------------------
###----------------------------------------------------------------------------
# 5: Modelling A(H1N1) Child Vaccination Prioritization with Antiviral Treatment

rm(list=ls()) 

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Child Vaccination Prioritization with Antiviral Treatment
# 5.1 Defining parameters (transition rates in units of days^-1)
child_priority_params_treatment <- list(
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
  N_total = 68300000,   # the total UK population (c + a)
  
  # new parameters
  gamma_T = 0.7,             # the recovery rate for treated individuals, which is 3.5x the normal recovery rate
  h = 4,                     # the rate at which antiviral treatment is administered
  target_vaccination = 0.95  # 95% of children will eventually be vaccinated (this is 20% of the total UK population)
)

# 5.2 Defining initial conditions
isv_child_priority_vaccination_treatment <- c(
  Sc = (child_priority_params_treatment$Nc - 1) * 0.05,  # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 95% of the child population; this is because 20% of 68,300,000 is 13,660,000, which is 95% of 14,400,000
  Vc = (child_priority_params_treatment$Nc - 1) * 0.95,  # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 95% of the child population; this is because 20% of 68,300,000 is 13,660,000, whihc is 95% of 14,400,000
  Ec = 0,                                                # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                                # the pandemic starts when an infected person is infectious
  Tc = 0,                                                # we assume no one has been treated
  Rc = 0,                                                # we assume no one has recovered in initial state values
  Dc = 0,                                                # we assume no one has died in isv
  
  Sa = (child_priority_params_treatment$Na * 1),         # because we are prioritizing children, all 20% of the allocatable vaccine is to be given to the children                       # we make the same assumtions for adults as we did for children, except for Ia
  Va = 0,
  Ea = 0,
  Ia = 0,                                                # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ta = 0,
  Ra = 0,
  Da = 0
)

# 5.3 Defining timesteps
child_priority_vaccination_treatment_time <- seq(from = 0, to = 50, by = 1)    # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 5.4 Defining model function
child_priority_vaccination_treatment_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 5.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia)) + (beta_cc * Tc / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ta / sum(Sa, Va, Ea, Ia))     # we must account for how antiviral treatment affects transmission dynamics, and the eliminate the assumption of an undynamic population
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_aa * Ta / sum(Sa, Va, Ea, Ia)) + (beta_ac * Tc / sum(Sc, Vc, Ec, Ic))
    
    # 5.4.2 Describing differential equations
    dSc <- -lambda_c * Sc                                                      # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- -((1 - nu_e) * lambda_c * Vc)
    dEc <- (lambda_c * Sc) + ((1 - nu_e) * lambda_c * Vc) - (sigma * Ec)
    dIc <- (sigma * Ec) - ((0.65 * h) * Ic) - ((1 - 0.65) * gamma * Ic) - (mu_c * Ic)
    dTc <- ((0.65 * h) * Ic) - (gamma_T * Tc)
    dRc <- ((1 - 0.65) * gamma * Ic) + (gamma_T * Tc)
    dDc <- mu_c * Ic
    
    dSa <- -lambda_a * Sa
    dVa <- 0
    dEa <- (lambda_a * Sa) - (sigma * Ea)
    dIa <- (sigma * Ea) - ((0.65 * h) * Ia) - ((1 - 0.65) * gamma * Ia) - (mu_a * Ia)
    dTa <- ((0.65 * h) * Ia) - (gamma_T * Ta)
    dRa <- ((1 - 0.65) * gamma * Ia) + (gamma_T * Ta)
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dTc, dRc, dDc, dSa, dVa, dEa, dIa, dTa, dRa, dDa)))
  })
}

###-----------------------------------------------------------------------------

# Showing children only
# 5.5 Solving differential equations using ODE
child_priority_vaccination_treatment_output <- ode(                            # the 'ODE' function can be used to solve differential equations
  y = isv_child_priority_vaccination_treatment,
  times = child_priority_vaccination_treatment_time,
  func = child_priority_vaccination_treatment_model,
  parms = child_priority_params_treatment,
  method = "ode45"
)

child_priority_vaccination_treatment_output <- as.data.frame(child_priority_vaccination_treatment_output)

# 5.6 Plotting the output
ggplot(child_priority_vaccination_treatment_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = "Sc")) +
  geom_line(aes(y = Vc, color = "Vc")) +
  geom_line(aes(y = Ec, color = "Ec")) +
  geom_line(aes(y = Ic, color = "Ic")) +
  geom_line(aes(y = Tc, color = "Tc")) +
  geom_line(aes(y = Rc, color = "Rc")) +
  geom_line(aes(y = Dc, color = "Dc")) +
  labs(title = "Pandemic A(H1N1) Model - Child Priority Vaccination with Treatment", x = "Time (days)", y = "Population")

###-----------------------------------------------------------------------------

# Alternative 5.5-5.6 showing adults only
child_priority_vaccination_treatment_output <- ode(                            # the 'ODE' function can be used to solve differential equations
  y = isv_child_priority_vaccination_treatment,
  times = child_priority_vaccination_treatment_time,
  func = child_priority_vaccination_treatment_model,
  parms = child_priority_params_treatment,
  method = "ode45"
)

child_priority_vaccination_treatment_output <- as.data.frame(child_priority_vaccination_treatment_output)

ggplot(child_priority_vaccination_treatment_output, aes(x = time)) +
  geom_line(aes(y = Sa, color = "Sa")) +
  geom_line(aes(y = Va, color = "Va")) +
  geom_line(aes(y = Ea, color = "Ea")) +
  geom_line(aes(y = Ia, color = "Ia")) +
  geom_line(aes(y = Ta, color = "Ta")) +
  geom_line(aes(y = Ra, color = "Ra")) +
  geom_line(aes(y = Da, color = "Da")) +
  labs(title = "Pandemic A(H1N1) Model - Child Priority Vaccination with Treatment", x = "Time (days)", y = "Population")

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# 6: Modelling A(H1N1) Adult Vaccination Prioritization with Antiviral Treatment

rm(list=ls()) 

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Adult Vaccination Prioritization with Antiviral Treatment
# 6.1 Defining parameters (transition rates in units of days^-1)
adult_priority_params_treatment <- list(
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
  N_total = 68300000,  # the total UK population (c + a)
  
  # new parameters
  gamma_T = 0.7,             # the recovery rate for treated individuals, which is 3.5x the normal recovery rate
  h = 4,                     # the rate at which antiviral treatment is administered
  target_vaccination = 0.95  # 95% of children will eventually be vaccinated (this is 20% of the total UK population)
)

# 6.2 Defining initial conditions 
isv_adult_priority_vaccination_treatment <- c(
  Sc = (adult_priority_params_treatment$Nc - 1),         # the research question specifies we only have enough vaccine to cover 20% of the population, which is equivalent to 25% of the adult population; We won't vaccinate children here
  Vc = 0,                                                
  Ec = 0,                                                # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                                # the pandemic starts when an infected person is infectious
  Tc = 0,                                                # we assume no one has been treated
  Rc = 0,                                                # we assume no one has recovered in initial state values
  Dc = 0,                                                # we assume no one has died in isv
  
  Sa = (adult_priority_params_treatment$Na * 0.75),      # because we are prioritizing adults, all 20% of the allocatable vaccine is to be given to the adults                       # we make the same assumtions for adults as we did for children, except for Ia
  Va = (adult_priority_params_treatment$Na * 0.25),
  Ea = 0,
  Ia = 0,                                                # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ta = 0,
  Ra = 0,
  Da = 0
)

# 6.3 Defining timesteps
adult_priority_vaccination_treatment_time <- seq(from = 0, to = 50, by = 1)    # specifying timesteps given the average duration of a A(H1N1) pandemic is 15 weeks (which is, 105 days)

# 6.4 Defining model function
adult_priority_vaccination_treatment_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 6.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia)) + (beta_cc * Tc / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ta / sum(Sa, Va, Ea, Ia))     # we must account for how antiviral treatment affects transmission dynamics, and the eliminate the assumption of an undynamic population
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_aa * Ta / sum(Sa, Va, Ea, Ia)) + (beta_ac * Tc / sum(Sc, Vc, Ec, Ic))
    
    # 6.4.2 Describing differential equations
    dSc <- -lambda_c * Sc                                                      # differential equations can be used to describe model compartments as visualized in the presentation slides
    dVc <- 0
    dEc <- (lambda_c * Sc) - (sigma * Ec)
    dIc <- (sigma * Ec) - ((0.65 * h) * Ic) - ((1 - 0.65) * gamma * Ic) - (mu_c * Ic)
    dTc <- ((0.65 * h) * Ic) - (gamma_T * Tc)
    dRc <- ((1 - 0.65) * gamma * Ic) + (gamma_T * Tc)
    dDc <- mu_c * Ic
    
    dSa <- -lambda_a * Sa
    dVa <- -((1 - nu_e) * lambda_a * Va)
    dEa <- (lambda_c * Sa) + ((1 - nu_e) * lambda_a * Va) - (sigma * Ea)
    dIa <- (sigma * Ea) - ((0.65 * h) * Ia) - ((1 - 0.65) * gamma * Ia) - (mu_a * Ia)
    dTa <- ((0.65 * h) * Ia) - (gamma_T * Ta)
    dRa <- ((1 - 0.65) * gamma * Ia) + (gamma_T * Ta)
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dTc, dRc, dDc, dSa, dVa, dEa, dIa, dTa, dRa, dDa)))
  })
}

###-----------------------------------------------------------------------------

# Showing children only
# 6.5 Solving differential equations using ODE
adult_priority_vaccination_treatment_output <- ode(                            # the 'ODE' function can be used to solve differential equations
  y = isv_adult_priority_vaccination_treatment,
  times = adult_priority_vaccination_treatment_time,
  func = adult_priority_vaccination_treatment_model,
  parms = adult_priority_params_treatment,
  method = "ode45"
)

adult_priority_vaccination_treatment_output <- as.data.frame(adult_priority_vaccination_treatment_output)

# 6.6 Plotting the output
ggplot(adult_priority_vaccination_treatment_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = "Sc")) +
  geom_line(aes(y = Vc, color = "Vc")) +
  geom_line(aes(y = Ec, color = "Ec")) +
  geom_line(aes(y = Ic, color = "Ic")) +
  geom_line(aes(y = Tc, color = "Tc")) +
  geom_line(aes(y = Rc, color = "Rc")) +
  geom_line(aes(y = Dc, color = "Dc")) +
  labs(title = "Pandemic A(H1N1) Model - Adult Priority Vaccination with Treatment", x = "Time (days)", y = "Population")

###-----------------------------------------------------------------------------

# Alternative 6.5-6.6 showing adults only
adult_priority_vaccination_treatment_output <- ode(                            # the 'ODE' function can be used to solve differential equations
  y = isv_adult_priority_vaccination_treatment,
  times = adult_priority_vaccination_treatment_time,
  func = adult_priority_vaccination_treatment_model,
  parms = adult_priority_params_treatment,
  method = "ode45"
)

adult_priority_vaccination_treatment_output <- as.data.frame(adult_priority_vaccination_treatment_output)

ggplot(adult_priority_vaccination_treatment_output, aes(x = time)) +
  geom_line(aes(y = Sa, color = "Sa")) +
  geom_line(aes(y = Va, color = "Va")) +
  geom_line(aes(y = Ea, color = "Ea")) +
  geom_line(aes(y = Ia, color = "Ia")) +
  geom_line(aes(y = Ta, color = "Ta")) +
  geom_line(aes(y = Ra, color = "Ra")) +
  geom_line(aes(y = Da, color = "Da")) +
  labs(title = "Pandemic A(H1N1) Model - Adult Priority Vaccination with Treatment", x = "Time (days)", y = "Population")

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# 7: Modelling A(H1N1) Child Vaccination Prioritization with Antiviral Treatment and Vaccination Rate

rm(list=ls())   

# Loading packages
library(deSolve)
library(reshape2)
library(ggplot2)

# Modelling A(H1N1) Child Vaccination Prioritization with Antiviral Treatment and Vaccination Rate
# 7.1 Defining parameters (transition rates in units of days^-1)
child_priority_params_T <- list(
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
  N_total = 68300000,  # the total UK population (c + a)
  
  # new parameters
  nu_r = 0.1,                # the rate at which vaccination is administered (10%, though an unrealistic figure)
  gamma_T = 0.7,             # the recovery rate for treated individuals, which is 3.5x the normal recovery rate
  h = 4,                     # the rate at which antiviral treatment is administered
  target_vaccination = 0.95  # 95% of children will eventually be vaccinated (this is 20% of the total UK population)
)

# 7.2 Defining initial conditions 
isv_child_priority_vaccination_T <- c(
  Sc = (child_priority_params_T$Nc - 1),       # in this scenario, we won't assume that 95% of children start off in compartment Vc; they must progress there over rate nu_r, given we can't vaccinate everyone all at once
  Vc = 0,                                      # the same logic follows here
  Ec = 0,                                      # due to A(H1N1)'s latency period, the pandemic does not start until an infected person is infectious
  Ic = 1,                                      # the pandemic starts when an infected person is infectious
  Tc = 0,                                      # we assume no one has been treated yet
  Rc = 0,                                      # we assume no one has recovered in initial state values
  Dc = 0,                                      # we assume no one has died in isv
  
  Sa = (child_priority_params_T$Na * 1),       
  Va = 0,
  Ea = 0,
  Ia = 0,                                      # we assume the pandemic starts when a child is infected and infectious, not an adult
  Ta = 0,
  Ra = 0,
  Da = 0
)

# 7.3 Defining timesteps
child_priority_vaccination_T_time <- seq(from = 0, to = 35, by = 1)            # i have reduced the timespan of the model to better visualize the plots

# 7.4 Defining model function
child_priority_vaccination_T_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 7.4.1 Defining force of infection (for 'c' and 'a')
    lambda_c <- (beta_cc * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ia / sum(Sa, Va, Ea, Ia)) + (beta_cc * Tc / sum(Sc, Vc, Ec, Ic)) + (beta_ca * Ta / sum(Sa, Va, Ea, Ia))     # we must account for how antiviral treatment affects transmission dynamics, and the eliminate the assumption of an undynamic population
    lambda_a <- (beta_aa * Ia / sum(Sa, Va, Ea, Ia)) + (beta_ac * Ic / sum(Sc, Vc, Ec, Ic)) + (beta_aa * Ta / sum(Sa, Va, Ea, Ia)) + (beta_ac * Tc / sum(Sc, Vc, Ec, Ic))
    
    # Vaccination process (gradually moving towards 95% child population coverage)
    nu_r_dynamic <- ifelse(Vc < target_vaccination * Nc, nu_r, 0)
    
    # 7.4.2 Describing differential equations
    dSc <- (-lambda_c * Sc) - (nu_r_dynamic * Sc)
    dVc <- (nu_r_dynamic * Sc) - ((1 - nu_e) * lambda_c * Vc)
    dEc <- (lambda_c * Sc) + ((1 - nu_e) * lambda_c * Vc) - (sigma * Ec)
    dIc <- (sigma * Ec) - ((0.65 * h) * Ic) - ((1 - 0.65) * gamma * Ic) - (mu_c * Ic)
    dTc <- ((0.65 * h) * Ic) - (gamma_T * Tc)
    dRc <- ((1 - 0.65) * gamma * Ic) + (gamma_T * Tc)
    dDc <- mu_c * Ic
    
    dSa <- (-lambda_a * Sa)
    dVa <- 0
    dEa <- (lambda_a * Sa) - (sigma * Ea)
    dIa <- (sigma * Ea) - ((0.65 * h) * Ia) - ((1 - 0.65) * gamma * Ia) - (mu_a * Ia)
    dTa <- ((0.65 * h) * Ia) - (gamma_T * Ta)
    dRa <- ((1 - 0.65) * gamma * Ia) + (gamma_T * Ta)
    dDa <- mu_a * Ia
    
    return(list(c(dSc, dVc, dEc, dIc, dTc, dRc, dDc, dSa, dVa, dEa, dIa, dTa, dRa, dDa)))
  })
}

# 7.5 Solving differential equations using ODE
child_priority_vaccination_T_output <- ode(
  y = isv_child_priority_vaccination_T,
  times = child_priority_vaccination_T_time,
  func = child_priority_vaccination_T_model,
  parms = child_priority_params_T,
  method = "ode45"
)


child_priority_vaccination_T_output <- as.data.frame(child_priority_vaccination_T_output)

# 7.6 Plotting the output
ggplot(child_priority_vaccination_T_output, aes(x = time)) +
  geom_line(aes(y = Sc, color = "Sc")) +
  geom_line(aes(y = Vc, color = "Vc")) +
  geom_line(aes(y = Ec, color = "Ec")) +
  geom_line(aes(y = Ic, color = "Ic")) +
  geom_line(aes(y = Tc, color = "Tc")) +
  geom_line(aes(y = Rc, color = "Rc")) +
  geom_line(aes(y = Dc, color = "Dc")) +
  labs(title = "Pandemic A(H1N1) Model - Child Priority Vaccination with Treatment", x = "Time (days)", y = "Population")

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
# END


