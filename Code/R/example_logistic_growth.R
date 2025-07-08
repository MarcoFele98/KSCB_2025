#_____________________________________________________________________________________________________#
#___ Noise induced order and fitting stochastic models _______________________________________________#
#___ https://royalsocietypublishing.org/doi/10.1098/rstb.2019.0381 ___________________________________#
#___ January 2024 - Marco Fele _______________________________________________________________________#
#_____________________________________________________________________________________________________#

# Set up ----
library(Rcpp)
library(ggplot2)
library(data.table)

theme_set(cowplot::theme_cowplot())
options(scipen = 10)

sourceCpp("gillespie_algorithm.cpp")


# Example logistic regression ----

# Initialize simulations
runs <- gillespie_simulation(states = c(10),
                            max_simulation_time = 25,
                            number_reagents = 1,
                            number_reactions = 3,
                            omega = 100,
                            reaction_rates = c(2, 1, 1),
                            stoichiometry_reagents = matrix(c(1,
                                                              1,
                                                              2), 
                                                            ncol = 1, 
                                                            byrow = T),
                            stoichiometry_products = matrix(c(2,
                                                              0,
                                                              1), 
                                                            ncol = 1,
                                                            byrow = T)) |>
  as.data.table() |>
  cbind(rho = 2,
        replicate = 0)

# Run stochastic simulations for different intrinsic growth rates and replicates
for(rho in seq(2, 3, by = 0.25)) {
  print(rho)
  
  for(replicate in 1:5) {
    runs <- runs |>
      rbind(gillespie_simulation(states = c(10), 
                                max_simulation_time = 25,
                                number_reagents = 1,
                                number_reactions = 3,
                                omega = 100,
                                reaction_rates = c(rho, 1, 1),
                                stoichiometry_reagents = matrix(c(1,
                                                                  1,
                                                                  2), 
                                                                ncol = 1, 
                                                                byrow = T),
                                stoichiometry_products = matrix(c(2,
                                                                  0,
                                                                  1), 
                                                                ncol = 1,
                                                                byrow = T)) |>
            as.data.table() |>
            cbind(rho = rho,
                  replicate = replicate))
    
  }
}

# Plot and save
(p <- ggplot(runs,
             aes(time, species_0, 
                 color = rho)) +
  geom_step(aes(group = interaction(replicate, rho)),
            alpha = 0.2) +
    geom_smooth(aes(group = rho),
                se = F) +
  scale_color_viridis_c(name = "\u03C1") +
  xlab("Time") +
  ylab("Number"))

ggsave("../../figures/logistic_growth.png",
       height = 3,
       width = 5,
       bg = "white")

# Example voter model ----

# Define the ODE as a function
voter_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
  dydt <- (alpha_1 - alpha_2) * y * (1 - y) - sigma_1 * y + sigma_2 * (1 - y)
  list(dydt)})
}

# Solve the ODE
out <- deSolve::ode(y = c(y = 0.01), 
                    times = seq(0, 500, by = 1), 
                    func = voter_ode, 
                    parms = c(alpha_1 = 1, 
                              alpha_2 = 1, 
                              sigma_1 = 0.01, 
                              sigma_2 = 0.01))

# Plot and save the result
(p <- ggplot(out) +
  geom_line(aes(time, y), linewidth = 1) +
  ylab("Proportion prefering option A") +
  ylim(c(0, 1)) +
  xlab("Time"))

ggsave("../../figures/voter_model_ode.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

# Initialize simulations
runs_voter <- gillespie_simulation(states = c(10, 90),
                                   max_simulation_time = 1000,
                                   number_reagents = 2,
                                   number_reactions = 4,
                                   omega = 50,
                                   reaction_rates = c(1, 1, 0.01, 0.01),
                                   stoichiometry_reagents = matrix(c(1, 1,
                                                                     1, 1,
                                                                     1, 0,
                                                                     0, 1), 
                                                                   ncol = 2, 
                                                                   byrow = T),
                                   stoichiometry_products = matrix(c(2, 0,
                                                                     0, 2,
                                                                     0, 1,
                                                                     1, 0), 
                                                                   ncol = 2,
                                                                   byrow = T))

data_l <- (runs_voter |>
  as.data.table() |>
  cbind(replicate = 1) |>
  as.data.table()
)[, ":="(duration = c(tail(waiting_time, -1), NA))]

# Run stochastic simulations for different intrinsic growth rates and replicates
for(replicate in 2:50) {
  print(replicate)
  
  data_l <- data_l |>
    rbind((gillespie_simulation(states = c(10, 90),
                               max_simulation_time = 1000,
                               number_reagents = 2,
                               number_reactions = 4,
                               omega = 50,
                               reaction_rates = c(1, 1, 0.01, 0.01),
                               stoichiometry_reagents = matrix(c(1, 1,
                                                                 1, 1,
                                                                 1, 0,
                                                                 0, 1), 
                                                               ncol = 2, 
                                                               byrow = T),
                               stoichiometry_products = matrix(c(2, 0,
                                                                 0, 2,
                                                                 0, 1,
                                                                 1, 0), 
                                                               ncol = 2,
                                                               byrow = T)) |>
            as.data.table() |>
            cbind(replicate = replicate) |>
            as.data.table()
    )[, ":="(duration = c(tail(waiting_time, -1), NA))])
  
}

# Plot and save results
(p <- ggplot(data_l[!is.na(duration) & time > 500]) +
  geom_histogram(aes(x = species_0, weight = duration, y = after_stat(density)),
                 bins = 51, fill = "black", color = "black", alpha = 0.1) +
  xlab("Number prefering option A") +
  ylab("Probability"))

ggsave("../../figures/voter_model_stochastic.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

# Non-linear opinion-dynamics model ----

# Long term behaviour
growth_decrese_nod <- data.table(x = seq(-2, 2, by = 0.01),
                               decrease_rate = seq(-2, 2, by = 0.01),
                               growth_rate = tanh(10 * seq(-2, 2, by = 0.01))) |>
  tidyr::pivot_longer(cols = contains("rate"),
                      names_pattern = "(.*)_rate",  
                      names_to = "type",
                      values_to = "rate")

(p <- ggplot(growth_decrese_nod) +
    geom_line(aes(x, rate, color = type), linewidth = 1) +
    annotate(geom = "point",
             x = -1, y = -1,
             size = 4)  +
    annotate(geom = "point",
             x = 1, y = 1,
             size = 4) +
    annotate(geom = "point",
             x = 0, y = 0,
             size = 2, shape = 21, stroke = 2,
             color = "black", fill = "white")  +
    ylab("Rate") +
    scale_color_discrete(name = "Type",
                         label = c("Decrease", "Growth")))

ggsave("../../figures/growth_decrese_nod.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

# Bifurcation diagram
equilibria_nod <- function(x, delta, alpha, u) {
  return(-delta * x + tanh(alpha * x + u))
}

parameter_space <- seq(0, 2, l = 1000)
results <- data.table(alpha = parameter_space, equ_1 = NA_real_, equ_2 = NA_real_, equ_3 = NA_real_)
for(i in 1:length(parameter_space)) { 
  alpha <- parameter_space[i]
  print(alpha)
  
  zeros <- mosaic::findZeros(equilibria_nod(x, 1, alpha, 0) ~ x, 
                     xlim = c(-1.5, 1.5))[[1]]
  results[i, alpha := alpha] 
  results[i,  1 + 1:length(zeros) := as.list(zeros)]    
}

(p <- results |>
  tidyr::pivot_longer(cols = contains("equ"),
               names_pattern = "equ_(.*)",  
               names_to = "branch",
               values_to = "value") |>
ggplot() +
  geom_point(aes(alpha, value), size = 2) +
  ylab("Equilibrium") +
  xlab("alpha"))

ggsave("../../figures/pitchfork.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

example_nod_1 <- data.table(x = seq(-2, 2, by = 0.01),
                                 decrease_rate = seq(-2, 2, by = 0.01),
                                 growth_rate = tanh(0.5 * seq(-2, 2, by = 0.01))) |>
  tidyr::pivot_longer(cols = contains("rate"),
                      names_pattern = "(.*)_rate",  
                      names_to = "type",
                      values_to = "rate")

(p <- ggplot(example_nod_1) +
    geom_line(aes(x, rate, color = type), linewidth = 1) +
    # annotate(geom = "point",
    #          x = -1, y = -1,
    #          size = 4)  +
    # annotate(geom = "point",
    #          x = 1, y = 1,
    #          size = 4) +
    # annotate(geom = "point",
    #          x = 0, y = 0,
    #          size = 2, shape = 21, stroke = 2,
    #          color = "black", fill = "white")  +
    ylab("Rate") +
    scale_color_discrete(name = "Type",
                         label = c("Decrease", "Growth")))

ggsave("../../figures/example_nod_1.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

example_nod_2 <- data.table(x = seq(-2, 2, by = 0.01),
                            decrease_rate = seq(-2, 2, by = 0.01),
                            growth_rate = tanh(1.5 * seq(-2, 2, by = 0.01))) |>
  tidyr::pivot_longer(cols = contains("rate"),
                      names_pattern = "(.*)_rate",  
                      names_to = "type",
                      values_to = "rate")

(p <- ggplot(example_nod_2) +
    geom_line(aes(x, rate, color = type), linewidth = 1) +
    # annotate(geom = "point",
    #          x = -1, y = -1,
    #          size = 4)  +
    # annotate(geom = "point",
    #          x = 1, y = 1,
    #          size = 4) +
    # annotate(geom = "point",
    #          x = 0, y = 0,
    #          size = 2, shape = 21, stroke = 2,
    #          color = "black", fill = "white")  +
    ylab("Rate") +
    scale_color_discrete(name = "Type",
                         label = c("Decrease", "Growth")))

ggsave("../../figures/example_nod_2.png",
       p,
       height = 3,
       width = 5,
       bg = "white")
