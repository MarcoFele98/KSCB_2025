#_____________________________________________________________________________________________________#
#___ Noise induced order and fitting stochastic models _______________________________________________#
#___ https://royalsocietypublishing.org/doi/10.1098/rstb.2019.0381 ___________________________________#
#___ January 2024 - Marco Fele _______________________________________________________________________#
#_____________________________________________________________________________________________________#

# Set up ----
library(Rcpp) # In first part of the tutorial, I simulate two models with the Gillespie algorithm implemented in Rcpp. If Rcpp is not set up, don't worry because it is not needed anywhere else 
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

# Time series 
# Define the ODE as a function
nod_ode <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    dydt <- -delta * y + tanh(alpha * y + u)
    list(dydt)})
}

# Solve the ODE
out_nod <- data.table(time = numeric(), y = numeric(), alpha = numeric(), start = numeric())
for(alpha in c(0.5, 1.5)) {
  for(start in seq(-1.2, 1.2, by = 0.1)) {
    out_nod <- out_nod |>
      rbind((deSolve::ode(y = c(y = start), 
                          times = seq(0, 10, by = 0.1), 
                          func = nod_ode, 
                          parms = c(delta = 1, 
                                    alpha = alpha, 
                                    u = 0)) |>
               as.data.table()
      )[, ":="(alpha = alpha,
               start = start)])
  }
}

# Plot and save the result
(p <- ggplot(out_nod) +
    geom_line(aes(time, y, color = start, group = start), linewidth = 1) +
    scale_color_viridis_c(option = "mako") +
    ylab("Decision") +
    facet_wrap(~alpha) +
    guides(color = "none") +
    xlab("Time"))

ggsave("../../figures/nod_time_series.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

# Elementary bifurcation ----
## Saddle node ----
saddle_node <-  data.table(
  x = seq(-2, 2, by = 0.01)
)[, ":="(p = x^2,
         is_stable = ifelse(x < 0, "yes", "no"),
         eigenvalues = 2 * x)
  ][order(-p)
    ][, ":="(appear = rleid(as.factor(round(p, 4))))]

(p <- ggplot(saddle_node) +
  geom_path(aes(p, x, linetype = is_stable), linewidth = 1) +
  xlim(c(-4, 4)) +
  scale_linetype_manual (values = c("dashed", "solid")) +
    xlab("Parameter p") +
    ylab("Equilibria x*") +
  guides(linetype = "none"))

ggsave("../../figures/saddle_node.png",
       p,
       height = 3,
       width = 5,
       bg = "white")


animation <- data.table(
  appear = max(saddle_node$appear)+1:max(saddle_node$appear)+25,
  no = NA_real_,
  yes = NA_real_) |>
  rbind(saddle_node[, c("eigenvalues", "is_stable", "appear")] |>
  tidyr::pivot_wider(names_from = "is_stable",
                     values_from = "eigenvalues")) |>
    ggplot() +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +  
  geom_hline(yintercept = 0, color = "gray70", linewidth = 1) +
  geom_point(aes(yes, 0), size = 6) +
  geom_point(aes(no, 0), size = 6) +
  ylab("Immaginary") +
  xlab("Real") +
  xlim(-4, 4) +
  ylim(-1, 1) +
  gganimate::transition_reveal(appear) 

anim <- gganimate::animate(animation, 
                fps = 20, 
                duration = 5)

gganimate::anim_save("../../figures/saddle_node.gif",
          width = 5,
          height = 3,
          bg = "white")

## Transcritical ----
transcritical <-  data.table(
  p = seq(-4, 4, by = 0.1)
)[, ":="(
  # find solutions of normal form xp+x^2
  x_1 = 0,
  x_2 = -p,
  # find eigenvalues of Jacobian (first derivative in respect of x in unidimensional case) of normal form at the equilibria
  eigenvalue_1 = p + 2 * 0, 
  eigenvalue_2 = p + 2 * -p)
][order(-p)
  ][, ":="(appear = rleid(as.factor(round(p, 4))))]

(p <- (transcritical |>
    melt(id.vars = "p",
         measure = list(x = c("x_1", "x_2"), eigenvalue = c("eigenvalue_1", "eigenvalue_2")),
         variable.name = "equ_number")
    )[, ":="(is_stable = eigenvalue > 0)] |>
    ggplot() +
    geom_path(aes(p, x, 
                  linetype = is_stable,
                  group = interaction(is_stable, equ_number)
    ), linewidth = 1) +
    xlim(c(-4, 4)) +
    #scale_linetype_manual (values = c("dashed", "solid")) +
    xlab("Parameter p") +
    ylab("Equilibria x*") +
    guides(linetype = "none"))

ggsave("../../figures/transcritical.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

trans_animation <- ggplot(transcritical) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +  
  geom_hline(yintercept = 0, color = "gray70", linewidth = 1) +
  geom_point(aes(eigenvalue_1, 0), size = 6) +
  geom_point(aes(eigenvalue_2, 0), size = 6) +
  ylab("Immaginary") +
  xlab("Real") +
  xlim(-4, 4) +
  ylim(-1, 1) +
  gganimate::transition_reveal(appear) 

trans_anim <- gganimate::animate(trans_animation, 
                           fps = 20, 
                           duration = 5)

gganimate::anim_save("../../figures/transcritical.gif",
                     width = 5,
                     height = 3,
                     bg = "white")

## Pitchfork ----
### Supercritical ----
pitchfork_super <-  data.table(
  p = seq(-4, 4, by = 0.1)
)[, ":="(
  # find solutions of normal form xp-x^3
  x_1 = 0,
  x_2 = -sqrt(p),
  x_3 = sqrt(p),
  # find eigenvalues of Jacobian (first derivative in respect of x in unidimensional case) of normal form at the equilibria
  eigenvalue_1 = p - 3 * 0^2, 
  eigenvalue_2 = p - 3 * (-sqrt(p))^2, 
  eigenvalue_3 = p - 3 * (sqrt(p))^2)
][order(-p)
][, ":="(appear = rleid(as.factor(round(p, 4))))]

(p <- (pitchfork_super |>
         melt(id.vars = "p",
              measure = list(x = c("x_1", "x_2", "x_3"), 
                             eigenvalue = c("eigenvalue_1", "eigenvalue_2", "eigenvalue_3")),
              variable.name = "equ_number")
)[, ":="(is_stable = eigenvalue > 0)] |>
    ggplot() +
    geom_path(aes(p, x, 
                  linetype = is_stable,
                  group = interaction(is_stable, equ_number)
    ), linewidth = 1) +
    xlim(c(-4, 4)) +
    #scale_linetype_manual (values = c("dashed", "solid")) +
    xlab("Parameter p") +
    ylab("Equilibria x*") +
    guides(linetype = "none"))

ggsave("../../figures/supercritical_pitchfork.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

super_pitch_animation <- ggplot(pitchfork_super) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +  
  geom_hline(yintercept = 0, color = "gray70", linewidth = 1) +
  geom_point(aes(eigenvalue_1, 0), size = 6) +
  geom_point(aes(eigenvalue_2, 0), size = 6) +
  ylab("Immaginary") +
  xlab("Real") +
  xlim(-10, 4) +
  ylim(-1, 1) +
  gganimate::transition_reveal(appear) 

super_pitch_anim <- gganimate::animate(super_pitch_animation, 
                                     fps = 20, 
                                     duration = 5)

gganimate::anim_save("../../figures/supercritical_pitchfork.gif",
                     width = 5,
                     height = 3,
                     bg = "white")

### Subcritical ----
pitchfork_sub <-  data.table(
  p = seq(-4, 4, by = 0.1)
)[, ":="(
  # find solutions of normal form xp+x^3
  x_1 = 0,
  x_2 = -sqrt(-p),
  x_3 = sqrt(-p),
  # find eigenvalues of Jacobian (first derivative in respect of x in unidimensional case) of normal form at the equilibria
  eigenvalue_1 = p + 3 * 0^2, 
  eigenvalue_2 = p + 3 * (-sqrt(-p))^2, 
  eigenvalue_3 = p + 3 * (sqrt(-p))^2)
][order(-p)
][, ":="(appear = rleid(as.factor(round(p, 4))))]

(p <- (pitchfork_sub |>
         melt(id.vars = "p",
              measure = list(x = c("x_1", "x_2", "x_3"), 
                             eigenvalue = c("eigenvalue_1", "eigenvalue_2", "eigenvalue_3")),
              variable.name = "equ_number")
)[, ":="(is_stable = eigenvalue > 0)] |>
    ggplot() +
    geom_path(aes(p, x, 
                  linetype = is_stable,
                  group = interaction(is_stable, equ_number)
    ), linewidth = 1) +
    xlim(c(-4, 4)) +
    #scale_linetype_manual (values = c("dashed", "solid")) +
    xlab("Parameter p") +
    ylab("Equilibria x*") +
    guides(linetype = "none"))

ggsave("../../figures/subcritical_pitchfork.png",
       p,
       height = 3,
       width = 5,
       bg = "white")

sub_pitch_animation <- ggplot(pitchfork_sub) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +  
  geom_hline(yintercept = 0, color = "gray70", linewidth = 1) +
  geom_point(aes(eigenvalue_1, 0), size = 6) +
  geom_point(aes(eigenvalue_2, 0), size = 6) +
  ylab("Immaginary") +
  xlab("Real") +
  xlim(-4, 4) +
  ylim(-1, 1) +
  gganimate::transition_reveal(appear) 

sub_pitch_anim <- gganimate::animate(sub_pitch_animation, 
                                 fps = 20, 
                                 duration = 5)

gganimate::anim_save("../../figures/subcritical_pitchfork.gif",
                     width = 5,
                     height = 3,
                     bg = "white")

## Universal unfolding of pitchfork ----
unfolding <- (read.csv("../../outputs/unfolding.csv",
                      header = T) |>
  as.data.table()
  )[, ":="(eigenvalues = p + 2 * a * x - 3 * x^2)
                   ][, ":="(is_stable = eigenvalues < 0)
                     ][order(a, b, p, is_stable, x)
                       ][, diff := p - shift(p), 
                         by = list(a, b)
                         ][, is_change := diff != shift(diff) & is_stable != shift(is_stable)]



unfolding <- (read.csv("../../outputs/unfolding.csv",
                       header = T) |>
                as.data.table()
)[, ":="(eigenvalues = p + 2 * a * x - 3 * x^2)
][, ":="(is_stable = eigenvalues < 0)
# identify branch id by takeing advantage of the fact that the stability switches every time a nullcline is crossed (from large x to small x) and the same x can be equilibria for multiple branches only when the branche changes stability
][order(a, b, x, is_stable)
  ][, id := rleid(is_stable),
    list(a, b)]

(p <- ggplot(unfolding) +
  geom_path(aes(p, x, group = id, 
                linetype = is_stable),
            linewidth = 1) +
  xlab("Parameter p") +
  ylab("Equilibria x*") +
    scale_linetype_manual(values = c("dashed", "solid")) +
  guides(linetype = "none") +
  facet_grid(cols = vars(a),
             rows = vars(b)))

ggsave("../../figures/unfolding.png",
       p,
       height = 9,
       width = 9,
       bg = "white")




# Selective ultra-sensitivity ----
dt <- 0.5
nZ <- 201

RA_alligned <- (read.csv("../../outputs/selective_ultrasensitivity/RA_alligned.csv",
                         header = F) |>
                       as.data.table() |>
                       setnames(old = 1:nZ,
                                new = paste0("z_", 1:nZ))
)[, ":="(time = ((1:.N)-1) * dt)]

RA_alligned_l <- (RA_alligned |>
                         tidyr::pivot_longer(cols = contains("z"),
                                             names_to = "z_id",
                                             values_to = "z") |>
                         as.data.table()
)[, ":="(theta_index = as.numeric(sub("z_", "", z_id)))
][, ":="(theta = seq(-pi, pi, by = 2*pi/nZ)[theta_index])]

RA_disalligned <- (read.csv("../../outputs/selective_ultrasensitivity/RA_disalligned.csv",
                         header = F) |>
                  as.data.table() |>
                  setnames(old = 1:nZ,
                           new = paste0("z_", 1:nZ))
)[, ":="(time = ((1:.N)-1) * dt)]

RA_disalligned_l <- (RA_disalligned |>
                    tidyr::pivot_longer(cols = contains("z"),
                                        names_to = "z_id",
                                        values_to = "z") |>
                    as.data.table()
)[, ":="(theta_index = as.numeric(sub("z_", "", z_id)))
][, ":="(theta = seq(-pi, pi, by = 2*pi/nZ)[theta_index])]

input_alligned <- (read.csv("../../outputs/selective_ultrasensitivity/alligned.csv",
                            header = F) |>
                     setnames(old = 1,
                              new = "input") |>
                     as.data.table()
                   )[, ":="(theta_index = 1:.N)
                              ][, ":="(theta = seq(-pi, pi, by = 2*pi/nZ)[theta_index])]

input_disalligned <- (read.csv("../../outputs/selective_ultrasensitivity/disalligned.csv",
                            header = F) |>
                     setnames(old = 1,
                              new = "input") |>
                     as.data.table()
)[, ":="(theta_index = 1:.N)
][, ":="(theta = seq(-pi, pi, by = 2*pi/nZ)[theta_index])]

(p <- ggplot(input_alligned) +
    geom_point(aes(input, theta)) +
    ylab("Direction \u03B8") +
    scale_x_continuous(breaks = c(-0.01, 0, 0.01),
                       limits = c(-0.01, 0.01)) +    xlab("Input"))

(p1 <- ggplot(RA_alligned_l) +
    geom_tile(aes(time, theta,
                  fill = z)) +
    ylab("Direction \u03B8") +
    xlab("Time") +
    theme(axis.title.y = element_blank()) +
    scale_fill_viridis_c(name = "Activation z",
                         limits = c(-0.1, 2)))

(part1 <- ggpubr::ggarrange(p, p1,
          widths = c(0.4, 1)))

(p <- ggplot(input_disalligned) +
    geom_point(aes(input, theta)) +
    ylab("Direction \u03B8") +
    scale_x_continuous(breaks = c(-0.01, 0, 0.01),
                       limits = c(-0.01, 0.01)) +
    xlab("Input"))

(p1 <- ggplot(RA_disalligned_l) +
    geom_tile(aes(time, theta,
                  fill = z)) +
    ylab("Direction \u03B8") +
    xlab("Time") +
    theme(axis.title.y = element_blank()) +
    scale_fill_viridis_c(name = "Activation z",
                         limits = c(-0.1, 2)))

(part2 <- ggpubr::ggarrange(p, p1,
                            widths = c(0.4, 1)))

(p <- ggpubr::ggarrange(part2, part1,
                            nrow = 2))

ggsave("../../figures/selective_ultrasensitivity.png",
       p,
       width = 8,
       height = 6,
       bg = "white")

