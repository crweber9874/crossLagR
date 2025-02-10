# # library(dplyr)
# #
# # # simulate_riclpm(waves =3)
# # # This works
# # # Now test on clpm
# #
# # # model_syntax_clpm(waves = 4, model_type = "ri- clpm")
# # #
# # # lavaan::lavaan(model_syntax_clpm(waves = 5, model_type = "ri- clpm"), simulate_riclpm(wave = 5, sample.nobs = 1000)$data)
# # #
# #
# #
# #
# simulate_riclpm(waves = 5,
#                 variance.p = 3,
#                 variance.q = 3)[[2]] %>%
#   reshape_long_sim_cr() %>%
#   group_by(id) %>%
#   mutate(i.average = mean(x),
#          within.x = x-i.average) %>%
#   na.omit() %>%
#   brms::brm(y ~  xlag + ylag,
#             data = .,
#             family = gaussian(),
#             cores = 4,
#             chains = 4,
#             iter = 3000,
#             warmup = 1000,
#             seed = 123,
#             control = list(adapt_delta = 0.99, max_treedepth = 15),
#             prior = c(
#               prior(normal(0, 0.001), class = Intercept)
#             ))
#
#
#
#
# simulate_riclpm(waves = 5,
#                 variance.p = 1,
#                 variance.q = 1)[[2]] %>%
#   reshape_long_sim_cr() %>%
#   group_by(id) %>%
#   mutate(i.average = mean(x),
#          within.x = x-i.average) %>%
#   na.omit() %>%
#   brms::brm(y ~  xlag + ylag + (1 | id),
#             data = .,
#             family = gaussian(),
#             cores = 4,
#             chains = 4,
#             iter = 3000,
#             warmup = 1000,
#             seed = 123,
#             control = list(adapt_delta = 0.99, max_treedepth = 15),
#             prior = c(
#               prior(normal(0, 0.001), class = Intercept)
#             ),
#             file = "my_compiled_model.rds")
#
#
#
# #
# loaded_model <- readRDS("my_compiled_model.rds")
# #
# newdat = simulate_riclpm(waves = 5)[[2]] %>%
#   reshape_long_sim_cr() %>%
#   group_by(id) %>%
#   mutate(i.average = mean(x),
#          within.x = x-i.average) %>%
#   na.omit()
#
# model_fit <- update(loaded_model,
#                     newdata = newdat, # Use the new simulated data here
#                     cores = 10,
#                     chains = 5,
#                     iter = 2000,
#                     warmup = 1000) # estimate brms  model multilegel model
#
# #
# # library(brms)
# # # Compile the Stan model once
# # first_model <- brm(y ~ xlag + ylag + (1 | id),
# #                    data = grouped_data,
# #                    family = gaussian(),
# #                    cores = 4,
# #                    chains = 4,
# #                    iter = 2000,
# #                    warmup = 1000,
# #                    seed = 123,
# #                    file = "my_compiled_model.rds") # Save the compiled model
# #
# # # Load the saved model for subsequent fits
# # loaded_model <- readRDS("my_compiled_model.rds")
# #
# # # Fit the model using the loaded object
# # model_fit <- update(loaded_model,
# #                     newdata = new_data, # Use the new simulated data here
# #                     cores = 4,
# #                     chains = 4,
# #                     iter = 2000,
# #                     warmup = 1000) # estimate brms  model multilegel model
# #
# # brms::brm(y ~ xlag + ylag + (1 | id), data = dat, family = gaussian(), cores = 10, iter = 2000, warmup = 1000, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))
# #
# #
# # # pivot the data frame to long format, where the wave corresponds to the post-fix variable name (e.g., x1, x2, x3), extract eh wave  number from the

#
#
# library(tidyverse)
# library(brms)
#
# N <- 200
#
# dat %>%
#   mutate(xo = as.numeric(NA),
#          p1 = as.numeric(NA),
#          p2 = as.numeric(NA),
#          p3 = as.numeric(NA)) -> df
#
# m1 <-
#   brm(
#     formula =
#       bf(x1 ~ 0 + mi(xo) + mi(p1)) +
#       bf(x2 ~ 0 + mi(xo) + mi(p2)) +
#       bf(x3 ~ 0 + mi(xo) + mi(p3)) +
#       bf(xo | mi() ~ 1) +
#       bf(p1 | mi() ~ 1) +
#       bf(p2 | mi() ~ 1 + mi(p1)) +
#       bf(p3 | mi() ~ 1 + mi(p2)) +
#       set_rescor(rescor = TRUE),
#     family = gaussian(),
#     prior =
#       prior(constant(1), class = "b", resp = "x1") +
#       prior(constant(1), class = "sigma", resp = "x1") +
#       prior(normal(0, 10), class = "b", resp = "x2") +
#       prior(constant(1), class = "sigma", resp = "x2") +
#       prior(normal(0, 10), class = "b", resp = "x3") +
#       prior(constant(1), class = "sigma", resp = "x3") +
#       prior(normal(0, 10), class = "Intercept", resp = "xo") +
#       prior(cauchy(0, 1), class = "sigma", resp = "xo") +
#       prior(normal(0, 10), class = "b", resp = "p1") +
#       prior(cauchy(0, 1), class = "sigma", resp = "p1"),
#     data = df,
#     backend = "cmdstanr",
#     cores = 10,
#     chains = 1,
#     threads = threading(2),
#     refresh = 5
#   )
