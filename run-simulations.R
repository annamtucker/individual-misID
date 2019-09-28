### simulations for TWS talk 2019
## comparing encounter type and frequency

# load packages and functions ----
# be sure the working directory is set to the same place as this script

source("load-packages.R")
source("zombie_function.R")
source("ghost_function.R")
source("helper-functions.R")


# simulation inputs ----
# all inputs except reps can be one or multiple values

# misread error rates
error.rates = 0.05

# study lengths 
yrs = 20

# true survival probability 
phi = 0.6

# true detection probability 
p = c(0.2, 0.5, 0.8)

# number of new individuals marked/encountered each year 
n.ind = 100

# number of replicates of each scenario
reps = 1000

# create unique scenarios
scenarios = as_tibble(expand.grid(years = yrs,
                                  n = n.ind,
                                  phi = phi,
                                  p = p,
                                  error = error.rates))
scenarios$scenario = c(1:nrow(scenarios))

# replicate scenarios
sims = as_tibble(expand.grid(scenario = c(1:max(scenarios$scenario)),
                             replicate = c(1:reps)))
sims <- sims %>%
  full_join(scenarios)


# run simulation and extract results ----
# these functions will simulate data, apply errors, 
# and fit CJS model phi(t)p(.) to both filtered and unfiltered datasets

# data filtering = removing first obs of each new individual (ghost) or 
# removing first obs of every ind in every occastion (zombie)

sim.results <- sims %>%
  mutate(zombie = pmap(list(years, n, error, phi, p), sim_zombie),
         ghost = pmap(list(years, n, error, phi, p), sim_ghost))

# extract results 
results <- sim.results %>% 
  unnest(zombie, ghost) %>% 
  gather(error.type, model.results, 8:9) %>% 
  mutate(type = rep(c("unfiltered", "filtered"), nrow(sims)*2)) 



# plot outputs ----

phis <-  sims %>% 
  mutate(ests = map(model.results, extract.reals)) %>% 
  select(-model.results) %>% 
  unnest(ests) %>% 
  filter(parm == "phi")


# compare detection probs

# choose to plot either filtered or unfiltered data estimates
plot.type = "filtered"

phis %>% 
  filter(type == plot.type & year < (yrs-1)) %>%
  group_by(error.type, year, p) %>% 
  summarize(lcl = quantile(estimate, 0.025),
            med = median(estimate),
            ucl = quantile(estimate, 0.975)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = med)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.4, fill = "gray60", col = NA) +
  geom_abline(lwd = 1.5, intercept = phi, slope = 0, lty = 2) +
  geom_line(lwd = 1.5, col = "gray40") +
  geom_smooth(method = "lm", se = F, lwd = 2, alpha = 0.2, col = "black") +
  facet_grid(p~error.type, scales = "free", labeller = labbeller(p = label_both)) +
  ylim(0.2, 1) +
  xlab("Year") +
  ylab("Survival estimate") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 24, margin = margin(4,4,4,4)),
        strip.text.y = element_text(size = 24, margin = margin(4,4,4,4), angle = 0),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)




# rmse 

rmse <- phis %>% 
  group_by(error.type, type, p) %>% 
  summarize(rmse = sqrt(sum((estimate - phi)^2)/(length(estimate)-1))) %>% 
  ungroup()

rmse %>% 
  ggplot(aes(x = p, y = rmse, col = as.character(type))) +
  geom_line(lwd = 1.5, lty = 2) +
  geom_point(size = 6) +
  facet_wrap(~error.type, scales = "free") +
  ylim(0, 0.5) +
  scale_color_manual(labels = c("Without filtering",
                                "With filtering"),
                     name = "",
                     values = c("palegreen3",
                                "dodgerblue3") ) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 24, margin = margin(4,4,4,4)),
        strip.text.y = element_text(size = 24, margin = margin(4,4,4,4), angle = 0),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        legend.position = "top",
        legend.text = element_text(size = 24)) +
  xlab("Detection probability (p)") +
  ylab("Root mean squared error (RMSE)")

