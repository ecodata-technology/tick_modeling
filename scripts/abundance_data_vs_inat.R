library(tidyverse)
library(brms)
library(tidybayes)
library(rinat)
library(USAboundaries)
library(sf)
library(patchwork)
library(ecodatamisc)

theme_set(theme_ecodata())

d <- read_csv("data/cleaned/nymphs_cdd.csv")

d %>% 
  ggplot(aes(x = nymphs_per_hour)) +
  geom_histogram()

m1 <- brm(nymphs_per_hour ~ s(cdd), 
          data = d, 
          chains = 4, cores = 4, 
          control = list(adapt_delta = 0.95))

m1_p <- data.frame(cdd = seq(from = 0, to = max(d$cdd, na.rm = T), by = 1)) %>% 
  add_fitted_draws(m1) %>% 
  median_hdci()

m1_peak_dist <- data.frame(cdd = seq(from = 0, to = max(d$cdd, na.rm = T), by = 1)) %>% 
  add_fitted_draws(m1) %>% 
  ungroup() %>% 
  group_by(.draw) %>% 
  slice_max(.value, n = 1) %>% 
  ungroup() %>% 
  median_hdci()

cdd_med <- m1_peak_dist$cdd
cdd_lo <- m1_peak_dist$cdd.lower
cdd_hi <- m1_peak_dist$cdd.upper

p1 <- m1_p %>% 
  ggplot(aes(x = cdd, y = .value, ymin = .lower, ymax = .upper)) +
  geom_vline(xintercept = cdd_med, color = "#FF0000", linetype = 2) +
  annotate(geom = 'rect', xmin = cdd_lo, xmax = cdd_hi, 
            ymin = -Inf, ymax = Inf,
            fill = "#FF0000", alpha = 0.1) +
  annotate(geom = 'text', label = paste("CDD: ", cdd_med), 
            x = cdd_med + 325, y = -75, color = "#FF0000", 
           family = formals(theme_ecodata)$base_family) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  xlim(0,8000) +
  labs(y = "Estimated Sampling Abundance",
       title = "Estimated <span style='font-family:LexendDecaLight;'>I. scapularis</span> densities in CT",
       subtitle = "Top data from drag sampling at 4 sites. Bottom data are iNat obs. from CT.<br><span style='color:#FF0000;'>Red line</span> is median peak, <span style='color:#ff000066
;'>pink band</span> is 95% credible interval.") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

conn <- us_states(states = "CT")

dt_dat <- get_inat_obs(taxon_name = "Ixodes scapularis", 
                       bounds = conn, maxresults=5000) %>% 
  select(datetime, scientific_name, common_name, observed_on, longitude, latitude) %>%
  add_yearday("datetime", datetime = T) %>%
  filter(year < 2022) %>% 
  add_siteday_cdd(lb = 32)

dt_dens <- make_cdd_density_tbl(dt_dat)

dt_dens_max <- slice_max(dt_dens, dens, n = 1)

p2 <- dt_dens %>% 
  ggplot(aes(x = cdd, y = dens)) +
  geom_line() +
  geom_vline(xintercept = dt_dens_max$cdd, color = "#FF0000", linetype = 2) +
  annotate(geom = 'text', label = paste("CDD: ", round(dt_dens_max$cdd)), 
            x = dt_dens_max$cdd + 325, y = 0.5e-4, color = "#FF0000", 
           family = formals(theme_ecodata)$base_family) +
  ylab("Density of iNat observations") +
  xlab("Cumulative Degree Days") +
  xlim(0,8000)


p1 / p2

ggsave("images/tick_brms_gam_vs_inat_density.jpg", 
       bg = "white", width = 11, height = 8.5, device = grDevices::jpeg)

dt_dat_r <- get_inat_obs(taxon_name = "Ixodes scapularis", 
                       bounds = conn, maxresults = 5000, quality = "research") %>% 
  select(datetime, scientific_name, common_name, observed_on, longitude, latitude) %>%
  make_julian_year("datetime", datetime = T) %>%
  filter(year < 2022) %>% 
  add_siteday_cdd(lb = 32)

p3 <- dt_dat_r %>% 
  ggplot(aes(x = cdd)) +
  geom_density(bw = "bcv") +
  ylab("Density of iNat observations") +
  xlim(0,8000)

p1 / p3

ggsave("images/tick_brms_gam_vs_inat_research_density.jpg", 
       bg = "white")
