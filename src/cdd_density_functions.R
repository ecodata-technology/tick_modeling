library(tidyverse)
library(daymetr)

#TODO: allow rounding of lat/lon for Daymet pulls
#TODO: parallelize Daymet calls
#TODO: can we try scaling temps to Kelvins and see if that affects CDD calculations?
# TODO: come up with better matching than findInterval. actually i'm maybe convinced that this is better. an observation of CDD 95 could not happen at CDD 100, and if CDD 100 is reached on day 24, an observation of CDD 95 would happen the day before. on the other hand, CDD is a tally that ENDS the day. so if the CDD for 100 is on day 24, then it probably passed THROUGH 95 on the way to 100. maybe mess with the bounding of findInterval a bit?

# convenience version of the download_daymet function
download_daymet_v <- function (site = "Daymet", lat = 36.0133, lon = -84.2625, start = 2000, end = as.numeric(format(Sys.time(), "%Y")) - 2, vars = "tmax,tmin,dayl,prcp,srad,swe,vp", path = tempdir(), internal = TRUE, silent = FALSE, force = FALSE, simplify = FALSE) 
{
  if (!silent & !internal & identical(path, tempdir())) {
    message("NOTE: by default data is stored in tempdir() ...")
  }
  url <- daymetr:::server()
  if (!force) {
    max_year <- as.numeric(format(Sys.time(), "%Y")) - 1
  }
  else {
    max_year <- as.numeric(format(Sys.time(), "%Y"))
  }
  if (start < 1980) {
    stop("Start year preceeds valid data range!")
  }
  if (end > max_year) {
    stop("End year exceeds valid data range!")
  }
  year_range <- paste(seq(start, end, by = 1), collapse = ",")
  query <- list(lat = lat, lon = lon, vars = vars, 
                year = year_range)
  daymet_file <- file.path(normalizePath(path), sprintf("%s_%s_%s.csv", 
                                                        site, start, end))
  daymet_tmp_file <- file.path(normalizePath(tempdir()), sprintf("%s_%s_%s.csv", 
                                                                 site, start, end))
  if (!silent) {
    message(paste("Downloading DAYMET data for: ", site, 
                  " at ", lat, "/", lon, " latitude/longitude !\n", 
                  sep = ""))
  }
  error <- httr::GET(url = url, query = query, httr::write_disk(path = daymet_tmp_file, 
                                                                overwrite = TRUE))
  if (httr::status_code(error) == 400) {
    file.remove(daymet_tmp_file)
    stop("Your requested data is outside DAYMET spatial coverage.\n\n            Check the requested coordinates.")
  }
  if (httr::status_code(error) > 400) {
    file.remove(daymet_tmp_file)
    stop("The server is unreachable, check your connection.")
  }
  if (!silent) {
    message("Done !\n")
  }
  if (internal) {
    tmp_struct <- read_daymet(daymet_tmp_file, site = site, 
                              simplify = simplify)
    return(tmp_struct)
  }
  else {
    if (!identical(daymet_tmp_file, daymet_file)) {
      file.copy(daymet_tmp_file, daymet_file, overwrite = TRUE, 
                copy.mode = FALSE)
      invisible(file.remove(daymet_tmp_file))
    }
    else {
      message("Output path == tempdir(), file not copied or removed!")
    }
  }
}

make_julian_year <- function(data, date_col, datetime = F, keep_date_col = T){
  
  if(datetime){
    data[[date_col]]<- as_datetime(data[[date_col]])
  } else {
    data[[date_col]]<- as_date(data[[date_col]])
  }
  
  data$julian_day <- yday(data[[date_col]])
  data$year <- year(data[[date_col]])
  
  if(!keep_date_col){
    data <- select(data, -all_of(date_col))
  }
  
  return(data)
}

# get dataframe of cdd values for unique site-year combinations in the given data
get_siteyear_cdd <- function(data, ub, lb){
  ad <- data %>% 
    distinct(latitude, longitude, year)
  
  ad <- ad %>% 
    mutate(data = pmap(., .f = possibly(function(latitude, longitude, year){
      
      if(year == year(today())){
        download_daymet_v(site = "DummySite", lat = latitude, lon = longitude,
                          start = year-1, end = year, 
                          internal = T, force = F, vars = "tmax,tmin") %>% 
          .$data %>% 
          filter(year == year) %>% 
          select(yday, starts_with("t"))
        
      }
      
      download_daymet_v(site = "DummySite", lat = latitude, lon = longitude,
                        start = year, end = year, 
                        internal = T, force = F,
                        vars = "tmax,tmin") %>% 
        .$data %>% 
        select(yday, starts_with("t"))
      }, otherwise = data.frame(NULL))))
  
  ad %>% 
    unnest(data) %>% 
    mutate(tavg = ((tmax..deg.c. + tmin..deg.c.)/2) * (9/5) + 32) %>% 
    select(-tmax..deg.c., -tmin..deg.c., julian_day = yday) %>% 
    mutate(tavg = if_else(tavg < lb | tavg > ub, lb, tavg),
           dd = tavg - lb) %>% 
    group_by(latitude, longitude, year) %>% 
    mutate(cdd = cumsum(dd)) %>% 
    select(-tavg, -dd) %>% 
    ungroup()
}

# add cdd to each row of the given dataframe
# data argument needs to have columns:
# SiteID
# latitude
# longitude
# year
# julian_day

add_siteyear_cdd <- function(data, ub = 120, lb = 50){
  
  left_join(data,
            get_siteyear_cdd(data = data, ub = ub, lb = lb),
            by = c("latitude", "longitude", 
                   "year", "julian_day"))
}


# density plot functions --------------------------------------------------

get_site_cdd_curves <- function(sites, ub, lb){
  s <- sites %>% 
    get_siteyear_cdd(ub = ub, lb = lb) %>% 
    left_join(sites %>% select(siteid, latitude, longitude),
              by = c("latitude", "longitude")) %>% 
    mutate(date = parse_date_time(
      paste(year, as.character(julian_day)), orders = "yj") %>% 
        as_date()
    ) %>% 
    select(siteid, date, cdd)
  
  ss <- s %>% 
    filter(cdd > min(s$cdd, na.rm = T), cdd < max(s$cdd, na.rm = T))
  
  s0 <- s %>% 
    filter(cdd == min(s$cdd, na.rm = T)) %>% 
    slice_max(date, n = 1)
  
  sm <- s %>% 
    filter(cdd == max(s$cdd, na.rm = T)) %>% 
    slice_min(date, n = 1)
  
  bind_rows(s0, ss, sm)
}

make_breaks <- function(site_curves, nbreaks){
  
  rows <- findInterval(
    seq(from = min(site_curves$cdd), 
        to = max(site_curves$cdd), 
        length.out = nbreaks),
    site_curves$cdd)
  
  site_curves %>%
    filter(row_number() %in% rows) %>%
    mutate(lab = format(date, "%b %d"))
  
}

make_cdd_density_tbl <- function(data){
  
  dens <- density(data[["cdd"]], na.rm = T)
  tibble(cdd = dens[["x"]], dens = dens[["y"]]) %>%
    filter(cdd > 0, cdd < max(data[["cdd"]], na.rm = T))
}

match_cdd_obs_site <- function(obs_data, site_curve){
  obs_data %>% 
    mutate(site_cdd_index = findInterval(cdd, site_curve[["cdd"]]),
           site_cdd_match = site_curve[["cdd"]][site_cdd_index],
           site_date_match = site_curve[["date"]][site_cdd_index] %>% 
             as_date()) %>% 
    select(-site_cdd_index)
}


