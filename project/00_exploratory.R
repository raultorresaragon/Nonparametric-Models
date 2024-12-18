library(tidyverse)
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# load data and light tidying #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
wd <- read_csv("project/wd_gundeaths.csv") 
lower_names <- str_to_lower(names(wd))
wd <- wd |> `colnames<-`(lower_names)
names(wd)[str_detect(names(wd), "percent.+")] <- 
  str_replace(names(wd)[str_detect(names(wd), "percent.+")], "percent[s]{0,1}","pct")
names(wd) <- str_replace_all(names(wd),"\\.","_")

# ~~~~~~~~~~~~~~~~ #
# outcome variable # 
# ~~~~~~~~~~~~~~~~ #

### Note: our outcome variable is "total firearm deaths per 10,000 per year"

# missingness
sum(is.na(wd$total_firearm_deaths))
names(wd)[names(wd) == "total_hunting_licenses_tags_permits_and_stamps"] <- "hunt_lic"

# even panel?
wd |> group_by(year) |> mutate(N = n()) %>% .$N |> summary()

# turn into ts object aggregated by year
gun_deaths <- wd |> group_by(year) |> 
              summarise(gun_deaths = sum(1e4 * total_firearm_deaths/population, 
                                         na.rm = TRUE)) %>% 
              .$gun_deaths |> ts()

xlabs <- unique(wd$year)
plot(gun_deaths, ylab = "Gun deaths per 10,000", xaxt = "n", 
     main = "Gun deaths\nin Texas per year")
axis(1, 
     at = seq(1, length(gun_deaths), by=3), 
     labels = xlabs[seq(1, length(gun_deaths), by=3)])

# acf / pacf
acf(gun_deaths)
pacf(gun_deaths)

# ~~~~~~~~~~~~~~~~~~~ #
# predictor variables #
# ~~~~~~~~~~~~~~~~~~~ #

miss_tbl <- tibble(state = character(),
                          permit_miss = numeric(),
                          hunt_lic_miss = numeric(),
                          poverty_rate_miss = numeric(),
                          unem_rate_miss = numeric()
)

for(s in unique(wd$state)) {
  subd <- wd[wd$state == s, 
                c("state","year","permit",
                  "hunt_lic","poverty_rate","unemployment_rate")]
  
  row <- as_tibble_row(list("state" = s,
                            "permit_miss" = sum(is.na(subd$permit)),
                            "hunt_lic_miss" = sum(is.na(subd$hunt_lic)),
                            "poverty_rate_miss" = sum(is.na(subd$poverty_rate)),
                            "unem_rate_miss" = sum(is.na(subd$unemployment_rate))
                            ))
  miss_tbl <- bind_rows(miss_tbl, row)
}
miss_tbl

wd <- wd |> filter(state == "Texas")

# permits
sum(wd$permit, na.rm = TRUE)
permits <- wd$permit |> ts()
xlabs <- unique(wd$year)
plot(permits, ylab = "Permits", xaxt = "n", main = "Total number permits issued per year")
axis(1, at = seq(1, length(permits), by=3), labels = xlabs[seq(1, length(permits), by=3)])


# hunting licenses
hunt_lic <- wd$hunt_lic |> ts()
xlabs <- unique(wd$year)
plot(hunt_lic, ylab = "Permits", xaxt = "n", main = "Total number hunting liceses issued per year")
axis(1, at = seq(1, length(hunt_lic), by=3), labels = xlabs[seq(1, length(hunt_lic), by=3)])


# poverty rate
poverty_rate <- wd$poverty_rate |> ts()
xlabs <- unique(wd$year)
plot(poverty_rate, ylab = "Poverty rate", xaxt = "n", main = "Texas poverty rate from 1991 to 2016")
axis(1, at = seq(1, length(poverty_rate), by=3), labels = xlabs[seq(1, length(poverty_rate), by=3)])

# unem rate
unem_rate <- wd$unemployment_rate |> ts()
xlabs <- unique(wd$year)
plot(unem_rate, ylab = "Poverty rate", xaxt = "n", main = "Texas unemployment rate from 1991 to 2016")
axis(1, at = seq(1, length(unem_rate), by=3), labels = xlabs[seq(1, length(unem_rate), by=3)])

# y,

# y, x1, x2, t
xlabs <- unique(wd$year)
plot(wd$year, rethinking::standardize(wd$total_firearm_deaths), 
     type = "l", col = "red4", xlab = "Year", ylab ="standardized values", ylim = c(-2.5, 2.5),
     main = "Standardized gun deaths, permits, and poverty rate \n in Texas from 1991 to 2016")
lines(wd$year, rethinking::standardize(wd$permit), type = "l", col = "blue")
lines(wd$year, rethinking::standardize(wd$poverty_rate), type = "l", col = "orange")
lines(wd$year, rethinking::standardize(wd$unemployment_rate), type = "l", lty = 3, col = "orange")
lines(wd$year, rethinking::standardize(wd$hunt_lic), type = "l", lty = 3, col = "blue")

legend("bottomleft", 
       c("gun deaths","permits","poverty rate","hunt licenses", "unemplyment rate"), 
       lty = c(1,1,1,1,3,3), 
       col = c("red4", "blue", "orange", "blue", "orange"))


