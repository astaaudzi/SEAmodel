list.of.packages <- 
  c("tidyverse",
    "cowplot", 
    "googleway", 
    "ggrepel", 
    "ggspatial", 
    "libwgeom", 
    "sf", 
    "rnaturalearth", 
    "rnaturalearthdata",
    "rgeos"
  )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = T)
rm(list.of.packages, new.packages)

theme_set(theme_bw())

load(file = "datasets/TasmModel_sites.RData")

world <- ne_countries(scale = "medium", returnclass = "sf")
#class(world)

#center_tas <- c(-41.4545, 145.9707)

#Tassie 
ggplot(data = world) +
  geom_sf(col="black", fill="grey") +
  coord_sf(xlim = c(143, 151), ylim = c(-39.5, -44.5), expand = FALSE) +
  geom_point(x= study_sites$long, y = study_sites$latt, size = study_sites$surveys/80, data = study_sites, col = 'red') 
#  theme_void()

#inset 
ggplot(data = world) +
  geom_sf(col="black", fill="grey") +
  coord_sf(xlim = c(110, 155), ylim = c(-10, -45), expand = FALSE) +
    theme_void()

