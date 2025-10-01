###### MAP OF SITES WITH PATTERNS ########
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)

cases <- read.csv(file = '~/2025-interactions/all_cases.csv', header = TRUE)
meta_mags <- read.csv("~/MAGs_BGCs_interactions/metadata.csv")

example <- recreate_table("meta_mOTU_v25_12843", "gcf_871", mags_by_sites, bgcs_by_sites)

#mOTU12843 <- meta_bgcs %>%
#  filter(gcf == "gcf_871")
#mOTU12843 <- meta_mags %>%
#  filter(mOTUs_Species_Cluster == "meta_mOTU_v25_12843")


#sitios con cordenadas
sites <- data.frame(
  station = meta_mags$station,
  latitude = meta_mags$latitude,
  longitude = meta_mags$longitude
)
sites <- sites %>%
  group_by(station) %>%
  summarise(
    longitude = first(longitude),
    latitude = first(latitude)
  )
#unir sitios a la tabla
example <- example %>%
  left_join(sites %>% select(latitude, longitude, station), by = "station")

example <- example %>%
  mutate(tipo = case_when(
    meta_mOTU_v25_12843 > 0 & gcf_871 == 0 ~ "meta_mOTU",
    meta_mOTU_v25_12843 == 0 & gcf_871 > 0 ~ "gcf",
    meta_mOTU_v25_12843 > 0 & gcf_871 > 0  ~ "ambos",
    TRUE                                   ~ "ninguno"
  ))
cols <- c(
  "meta_mOTU" = "#1f77b4",  # azul
  "gcf"       = "#d62728",  # rojo
  "ambos"     = "purple",  # morado
  "ninguno"   = "grey60"    # gris
)

# convertir a objeto espacial para terra
pts <- vect(example, geom = c("longitude", "latitude"), crs = "EPSG:4326")
mundo <- ne_countries(scale = "medium", returnclass = "sf")
oceanos <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = oceanos, fill = "#b4d3e2", color = NA) +
  geom_sf(data = mundo, fill = "#5b6684", color = "#5b6684") +
  geom_point(data = example, aes(x = longitude, y = latitude, color = tipo),
             alpha = 0.9, size = 2) +
  scale_color_manual(values = cols, name = "Tipo de interacción") +
  theme_minimal() +
  theme(legend.position = "bottom")
map









