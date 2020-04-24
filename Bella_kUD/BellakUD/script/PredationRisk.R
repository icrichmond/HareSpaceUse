# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: April 22, 2020

# This script is for estimating the predation risk of snowshoe hare habitat using structural
# complexity of the environment

# source code for the kUD file in this code is KernelEstimation.R 

# load required packages 
easypackages::packages("tidyr", "lubridate", "dplyr", "ggplot2", "ggcorrplot", "sf", "raster")


# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #

# import structural complexity datasets
canopyclosure <- read.csv("input/CSRawData_CanopyClosure.csv")
head(canopyclosure)
horizcomplex <- read.csv("input/CSRawData_HorizontalComplexity.csv")
head(horizcomplex)
# load lowland blueberry C:N stoich raster (as per Juliana's results)
vaancn <- raster("input/VAAN_CN.tif")
# clip raster to study area 
e <- extent(860000, 863000, 5383000, 5386000)
vaancnclip <- crop(vaancn, e)
image(vaancnclip)
# load complexity sampling locations shapefile
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- st_transform(bl_cs_pts, crs = st_crs(vaancn))
# load vUD raster layer 
vUDBrick <- raster("output/vUDRaster.tif")

# convert datasets into tibbles, rename plots, code dates to make manipulation easier 
cc <- as_tibble(canopyclosure, .name_repair = "universal")
cc <- cc %>% rename(Plot = ï..Plot) %>%
  mutate(Date = lubridate::ymd(Date))

hc <- as_tibble(horizcomplex, .name_repair = "universal")
hc <- hc %>% rename(Plot = ï..Plot) %>%
  mutate_all(factor) %>%
  mutate(Date = lubridate::ymd(Date))
hc <- drop_na(hc)


# horizontal complexity data was measured using a Nudds board (Nudds, 1977)
# to choose what distance we use to measure horizontal complexity, we need to determine
# which distance has the greatest variation in scores 
hcsum <- hc %>%
  group_by(Distance, Score) %>%
  tally()

ggplot(hcsum, aes(x=Distance,y=n, fill=Score))+
  geom_col()+
  scale_fill_brewer(palette = 1)+
  theme(panel.background = element_rect(fill = "darkgrey"))+
  labs(fill="Score")
# inspection of the graph shows that a distance of 10 metres experiences the most
# variation in horizontal complexity scores
# therefore, use values from 10 m only going forward 
hc <- hc %>% filter(Distance == 10)
# terrestrial predators are all under 1 m tall (coyote & lynx), subset these heights
# from the dataset and average for each plot 
hc <- hc %>% filter(Height %in% c(0.5,1))
# calculate an average horizontal complexity score for each plot
# to do this we have to change the Score variable to an integer
hc <- hc %>% mutate(score = as.numeric(Score))

hcmean <- hc %>%
  group_by(Plot) %>%
  summarise(meanhc = mean(score))

cchc <- inner_join(cc, hcmean, by = "Plot")

# test for correlation between these two variables 
canopyhorizcorr <- tibble(cchc$CoverValue, cchc$meanhc)
cchccorr <- (cor(canopyhorizcorr))
ggcorrplot(cchccorr, hc.order = TRUE, lab = TRUE)
# correlation is not highly significant (-0.32)

# make canopy cover and horizontal complexity data spatially explicit by combining
# coordinate data and cchc data
cchc.spatial <- inner_join(cchc, bl_cs_pts, by = "Plot")

# convert VAAN C:N raster to dataframe to plot in ggplot
vaancnclip.df <- as.data.frame(vaancnclip, xy=TRUE)
# plot the stoich layer with the 
ggplot(cchc.spatial, aes(x = POINT_X, y = POINT_Y))+
  geom_raster(aes(x=x, y=y, fill = VAAN_CN), data = vaancnclip.df)+
  geom_point()

# extract the stoich and kUD values at each sampling point
# create raster brick of kUD values and stoich values for easier extraction 
stoichkUD <- brick(list(vaancnclip, vUDBrick))
csValues <- extract(stoichkUD, cchc.spatial)
