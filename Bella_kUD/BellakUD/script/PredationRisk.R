# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: April 22, 2020

# This script is for estimating the predation risk of snowshoe hare habitat using structural
# complexity of the environment

# source code for the kUD file in this code is KernelEstimation.R 

# load required packages 
easypackages::packages("maptools", "sf", "matrixStats", "tidyr", "lubridate", "dplyr", "ggplot2", "ggcorrplot", "sf", "raster", "lme4", "RCurl", "tibble", "rgdal")

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
vUDBrick <- brick("output/vUDRaster.grd")

# convert datasets into tibbles, rename plots, code dates to make manipulation easier 
cc <- as_tibble(canopyclosure, .name_repair = "universal")
cc <- cc %>% dplyr::rename(Plot = ï..Plot) %>%
  mutate(Date = lubridate::ymd(Date))

hc <- as_tibble(horizcomplex, .name_repair = "universal")
hc <- hc %>% dplyr::rename(Plot = ï..Plot) %>%
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
  dplyr::summarise(meanhc = mean(score))

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
# plot the stoich layer with the sampling points
ggplot(cchc.spatial, aes(x = POINT_X_x, y = POINT_Y_y))+
  geom_raster(aes(x=x, y=y, fill = VAAN_CN), data = vaancnclip.df)+
  geom_point()
# convert cchc.spatial into a SpatialPointsDataFrame
pts <- cbind(cchc.spatial$POINT_X_x, cchc.spatial$POINT_Y_y)
cchc.spatial <- SpatialPointsDataFrame(pts, cchc.spatial, proj4string = CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
# extract the stoich and kUD values at each sampling point
# create raster brick of kUD values and stoich values for easier extraction 
stoichkUD <- brick(list(vaancnclip, vUDBrick))
csValues <- raster::extract(stoichkUD, cchc.spatial, df = TRUE)
# want to combine cchc and csValues dataframes so that there are all values for 
# each sampling plot 
cchc <- as_tibble(cchc) %>%
  mutate(Plot = as.character(Plot))
csValues <- as_tibble(csValues) %>%
  mutate(ID = as.character(ID))
# use bind and not join because the plots are in the same order 
cscchc <- bind_cols(csValues, cchc)
# get median value of all collars for each sampling plot 
cscchc <- as_tibble(cscchc) %>%
  rowwise() %>% 
  mutate(kUDmedian =  median(c(!!! rlang::syms(grep('X', names(.), value=TRUE)))))

# --------------------------------------- #
#   Generalized Linear Models with kUD    #
# --------------------------------------- #
# now we have plots with average horizontal complexity values, 
# canopy closure values, kernel utilization values for all 33 
# hares, median kUD values for every plot and stoich values for 
# lowland blueberry C:N

# first going to test if there is a relationship between median values and 
# predation risk/food quality
medianmodel <- glm(kUDmedian ~ CoverValue + meanhc + VAAN_CN, family = Gamma, data = cscchc)
plot(medianmodel)
# residuals are extremely non-normal due to skew in the data


# want to test to see if intensity of use is explained by 
# predation risk and/or food quality
# also include individual and plot as random effects 
# using lme4 as these are mixed effect models
# stack the data so it is ready for analysis
cscchc_stack <- pivot_longer(cscchc, cols = starts_with("X"), names_to = "CollarID", names_prefix = "X", values_to = "kUD")
# fix dataset so appropriate columns are factors
cscchc_stack <- cscchc_stack %>% mutate(Plot = as.factor(Plot)) %>%
  mutate(CollarID = as.factor(CollarID))
head(cscchc_stack)
# standardize the explanatory variables
cscchc_stack <- cscchc_stack %>%
  add_column(kUD_s = scale(cscchc_stack$kUD, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CN_s = scale(cscchc_stack$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(CoverValue_s = scale(cscchc_stack$CoverValue, center = TRUE, scale = TRUE)) %>%
  add_column(meanhc_s = scale(cscchc_stack$meanhc, center = TRUE, scale = TRUE))
# in order to make the model work, nAGQ is set to zero. This 
# indicates less precision with respect to the effects of the 
# random variables
indmod <- glmer(kUD ~ CoverValue_s + meanhc_s + VAAN_CN_s + (1|CollarID) + (1|Plot), nAGQ=0, data=cscchc_stack, family = inverse.gaussian)
plot(indmod)
qqnorm(residuals(indmod))
qqline(residuals(indmod))
# model is still extremely non-normal due to intensely skewed data 

# --------------------------------------- #
#     Generalized Linear Models with      #
#             binomial kUD                #
# --------------------------------------- #
# Going to convert continuous response variable (kUD and median kUD) to a binomial 
# to address this skew. If kUD > 0.90, considered "high" use (1), and if kUD < 0.90m 
# considered "low" use (0)
bicscchc <- as_tibble(cscchc) %>% 
  mutate_each(funs(ifelse(. >= 90,1,0)), contains(c("X", "median")))

bimedianmodel <- glm(kUDmedian ~ CoverValue + meanhc + VAAN_CN, family = binomial, data = bicscchc)
qqnorm(residuals(bimedianmodel))
qqline(residuals(bimedianmodel))
# still extremely non-normal

bicscchc_stack <- pivot_longer(bicscchc, cols = starts_with("X"), names_to = "CollarID", names_prefix = "X", values_to = "kUD")

# fix dataset so appropriate columns are factors
bicscchc_stack <- bicscchc_stack %>% mutate(Plot = as.factor(Plot)) %>%
  mutate(CollarID = as.factor(CollarID))
head(bicscchc_stack)
# standardize the explanatory variables
bicscchc_stack <- bicscchc_stack %>%
  add_column(VAAN_CN_s = scale(bicscchc_stack$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(CoverValue_s = scale(bicscchc_stack$CoverValue, center = TRUE, scale = TRUE)) %>%
  add_column(meanhc_s = scale(bicscchc_stack$meanhc, center = TRUE, scale = TRUE))
# in order to make the model work, nAGQ is set to zero. This 
# indicates less precision with respect to the effects of the 
# random variables
biindmod <- glmer(kUD ~ CoverValue_s + meanhc_s + VAAN_CN_s + (1|CollarID) + (1|Plot), nAGQ=0, data=bicscchc_stack, family = binomial)
plot(biindmod)
qqnorm(residuals(biindmod))
qqline(residuals(biindmod))
head(cscchc_stack)
# still extremely non-normal

# --------------------------------------- #
#     Generalized Linear Models with      #
#         kUD home range overlap          #
# --------------------------------------- #
# because the kUD data is so skewed - space use is generally high, throughout the grid
# I am going to use the 50% home ranges or "core usage areas" of the individuals 
# and do an analysis with degree of overlap as the response variable
# load the 50% kernel home range
hares.kUDhr.50 <- readOGR("output", "hares.kudhr.50")
hares.kUDhr.50 <- spTransform(hares.kUDhr.50, CRSobj = "+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs")
# map the hare home ranges with the complexity sampling points 
hares.kUDhr.50.points <- fortify(hares.kUDhr.50, region = "id")
hares.kUDhr.50.df <- plyr::join(hares.kUDhr.50.points, hares.kUDhr.50@data, by = "id")

ggplot(hares.kUDhr.50.df, aes(long, lat))+
  geom_polygon(aes(group = id, fill = id, alpha = 0.25))+
  geom_point(aes(x = POINT_X_x, y = POINT_Y_y),bl_cs_pts)+
  coord_equal()+
  scale_fill_discrete()
# extract the home ranges that overlap each complexity sampling point
overlapValues <- raster::extract(hares.kUDhr.50, cchc.spatial,fun=sum,df = TRUE)
# count the number of overlapping home ranges  
overlapcount <- overlapValues %>% as_tibble(overlapValues) %>%
  group_by(point.ID) %>%
  summarise(n())
overlapcount <- rename(overlapcount, overlap = "n()")
overlapValues <- inner_join(overlapValues, overlapcount, by = "point.ID") 
# make sure that the overlap count did not count NA rows as 1 
overlapValues$overlap <- ifelse(is.na(overlapValues$poly.ID), 0.001, overlapValues$overlap)
# removing many redundant rows and keeping one row per complexity sampling point
overlapcount <- overlapValues %>% distinct(overlapValues,point.ID,.keep_all = TRUE)
# join with explanatory variables df so that we can run analysis 
overlapdata <- bind_cols(overlapcount, cscchc)
# now run the analysis 
overlapmodel <- glm(overlap ~ VAAN_CN + CoverValue + meanhc, data = overlapdata)
summary(overlapmodel)
plot(overlapmodel)
