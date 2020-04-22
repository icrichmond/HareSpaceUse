# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: April 22, 2020

# This script is for estimating the predation risk of snowshoe hare habitat using structural
# complexity of the environment



# load required packages 
easypackages::packages("tidyr", "lubridate")


# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #

# import structural complexity datasets
canopyclosure <- read.csv("input/CSRawData_CanopyClosure.csv")
head(canopyclosure)
horizcomplex <- read.csv("input/CSRawData_HorizontalComplexity.csv")
head(horizcomplex)
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
  scale_fill_brewer(palette = 2)+
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
canopyhorizcorr <- tibble(cc$Average, hc$Score)
abbacorr <- (cor(abbacorrdata))
ggcorrplot(abbacorr, hc.order = TRUE, lab = TRUE)
ggsave("graphics/StoichModels_2Step/Correlations/ABBAcorr.jpg")

