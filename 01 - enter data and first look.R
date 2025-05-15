# Thu Apr 10 16:03:49 2025 ------------------------------

# bring in data supplied by Pierre-Paul Bitton on number if times of marked
# puffins resighted on Great Island, 2021-2023. 

# Using Riecke et al. 2022 approach to incorporate number of times seen as
# an overdispersed Poisson instead of a binary response

# this method brings in the state of available for resight, similar to the 
# robust design approach, but without the need to sub-set resight intervals
# within year

library(here)

# fetch the data - assumes only one .csv file in data folder
myfile <- dir(here('data'), 'csv')

# call it df
df <- read.csv(here('data',myfile))

# create year variable
df$year <- substr(df$date, 1,4)

# what do we have?
table(df$year)
table(df$ind)
length(table(df$ind))
table(df$ind, df$year)
table(df$pop)
table(df$Activity, df$year)

# shows 75 banded in 2022, but only 74, 1255-03430 is shown as banded in 
# 2022 and 2021
# the activity column gets ignored in the ch-array creation, so no need to fix


sumdf <- as.data.frame(table(df$ind, df$year))

names(sumdf) <- c('ind', 'year', 'count')

# re-arrange data frame to get it ready for capture-history output
mytable <- table(df$ind, df$year)
dimnames(mytable)
mytab <- data.frame(ind = dimnames(mytable)[[1]],
                    byear = '2021')

mytab$byear[51:124] <- 2022

mine <- cbind(mytab, matrix(table(df$ind, df$year)[,1:3], nrow = 124))


# some quick and dirty histograms, only the years post-banding
hist(mine$`2`[mine$byear == 2021], breaks = max(mine$`2`[mine$byear == 2021]))
hist(mine$`3`[mine$byear == 2021], breaks = max(mine$`3`[mine$byear == 2021]))
hist(mine$`3`[mine$byear == 2022], breaks = max(mine$`3`[mine$byear == 2022]))

mean(mine$`2`[mine$byear == 2021])
mean(mine$`3`[mine$byear == 2021])
mean(mine$`3`[mine$byear == 2022])

var(mine$`2`[mine$byear == 2021])
var(mine$`3`[mine$byear == 2021])
var(mine$`3`[mine$byear == 2022])

# looks like more observations and more variance in 2023 compared to 2022
# don't see a cohort effect thought, about the same distributions in 2023
# for 2021 and 2022 cohort

ch <- mine[,2:5]

# save the ch-array with counts as ch

save(ch, file = here('data', 'chout.Rdata'))

