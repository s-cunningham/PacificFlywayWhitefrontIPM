
library(tidyverse)


# Summarize the mean values, but also keep the station names
noaaMeans <- function(x) {
  un.day <- unique(x$date)
  means <- tibble()
  for (i in 1:length(un.day)) {
    temp <- subset(x, date==un.day[i])
    temp <- as.data.frame(temp)
    n.stations <- length(unique(temp$id))
    sta.list <- list(temp$id)
    day <- as.character(unique(temp$date))
    values <- temp[,3]
    mean.value <- mean(values, na.rm=TRUE)
    sd.value <- sd(values, na.rm=TRUE)
    new <- tibble(date=day, n.stations=n.stations, stations=sta.list, mean.value, std.dev=sd.value)
    means <- rbind(means, new)
  }
  return(means)
}


