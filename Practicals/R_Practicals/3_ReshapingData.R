library(knitr)
opts_chunk$set(echo = TRUE,
               eval = FALSE,
               results = "hide",
               message = FALSE, warning =  FALSE)

library(tidyverse)
csvFile <- file.path("~", "data", "intro_r", "transport.csv")
file.exists(csvFile)
data <- read_csv(csvFile)
data

dim(data)
nrow(data)
ncol(data)

data[1, 2]

data[1:5, "name"]

select(data, name, transport)

select(data, name, everything())
select(data, ends_with("ght"))
select(data, contains("a"))

select(data, -name)

data <- select(data, -X1)

mutate(data, height_m = height/100)

mutate(data, height_m = height/100, BMI = weight / height_m^2)

data <- mutate(data, height_m = height/100, BMI = weight / height_m^2)

data <- rename(data, height_cm = height)
data

filter(data, gender == "male")

filter(data, gender != "male")

filter(data, gender == "male", height_cm > 175)
filter(data, transport == "car", gender == "female")

arrange(data, weight)
arrange(data, desc(weight))

arrange(data, transport, height_cm)

data %>% arrange(transport, height_cm)

data %>% filter(transport == "bike") %>% arrange(weight)

data %>% # Take our original dataset
  filter(transport == "bike") %>%  # Find the cyclists
  arrange(weight) # Arrange by weight

data %>% summarise(mean(weight), mean(height_cm))

data %>%
  filter(gender == "female",
         transport == "bike") %>%
  summarise(max_BMI = max(BMI), 
            mn_height = mean(height_cm))

data %>%
  group_by(gender, transport) %>%
  summarise(mn_weight = mean(weight), 
            mn_height =mean(height_cm),
            mn_BMI = mean(BMI),
            n = n())

## library(tidyverse)
## pcrFile <- file.path("~", "data", "intro_r", "pcr.csv")
## file.exists(pcrFile)
## pcrData <- read_csv(pcrFile)

pcrData %>% pivot_longer(-Gene, names_to = "variable", values_to = "value")

pcrLong <- pcrData %>% pivot_longer(-Gene, names_to = "Treatment", values_to = "Ct")
head(pcrLong)

pcrLong %>%
  separate(Treatment, into = c("Treatment", "Timepoint"), sep = "_")

pcrLong %>% pivot_wider(names_from = "Treatment", values_from = "Ct")
