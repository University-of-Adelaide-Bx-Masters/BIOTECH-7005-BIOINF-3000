library(knitr)
opts_chunk$set(echo = TRUE, include = TRUE, eval = FALSE, 
               warning = FALSE, message = FALSE, fig.align = "center",
               results = 'hide', fig.show = 'hide')

## ?plot
## ?boxplot
## ?hist

library(tidyverse)

## csvFile <- file.path("~", "data", "intro_r", "transport.csv")
## file.exists(csvFile)
## data <- read_csv(csvFile) %>%
##   mutate(height = height/100, BMI = weight / height^2) %>%
##   select(-X1)
## data

ggplot(data, aes(x = weight, y = height))

ggplot(data, aes(x = weight, y = height)) +
  geom_point()

## ?geom_point

ggplot(data, aes(x = weight, y = height, colour = transport)) +
  geom_point()

ggplot(data, aes(x = weight, y = height, colour = transport, shape = gender)) +
  geom_point()

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender))

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender), size = 3)

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender)) +
  geom_smooth(method = "lm")

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender)) +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Weight (kg)", y = "Height (m)", shape = "Gender", colour = "Transport")

ggplot(data, aes(x = weight, y = height)) +
  geom_point(aes(colour = transport, shape = gender)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Weight (kg)", y = "Height (cm)", shape = "Gender", colour = "Transport") +
  facet_wrap(~gender) 

ggplot(data, aes(x = height, fill = gender)) +
  geom_density(alpha = 0.5)

ggplot(data, aes(x = gender, y =height, fill = gender)) +
  geom_boxplot()

ggplot(data, aes(x = gender, y =height, fill = gender)) +
  geom_boxplot() +
  facet_wrap(~transport)

ggplot(data, aes(x = height)) +
  geom_histogram(bins = 20)

ggplot(data, aes(x = height)) +
  geom_histogram(bins = 20, fill = "grey50", colour = "black")

data %>%
  group_by(transport, gender) %>%
  summarise(mn_height = mean(height), sd_height = sd(height)) %>%
  ggplot(aes(x = transport, y = mn_height, fill = transport)) +
  geom_bar(stat = "identity") +
  facet_wrap(~gender) +
  guides(fill = FALSE)

data %>%
  group_by(transport, gender) %>%
  summarise(mn_height = mean(height), sd_height = sd(height)) %>%
  ggplot(aes(x = transport, y = mn_height, fill = transport)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    aes(ymin = mn_height - sd_height, ymax = mn_height + sd_height),
    width = 0.6
  ) +
  facet_wrap(~gender) +
  guides(fill = FALSE)

ggplot(data, aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  theme_bw()

?theme

?element_text

ggplot(data, aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(data, aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(data, aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = c(0.85, 0.15))

ggplot(data, aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = c(0.85, 0.15)) +
  scale_y_continuous(limits = c(1.5, 2)) +
  scale_fill_manual(values = c("grey70", "lightblue"))
