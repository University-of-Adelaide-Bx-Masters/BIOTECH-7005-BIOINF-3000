library(knitr)
library(tidyverse)
library(pander)

include_graphics("images/NewRMarkdown.png")

## ?PlantGrowth

## PlantGrowth %>%
##   group_by(group) %>%
##   summarise(n = n(), Mean = mean(weight))

## PlantGrowth %>%
##   group_by(group) %>%
##   summarise(n = n(), Mean = mean(weight)) %>%
##   pander(caption = "Sample Sizes and average weights for each group")

## ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) +
##   geom_boxplot() +
##   theme_bw() +
##   labs(x = "Treatment Group", y = "Dried Weight (g)")

model_fit <- lm(weight ~ group, data = PlantGrowth)

summary(model_fit)
anova(model_fit)

## model_fit %>% summary() %>% pander()

## model_fit %>% anova() %>% pander()
