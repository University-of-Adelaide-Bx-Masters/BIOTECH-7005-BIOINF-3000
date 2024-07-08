library(knitr)
opts_chunk$set(
  echo = TRUE,
  eval = TRUE,
  message = FALSE, 
  warnings =  FALSE
)

library(tidyverse)

knitr::include_graphics("images/importGUI_preview.png")

knitr::include_graphics("images/importGUI_code.png")

## ?read_csv

## transport <- read_csv("~/data/intro_r/transport.csv")

## transport <- read_csv(file = "~/data/intro_r/transport.csv")

## transport <- read_csv("path/to/file", skip = 3)

opts_chunk$set(eval = FALSE)

## View(transport)
## transport
## head(transport)

no_header <- read_csv("~/data/intro_r/no_header.csv")
no_header

no_header <- read_csv("~/data/intro_r/no_header.csv", col_names = FALSE)

no_header <- read_csv("~/data/intro_r/no_header.csv", col_names = FALSE, col_types = "-ccnnc")

no_header <- read_csv("~/data/intro_r/no_header.csv", col_names = FALSE, col_types = "-cnnnc")

comments <- read_csv("~/data/intro_r/comments.csv")

comments <- read_csv("~/data/intro_r/comments.csv", comment = "#")

bad_colnames <- read_csv("~/data/intro_r/bad_colnames.csv")

bad_colnames <- read_csv("~/data/intro_r/bad_colnames.csv",  skip =  1, col_names = FALSE)
colnames(bad_colnames)
colnames(bad_colnames) <- c("rowname", "gender", "name", "weight", "height", "transport")

c()

missing_data <- read_csv("~/data/intro_r/missing_data.csv")

missing_data <- read_csv("~/data/intro_r/missing_data.csv", na = "-", col_types = "-ccnnc")
