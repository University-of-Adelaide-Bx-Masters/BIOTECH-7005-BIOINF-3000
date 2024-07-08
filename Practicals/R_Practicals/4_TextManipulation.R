library(knitr)
opts_chunk$set(
  echo = TRUE,
  eval = FALSE,
  results = 'hide',
  message = FALSE, 
  warnings =  FALSE
)

## library(tidyverse)
## pcrFile <- file.path("~", "data", "intro_r", "pcr.csv")
## file.exists(pcrFile)
## pcrData <- read_csv(pcrFile)

pcrData$Gene

str_to_title(pcrData$Gene)

pcrData <- mutate(pcrData, Gene = str_to_title(Gene))

pcrData$Gene <- str_to_title(pcrData$Gene)

str_replace(string = "Hi Mum", pattern = "Mum", replacement =  "Dad")

str_replace("Hi Mum", "Mum", "Dad")

str_replace(string = "Hi Mum", pattern = "M..", replacement = "Dad")

str_replace(string = "Hi Mum", pattern = "M.+", replacement = "Dad")
str_replace(string = "Hi Mother", pattern = "M.+", replacement = "Dad")

str_replace(string = "Hi Mother", pattern = "(H.+) (M.+)", replacement = "\\2! \\1!")

str_replace("Hi Mum", "[Mm]", "b")
str_replace_all("Hi Mum", "[Mm]", "b")
str_replace_all("Hi Mum", "[aeiou]", "o")
str_replace_all("Hi Mum", "[a-z]", "o")

str_replace_all("Hi Mum", "(Mum|Dad)", "Parent")
str_replace_all("Hi Dad", "(Mum|Dad)", "Parent")
str_replace_all("Hi Dad", "Hi (Mum|Dad)", "Dear Beloved Parent")

pcrData$Gene %>% str_detect("A")
pcrData$Gene %>% str_count("[A-C]")
colnames(pcrData) %>% str_extract("(Resting|Stim)")
colnames(pcrData) %>% str_remove_all("hr")

colnames(pcrData) %>% str_split_fixed("_", n = 2)
