library(knitr)
opts_chunk$set(echo = TRUE,
               eval = FALSE,
               message = FALSE, warnings=  FALSE)

## 1 + 1
## 2 * 2
## 1 / 4
## 3 ^ 2

knitr::include_graphics("images/RStudio.png")

x <- 5

x
print(x)

knitr::include_graphics("images/EnvironmentTab.png")

1 + x
x^2

sqrt(x)
log(x)

# Create an object called x
x <- 1:5

# What do we have in the object `x`
print(x)

# I'm not sure. Which values are greater than one?
x > 1

# Let's square every value in `x`
x^2
# And we can find the square root of every value
sqrt(x)

opts_chunk$set(include = TRUE,
               eval = TRUE,
               fig.show = 'asis')
