---
title: "Reproducible Science Practical"
author: "David Lawrence"
date: "2024-09-12"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# R studio

This is generated from an R markdown file, open it in R studio:

BIOTECH-7005-BIOINF-3000/Practicals/reproducible_science_practical/reproducible_science_practical.rmd

# Hashes

### Command line

```{bash}
echo "Hello, World!" | md5sum
```


```{bash}
md5sum data/hashing/hello_world.txt
```

### File integrity

We copy billions of bytes all the time. How can we be sure a bit hasn't been flipped somewhere?

Hashes can be used to verify files moved across networks and disks. Sometimes you'll see hashes on FTP sites, or download pages. You can use these to verify your files haven't been corrupted.

md5sum has a short cut to perform checks across multiple files:


```{bash}
md5sum -c md5sum_hashes.txt || true # don't die on error
```

### Hashes in R

Most languages have a way to calculate the common hash functions.

In R, install the "digest" package:


```{r}
install.packages("digest")
```

Then you can calculate it from a string, or a file:


```{r}

library(digest)

# Using same file as md5sum above
file_path <- "data/hashing/hello_world.txt"
md5_file_hash <- digest(file = file_path, algo = "md5")
print(md5_file_hash)
```

Or with a string:

```{r}
string <- "Hello, World!\n"
md5_hash <- digest(string, algo = "md5", serialize=FALSE)
print(md5_hash)

```


# Diff

The Unix diff command compares two files *line by line* and displays the differences between them. It highlights which lines need to be added, deleted, or changed to make the files identical. The character at the start of each line means:

"<" lines in first file 
">" lines in the second file

It works great on simple files:

```{bash}
#!/bin/bash

echo "List 1:"
cat data/simple_diff/shopping_list1.txt

echo
echo "List 2:"
cat data/simple_diff/shopping_list2.txt

echo
echo "The difference is: "
diff data/simple_diff/shopping_list1.txt data/simple_diff/shopping_list2.txt || true  # we don't want to return an error here



```

However, compare the same data that's unsorted:


```{bash}
#!/bin/bash

echo "List 1:"
cat data/simple_diff/shopping_list1_unsorted.txt

echo
echo "List 2:"
cat data/simple_diff/shopping_list2_unsorted.txt

echo
echo "The difference is: "
diff data/simple_diff/shopping_list1_unsorted.txt data/simple_diff/shopping_list2_unsorted.txt || true  # Don't want to return error code

```

While this output is technically correct - it could be used to perform modifications to make the files identical, it's probably not what you want.

### Diff exercise

I generated 100 fake scientific papers, and put them in a text file in "data/scientific_papers"

There's a second file there, with 99 of the same papers, but 1 is missing. *Can you find it?*

```{bash}
#!/bin/bash

# Show that the same lines are in both 
grep -n holistic data/scientific_papers/*
```

# Determinism

The code you have written so far has all been deterministic - you could run the program a million times and it would always produce exactly the same answer.

Determinism is useful for reproducibility, but also for debugging. Remember the sorted vs unsorted diff above? If a code change causes a single line to change, it's much easier to find it if the data appears in the same order each run.

### Race conditions

When running multiple things in parallel, things like CPU scheduling or disk I/O can affect run time, so the order wil be non-deterministic.

Consider the Dad-joke bash script below, which executes 5 jobs taking between 0 and 5 seconds to execute:


```{bash}
#!/bin/bash

function joke1 {
  sleep $[ ( $RANDOM % 5 )  + 1 ]s
  echo "1. Knock Knock"
}

function joke2 {
  sleep $[ ( $RANDOM % 5 )  + 1 ]s
  echo "2. Who's there?"
}

function joke3 {
  sleep $[ ( $RANDOM % 5 )  + 1 ]s
  echo "3. Race Condition"
}


function joke4 {
  sleep $[ ( $RANDOM % 5 )  + 1 ]s
  echo "4. Race condition who?"
}

function joke5 {
  sleep $[ ( $RANDOM % 5 )  + 1 ]s
  echo "5. Race condition who depends on execution order!"
}

{ 
  joke1 &
  joke2 &
  joke3 &
  joke4 &
  joke5 &
  wait
} # | sort
```

### Random numbers (RNG)

Some algorithms use randomness, eg: Monte-Carlo simulations, stochastic gradient descent, randomly sub-sampling a file, or an aligner randomly assigning a read between two equally scoring locations. 

While you could sort this as per above, there's another way, that takes into account how computers generate random numbers.

CPUs generate "pseudo" random numbers, they use functions that take a number as input, and produce another number that looks really different. You can produce a sequence of numbers by feeding the previous number back in.

But you need a starting number - this is called a random "seed"

What's a good seed? Something a CPU can get that's random - current time in milliseconds, temperature of a CPU. Someone even took photos of a lava lamp...

R, like most RNGs, set the seed for you, so calls to random are non-deterministic by default. Here's how to generate 5 random numbers in R:


```{r}
runif(5)
```

Now run that multiple times. Are you getting the same number?

```{r}
runif(5)
```


```{r}
runif(5)
```



```{r}
random_numbers <- runif(100000)  # 100k random numbers
hist(random_numbers, main = "Histogram of Uniform Random Numbers")
```

```{r}
# Now we set the seed
set.seed(42)
runif(5)
```

```{r}
# Now we set the seed
set.seed(42)
runif(5)
```


What does an RNG look like?  You are not expected to understand this - just see that it's not super complicated to make somewhat random numbers...


```{r}
# Linear Congruential Generator (LCG) RNG
simple_rng_lcg <- function(n, seed = 123) {
  a = 1664525
  c = 1013904223
  m = 2^32
  random_numbers <- numeric(n)  # Initialize a vector to store random numbers
  X <- seed  # Starting with the seed
  
  for (i in 1:n) {
    X <- (a * X + c) %% m  # LCG formula
    random_numbers[i] <- X / m  # Normalize the result to get a number between 0 and 1
  }
  
  return(random_numbers)
}

# Example: Generate 5 random numbers using the LCG method
random_numbers <- simple_rng_lcg(5)
print(random_numbers)
```



```{r}
# This isn't as good as the default R one, but it's not too bad...
random_numbers <- simple_rng_lcg(100000)
hist(random_numbers)
```


# Tool / Library versions

To build modern software you often use large number of 3rd party libraries (eg built by random people and put on the internet). Examples include CRAN (R) and Pip (Python)

These libraries can change over time. If you don't specify what version to use, your script can suddenly stop working upon any one of the libraries you use (or libraries they use...) being updated.

A way to fix this is to explicitly list the version of a library used.

### Finding versions

For a tool - usually "--version". You may want to clean it a bit:


```{bash}
R --version | grep "R version"
```

In R, you can call packageVersion:

```{r}
packageVersion("digest")
```

### Pinning versions

How do we ensure R is loading the correct version?

See [renv](https://rstudio.github.io/renv/articles/renv.html) which allows you to load specific versions of packages


```{r}
install.packages("renv")
```


*Same code, different outputs!*

An example of a popular library function changing behavior across versions is dplyr

In pre-1.0 dplyr, after calling summarize(), the data is automatically ungrouped.

After 1.0 - the data remains grouped after summarization unless you explicitly ungroup it using ungroup()

*note: The below takes a while so you probably don't want to do it in class*

*note: It's generally a bad idea to mix different versions in the same script!*



```{r}
renv::install("dplyr@0.8.5", prompt=FALSE)

library(dplyr)
ddplyr_version = packageVersion("dplyr")
print(ddplyr_version)

```

Now, install a newer version:


```{r}
renv::install("dplyr@1.0.5", prompt=FALSE)

library(dplyr)
ddplyr_version = packageVersion("dplyr")
print(ddplyr_version)

```


On your own project (not here) sometime, try running:

```{r}
# renv::init()
# renv::snapshot()
```

This will dump out all the versions currently loaded in the systme into a file called "renv.lock" which looks like eg:

```
    "digest": {
      "Package": "digest",
      "Version": "0.6.37",
      "Source": "Repository",
      "Repository": "CRAN",
      "Requirements": [
        "R",
        "utils"
      ],
      "Hash": "33698c4b3127fc9f506654607fb73676"
    },
```

# Unit testing

Can you write software right on the first go? I can't!

How do you know if it gives you the right answer? You test it.

"Unit testing" is a form of testing where you write code to run and execute your code. It is aimed at testing small components (units) like a few functions.

Running automated unit tests regularly gives you confidence you haven't broken everything by making a change somewhere.

Some people to "test first development" where you define the expected answers before you even 


```{r}
library(testthat)

# Define the function that counts nucleotides
count_nucleotides <- function(dna_string) {
  nucleotides <- c("A", "C", "G", "T")
  counts <- sapply(nucleotides, function(n) {
    matches <- gregexpr(pattern = n, text = dna_string)[[1]]
    if(matches[1] == -1) {
      return(0)
    } else {
      return(length(matches))
    }
  })
  names(counts) <- nucleotides
  return(counts)
}

# Now let's write a test for this function
test_that("count_nucleotides correctly counts nucleotides", {
  dna_string <- "ACGTACGTACGTACGTACGT"
  expect_equal(count_nucleotides(dna_string), c(A=5, C=5, G=5, T=5))

  dna_string <- "AAACCCGGGTTT"
  expect_equal(count_nucleotides(dna_string), c(A=3, C=3, G=3, T=3))

  dna_string <- ""
  expect_equal(count_nucleotides(dna_string), c(A=0, C=0, G=0, T=0))
})
```

# Dependencies

Many things in life have "dependencies" - you need to do one step before another.

In computing, examples include:

* You need to install a library for a program to work
* Map reads before calling variants (run pipeline step as input to later one)
* You need to build a sub component (code module) before another that uses it

These dependencies can be represented by a directed acyclinc graph.

### Directed Acyclic Graphs

If you don't already have the package:

```{r}
# install.packages("igraph")
```

Show DAG

```{r}
# Load the igraph package
library(igraph)

# Create a DAG representing a bioinformatics pipeline
edges <- c(#"Raw Reads", "Quality Control", 
           # "Quality Control", "Alignment",
           "Raw Reads", "Alignment",
           "Alignment", "BAM Indexing",
           "BAM Indexing", "Calculate gene coverage",
           "BAM Indexing", "Variant Calling",
           "BAM Indexing", "CNV detection",
           "Variant Calling", "Merge calls",
           "CNV detection", "Merge calls",
           "Merge calls", "Annotation",
           "Annotation", "Filtering and prioritisation")

# Create a graph object
g <- graph(edges, directed = TRUE)

# Define the layout
l <- layout_as_tree(g, root=1, mode="out")

# Plot the graph
plot(g, layout=l, 
     edge.arrow.size = 0.5, 
     vertex.color = "lightblue", 
     vertex.size = 20, 
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     edge.curved=0.2)

```

You'll notice we didn't explicitly tell the graphing library where to draw the nodes and edges.

We *defined the dependencies* then allowed the library to process the directed acyclic graph, then draw the image.

### Workflows

Some workflow tools allow you to define dependencies, then it will auto manage things for you, such as:

* On errors - run from last successful step, rather than the whole pipeline again
* Know what files need to be updated if a step changes (eg if you switch out variant caller - everything below that needs to be re-run but not above)

Examples include:

* [SnakeMake](https://snakemake.readthedocs.io/en/stable/) - which is based on Python
* [NextFlow](https://www.nextflow.io/) - which is based on Groovy (scripting language based on Java) - see [nf core pipelines](https://nf-co.re/pipelines/) - great collection of pre-assembled Bioinfo pipelines
* Make - you may run into this building Unix tools but probably don't want to use it yourself

# Managing environments

[Anaconda](https://docs.anaconda.com/working-with-conda/environments/)

Containers - [What is Docker - 15m video](https://www.youtube.com/watch?v=rOTqprHv1YE)

VMs - You are already using a virtual machine! There are free tools to set one up yourself! See [VirtualBox](https://www.virtualbox.org/) among others

# Further reading

It'll take a few days each to understand workflow languages, or containers, but you can get very quick and easy wins just by:

* Pinning your software versions  [renv](https://rstudio.github.io/renv/articles/renv.html) 
* Clear and re-run your notebooks before trusting the results or sharing them. [paper](https://link.springer.com/article/10.1007/s10664-021-09961-9)
