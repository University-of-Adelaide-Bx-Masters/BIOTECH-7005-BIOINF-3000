---
title: "Code comments"
author: "David Lawrence"
date: "2022-09-05"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


**"Programs must be written for people to read, and only incidentally for machines to execute."**

― Harold Abelson, Structure and Interpretation of Computer Programs 

Comments are parts of programming languages that are *removed* from the source code before execution. Comments are just for people.

In the early 80s a famous computer scientist named Donald Knuth took code documentation to the extreme and invented a programming paradigm called "Literate programming" - where you describe in detail what you are doing in English, and intersperse the human language with small amounts of code.

If this seems familiar, it's because it has been adopted by scientific computing in things like R notebooks / R markdown. In those type of environments, you can be as verbose as you like, and use embedded images, tables etc to lead readers through a story.

For the rest of the tutorial, I am only going to talk about comments in standard programming language source code.

## R

In the R language, comments start with "#" and last until the end of the line.

```
# This is a comment in R!
# Some programming languages have multi-line comments
# but with R, you just keep making more single-line comments
```

## Good comments

Good comments help a reader understand the code more quickly and correctly.

When creating code - you have all of this extra information in your head - background information, what you want to achieve, tradeoffs you made - that led you to write it the way you did.

If someone else (or you 6 months later) comes back to this code (say to look for a bug), they may not have all of this information.

Comments can help add context, and help the reader understand **why** a program is written the way it is.

Think - if you came back 6 months later with no memory of writing this - what would you wish to know?

## Bad comments

Bad comments can be useless - if they add nothing to the code.

```
i = 1
i = i + 1  # Add 1 to i
```

Bad comments can add clutter / noise to a program, without adding an useful information.

Maybe when you are just starting out it's ok to write notes for yourself, but as you get better and share things, aim for someone that is at least basically familiar with the language and standard libraries. They can look them up elsewhere, don't clutter up your code!

Comments can get out sync with the code. Comments can lie!

```
i = 1
i = i + 2  # Add 1 to i
```

A wrong comment is far worse than no comment.

Comments need to 'pull their weight' to be worth their costs in space and maintenance.

## Some R code

Here's some R code. What does it do?

```{r}

a = "GATTACA"
b = toupper(strsplit(a, '')[[1L]])
c = ifelse(b == 'A', 'T',
    ifelse(b == 'C', 'G',
    ifelse(b == 'G', 'C',
    ifelse(b == 'T', 'A', 'N'))))
d = paste(rev(c), collapse = '')
d
```

This code produces the reverse complement of a DNA sequence

```
GATTACA
|||||||
CTAATGT
```

The code is correct, but not very clear. Perhaps the biggest problem is it's not obvious what the variables a,b,c,d etc are used for.

## Adding some comments

Perhaps we can make it more clear, by adding some comments. The more the better!

```{r}
# this section takes a DNA sequence 'a', obtains the reverse complement and stores in 'd'
a = "GATTACA"  # DNA sequence
b = toupper(strsplit(a, '')[[1L]])  # split sequence into individual bases
c = ifelse(b == 'A', 'T',  # convert to complement
    ifelse(b == 'C', 'G',  # G is complement of C
    ifelse(b == 'G', 'C',  # C is complement of G
    ifelse(b == 'T', 'A', 'N'))))  # N is used as unknown base
d = paste(rev(c), collapse = '') # reverse and join with no separators, d = reverse complement
d # Print variable
```


I also added a comment about what the whole section does

## Better names

It looks like we are trying to fix bad names and structure with comments.

* The "section" could be made into a function - and given a name.
* The variable names could be made more descriptive

```{r}

# Code taken from https://stackoverflow.com/a/69991965
reverse_complement <- function (dna) {
    UNKNOWN_BASE = 'N'
    bases = toupper(strsplit(dna, '')[[1L]])
    complement = ifelse(bases == 'A', 'T',
        ifelse(bases == 'C', 'G',
        ifelse(bases == 'G', 'C',
        ifelse(bases == 'T', 'A', UNKNOWN_BASE))))

    paste(rev(complement), collapse = '')
}

reverse_complement("GATTACA")
```


This is not only much clearer but the function can be re-used. The complexity is moved into the function (you only need to know the name and call it, not manage the complexity and variables etc yourself)

The lesson here is to not use comments to cover up bad code.

If you see people making "sections" out of a script, eg "init" and "graphing" etc - they've already identified components that can be broken up into functions, or packages - then use those as names.

## Dead code

It's ok to comment out code while you're working on it (say to try out different implementations) - but generally it's a bad idea to keep "dead code" around. It takes up space and will slowly get out of date. 

Since you're using source control (right!) - don't be afraid you can always get it back (though you probably never will)

## Simplifying formulas

A common example is where you may be tempted to comment a formula due to its complexity.

```{r}
cars[cars$speed > 20 | cars$dist > 70, ]  # Find cars going too far or fast
```

But remember - you can always pull apart a formula and assign it to variables

This increases the number of lines, but makes each line simpler. You can then add information to each line (the **WHY**)

It is incredibly important to document "magic numbers" - where did they come from? How did you choose them?

```{r}
SPEED_LIMIT = 20  # town limit mph
INSURANCE_MAX_DISTANCE = 70  # standard cover
too_fast = cars$speed > SPEED_LIMIT
too_far = cars$dist > INSURANCE_MAX_DISTANCE
cars[too_fast | too_far, ]
```


## Always include an URL if you copy/paste code from online

It's good form to cite, provides context, and it's possible someone may point out a bug or provide a better answer in the future!

For instance, here's some lines that contain "stackoverflow" in my code:

```
$ find . -name "*.py" | xargs grep stackoverflow
./upload/views/views.py:            @see https://stackoverflow.com/a/65596550/295724 """
./library/file_utils.py:    """ By Michael Anderson: http://stackoverflow.com/a/12593795 """
./library/celery_utils.py:From: http://stackoverflow.com/a/23908345
./library/django_utils/django_file_system_storage.py:        From https://stackoverflow.com/a/32349636/295724 """
./library/django_utils/django_postgres.py:    """ From https://stackoverflow.com/a/63715608/295724 """
./library/django_utils/unittest_utils.py:    @see https://stackoverflow.com/a/46079090/295724
./library/django_utils/unittest_utils.py:        and contain the file asked. @see https://stackoverflow.com/a/51580328/295724 """
./library/postgres_utils.py:Code by Mike T - http://stackoverflow.com/a/8150329
./library/jqgrid.py:        # normally, but can use a trick see https://stackoverflow.com/a/45369944/295724
./library/utils.py:# From https://stackoverflow.com/a/26496899
./library/utils.py:    """ From https://stackoverflow.com/a/17246726 """
./library/utils.py:    """ https://stackoverflow.com/a/9754466/295724 """
./library/utils.py:    """ https://stackoverflow.com/a/22045226 """
./library/utils.py:        From https://stackoverflow.com/a/18035135 """
./library/database_utils.py:        From: https://stackoverflow.com/a/47542953
./snpdb/templatetags/model_helpers.py:By Fydo from http://stackoverflow.com/a/26614950
./snpdb/admin_utils.py:    From https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
./snpdb/tests/test_urls.py:        and contain the file asked. @see https://stackoverflow.com/a/51580328/295724
```

## Non-standard code

You should aim to write straight-forward code, ideally so that if someone else (or you in 6 months time) comes back - they think "that's how I would have done it" and can quickly understand and work with it.

If the code is overly-complex, or appears to do strange things, the reader is going to have to work at "getting in the head" of the writer. Or they may decide to re-write it.

This is why if you ever write non-standard/unexpected code you need to write why, for instance working around a broken library, a non-obvious fact about the biology or methods, or benchmark results showing why the obvious solution is too slow, etc.

## TODO comments

I use these when writing code, eg you're in the flow and write eg # TODO: Handle edge case error here. Ideally, you'll get them by the end of the day but maybe you decide to cut corners to move fast. This is taking on "technical debt" and the TODO documents this.

If it matters, raise an issue to fix the TODO in the issue tracker. I don't generally do this though if it's a purely technical thing and doesn't affect the correctness gg "this function and function X could probably be merged".

If I'm writing test code etc, I sometimes write "TODO: Don't check in!!!" - and scan for that when reviewing changes in Git. This has saved me a few times!

## Documenting functions

You should spend time and care writing comments at interfaces - ie edges between components such as functions or APIs. This is where you pass on information to people that don't know and don't want to know the internal implementation details of your code.

A lot of languages have a special feature when using an IDE (code editor) where they show a special kind of help text when calling a function.

In R - you can only do this in packages, where you put comments like this above the functions:

```
#' Reverse complement a DNA sequence
#'
#' This function returns the reverse complement to the passed sequence
#' @param dna string of bases ('AGCT')
#' @keywords dna nucleotides
#' @examples "GATTACA" => "CTAATGT"
#' reverse_complement()

reverse_complement <- function (dna) {
    # TODO
}
```

You can then use special tools (roxygen2) to read your source code and generate nice looking documentation. This is a great way to keep your documentation and code in sync.

See https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
