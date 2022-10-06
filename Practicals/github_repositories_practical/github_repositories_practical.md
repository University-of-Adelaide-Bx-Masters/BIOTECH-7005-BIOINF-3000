# GitHub repositories practical

## GitHub

Git is the command line tool (and version control database)

GitHub is company that runs git on the cloud to help you host a repository

A centrally hosted repo on the internet (with authorization) is useful as you can access it from many computers

They also provide many other features to help you run software projects like issue trackers, wikis and more

## Projects / Repositories

Here are some repositories of famous bioinformatics tools.

Go and view the pages, skim the README and the issue trackers. Open some random issues

[IGV](https://github.com/igvteam/igv) - paper in 2013, last commit 10 hours ago

[Samtools](https://github.com/samtools/samtools) - paper in 2009, last commit 21 hours ago

[BWA](https://github.com/lh3/bwa) - paper in 2009, last commit 14 days ago

Compare the release date of the paper vs the last commit.

Consider how different it must be from the first release, and whether the authors realised someone (or they themselves) would be working on the code 13 years later!

## Licences

Look to see if you can quickly spot the licence against the bioinformatics projects above.

Click on the link, and open it up and skim it. Open [TLDR legal](https://tldrlegal.com/) and search for GPL or MIT licences.

## Documentation

Github has excellent help -[GitHub docs](https://docs.github.com/en)

However, I tend to just Google "GitHub X" to find what I'm looking for.

## User security

GitHub handles [user authentication](https://docs.github.com/en/get-started/learning-about-github/access-permissions-on-github)

You can make organisations (which are basically groups that own projects rather than individual users).

This is a matter of choice (or lab policy) but I generally release work projects under my organisation as it's more professional, makes it easier for my coworkers to maintain it, and our lab seems more permanent than an individual job. You can always transfer it, and have the old URL redirect, but generally it's better to make URLs permanent.

## Access tokens

GitHub uses [Personal access tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) rather than passwords. To generate one of these for the VM:

* Click on your user icon in the top right
* Click "Settings"
* In bottom left menu option - "Developer Settings", then "Personal Access token"
* Generate new token, label it something like "bioinformatics VM", it needs repository access

If you have password access to your VM, and you are happy to store your password in plain text on the VM, then you can do so via:

```
git config --global credential.helper store
```

Otherwise, just leave the tab open so you can copy/paste it, or copy/paste it into a temporary text editor

## Projects

Create a new project
* 
* In the top right hand corner, click the plus symbol, then “new repository”
* Call it eg 'week9_github_practical'
* Choose a licence, add a README. Probably make it private

## Issues - searching

Go to [IGV issues](https://github.com/igvteam/igv/issues)

Notice how all you can see are the titles. Titles are extremely important.

Before raising a bug, you should always check it see if it has already been reported (ie don't create a "dupe")

Run a few searches. Experiment with clicking on Open/Closed. Notice how special tokens in the search form change. 

Again, notice how the title is critical for finding things as you can't see the contents (though you can search for it, eg "font" will bring up titles without "font" in it, as they have it in the issue)

## Issues - labels

As well as titles, labels can be used to organise issues into related themes. Read the documentation [Github - managing labels](https://docs.github.com/en/issues/using-labels-and-milestones-to-track-work/managing-labels)

* Create some labels for your project, with distinct colours - if you're unimaginative you could use "red", "green", "blue" with appropriate colours
* Raise some fake issues on your project, don't worry about content just make a lot of different titles (ie "problem X with red", "problem Y with red", "problem X with green", "problem Y with green")
* Tag the issues with appropriate keywords (view each individual issue, click label cog on right hand toolbar, select from autocomplete)
* Perform a search for keywords, and notice how the labels make it easier to find related 
* Try clicking on a label (this does a search by label)

## Issues - lifecycle

Issues move througha  lifecycle, ie:

raised -> fixed -> tested -> closed
raised -> tagged as "duplicate", "wontfix" -> comment about which other issue is a duplicate + close

* Find an issue, 
* Add "wontfix" label
* Add a comment at the bottom saying this issue is too hard / out of scope of the project and close the issue
* Go back to the search, notice it's gone

## Issues - writing tips

Even for your own projects - try and write it with as much information as possible (aiming to be standalone, so you don't need much outside context), so that a 3rd party can come in and make sense of things.

If it's clear - maybe it'll stop them raising a dupe. Maybe someone will fix it for you!

Depending on the type of bug, there's certain information you want to include

For crashes: stack traces and logs

For ones that involve complicated user steps - Try to find a way to “reproduce” the error consistently and produce clear step by step instructions. Ideally your instructions should reproduce the issue 100% of the time with a minimal amount of steps

For ones involving input files, perhaps try and produce a minimal, cut down one.

For unexpected behavior, or outputs, [Actual vs Expected](https://medium.com/we-are-testers/chapter-2-how-to-write-useful-actual-and-expected-results-details-in-your-bug-report-10b83e5aaa75) is very useful

## Milestones

Milestones are a way of managing releases and deadlines

For instance, you may have 20 issues, but your lab head says you have to stop coding and start writing up your paper in 2 weeks.

* Create a milestone called "Code complete for paper" with a deadline in 2 weeks time
* Assign a few issues to this (you will need to go into each individual issue, click milestone cog, select from autocomplete)
* Go back to issues, search by milestone

## Wiki

* Create a few wiki pages for your practical project (the one you made not forked)
* Use 2 levels of headers "#" and "##"
* Go to the wiki page, notice that it has built a hierarchical wiki browser for you out of your headers

# Simulated open source project repo

## Practical week 9 simulated open source repo

There's a [Week 9 GitHub repo](https://github.com/davmlaw/practical_week_9_github_repositories) which is a simulation of an open source project with bugs.

It has randomly generated data, so has the same structure (but different data) as [Assignment 5](https://university-of-adelaide-bx-masters.github.io/BIOTECH-7005-BIOINF-3000/Assignments/Assignment5.html)

## File history

From the project page, click on "data" then click on your "VM file". You can see the latest content, and the last commit.

* Click "history" in the top right hand corner to see all the commits.
* Go back and click "blame" to see commits that last touched a row

## Branches

Branches allow you to keep different versions of code (say, a "stable" vs "experimental/active development") and keep work separate.

You can split off branches, and merge branches into each other.

[Git - branches in a nutshell](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell)
[GitHub - about branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches)

* Click on the "branches" links on the practical repo

## Raising an issue

Raise an issue in the [Week 9 GitHub repo](https://github.com/davmlaw/practical_week_9_github_repositories) for the high severity bug (XXXXXXX)

Make sure the title includes your filename and what's wrong with it, ie "file doc/BIOINF_3000_2022-1 contains 'XXXXXXX'"

Remember to be nice - perhaps start the issue with "Hi, thanks for your project! I think I may have found an issue" - be humble and kind. You may be adding more stress to an overworked person trying to do altruistic work...

A bad issue raises stress but a good issue helps people, and makes it more easy for them to fix, and happier to do so.

Write [Actual vs Expected](https://medium.com/we-are-testers/chapter-2-how-to-write-useful-actual-and-expected-results-details-in-your-bug-report-10b83e5aaa75) for your file, only mentioning the high severity bug (XXXXXXX)

## Fork the repository to your personal account.

* You can now clone your private fork of the repo to your local machine.

## Fix the high severity issue (XXXXXXX), in the "main"" and "stable" branches

* Edit your data/ file, and remove the high severity issue (XXXXXXX) ONLY (ie DO NOT change the low severity (UPPERCASE) at the same time). There should only be 1 space between the words after removing the XXXXX
* Write a descriptive commit message, which references the high severity issue above
* This bug needs to be fixed in the "stable" branch as well, so switch to that branch and [cherry-pick](https://git-scm.com/docs/git-cherry-pick) the commit you just made
* Push both branches to your repo (you may need to handle access tokens here (see top of this file))

## Your forked repo

Go to your forked repo page.

* Notice your files last updates have changed
* Notice it says "X commits ahead of the old repo"

## Submit a pull request

* Comments should link to both issues and give a brief description of the change
* Go to the original repo (not your fork) - notice that the "pull request" tab now shows some numbers. Click it and find your pull request
* Go to your issue - you should be able to see the linked pull request



