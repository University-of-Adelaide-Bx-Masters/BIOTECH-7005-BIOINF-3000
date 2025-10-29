---
title: "GitHub Extra"
author: "David Lawrence"
date: "2025-10-29"
output: html_document
---

Some extra material on GitHub / Open source for students that finish early

## Licences

On one of the bioinformatics project pages above, look in the top right for the licence link. If LICENCE is a standard file (as per IGV), it will show eg "MIT Licence"

Click on the [IGV licence link](https://github.com/igvteam/igv/blob/master/license.txt), and open it up and skim it quickly.... it's a lot of legalese

There's a great site called [TLDR legal](https://tldrlegal.com/) that simplifies things, eg:

[MIT](https://www.tldrlegal.com/license/mit-license) - A short, permissive software license. Basically, you can do whatever you want as long as you include the original copyright and license notice in any copy of the software/source

[GNU General Public License v3](https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3) (GPL-3) - You may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.

**Task:** Please open the links above, and read through them

## GPL is a viral licence

You can use GPL libraries in your code without worrying about the licence, however

* You can only copy/paste code from a GPL project if your project is also GPL
* You can’t make it more restrictive (not release source code)
* You also can’t downgrade the licence  (release as MIT!)

Both Android (Linux) and Apple iPhone (BSD) are based on open source kernels

Linux is GPL - so any changes to the kernel by Google or Samsung must be released to the world.

Apple can take BSD, add their own work and not share it.

Notice for both that companies are making money reselling open source code (this is totally fine!)

Generally: If all you care about is use (don’t care if people spin off new projects, or if they give back) then go for MIT. Due to academic norms, they would be expected to cite it.

If you don't like the idea of people using your code without contributing, use GPL

**Discussion:** Licences - your thoughts?

## Documentation

Github has excellent help -[GitHub docs](https://docs.github.com/en)

However, I tend to just Google "GitHub X" to find what I'm looking for.

# Open Source software

Github is a place to find software, but if you're doing it right, GitHub is a social media site for programmers

Some use the term "free software" but you should think of free as in "free speech" not "free beer" - the most important part is that you can free to see how it works, free to modify it, if it's broken.

You can get paid a lot of money writing code... so why do people work on code then give it away?

* Usually start from "Scratching your own itch" - ie solve a problem they themselves had.
* "I’m sick of doing this repetitive task - I’m going to write a way to automate it"
* "All of the existing software sucks, I’m going to write something that does it properly!"

Often corporate programming is not that fun, and they want to do something different, solve a problem “the right way” (where engineering effort does not have an immediate source of income) or just work with personal freedom.

* Almost all open source programmers have “day jobs”
* Some companies sponsor development, eg some Linux developers are funded by Intel, RedHat, Samsung, Google, IBM etc
* Some developers consult, provide corporate training, write books, give talks etc
* Academic (while academia provides grants for projects which can spin off open source software, it is very rare to provide funding to maintain software)

![xkcd cartoon](https://www.explainxkcd.com/wiki/images/d/d7/dependency.png)

"There’s a reason why open source projects are called “projects,” rather than just code. While code is the final output of a project, the term “project” refers to the entire bundle of community, code, and communication and developer tools that support its underlying production."

- Nadia Eghbal, "Working in Public"

Software is very iterative. If you make something and put it out there. People use it (not in the way you expect!) then start making requests, asking for help, etc

If your project is successful you can often end up spending a lot of time managing it, and doing boring things like user support rather than the fun stuff (programming cool new features and doing science)!

"Running a successful open source project is just Good Will Hunting in reverse, where you start out as a respected genius and end up being a janitor who gets into fights”

- Byrne Hobart

But it is rewarding, and from a scientific career perspective - people are more likely to use and cite good, working tools!

Professionally - it also helps with exposure and networking. When hiring for jobs, I love to see someone list their GitHub account, you can see what they've done, and also have a glimpse into how they work (are they helpful and polite?)

## Advantages of Open Source

* Open (no hidden spyware, security issues etc).
* In science it’s important that you understand what’s happening!
* Code is the ultimate documentation (papers are out of date at submission time)
* Innovation - the corporate world is often somewhat stifling. You can create an open source project just because you have a cool idea

## Problems in Open Source

* Projects can get abandoned when the fun is gone
* Maintenance and bug fixes are less fun than writing new code 
* Writing documentation is less fun than writing new code

Generally - fun driven development means boring jobs don’t get done.

**Discussion:** Has anyone had good or bad experience with open source software? 

## Issues

If you find a bug in some software, please report it by raising an issue.

Check that hasn't been fixed (ie try latest version) and hasn't already been reported. Ideally work for a bit to narrow it down and make it simple to reproduce.

Raising an issue if you find a bug is a useful contribution!

Suggesting a new feature you'd like is also useful, though it's a lot easier to think of new features to do them, and not everyone shares your priorities.

There are a lot more people raising issues than fixing bugs and closing them

Thus, number of issues continues to grow while the software is being used (and it will probably be abandoned before they are all closed)

If a developers says "PR welcome", that's a polite way of saying "I have other things to do, if you want this, do it yourself!"

## Bug fix / feature timeline

In order for a bug to be fixed, the following has to happen:

* Bug is discovered
* Issue raised
* Developer reproduces bug or Agrees feature is worthwhile
* Developer fixes bug or writes feature
* Code is merged

The first steps are valuable, but developers are busy.

You can wait a loooooong time for a bug fix.

[Here's a bug I raised](https://github.com/Ensembl/ensembl-vep/issues/1023), that took 812 days to fix. I am still super pumped.

Why? I haven't written Perl in over 10 years, so had to wait for someone else to do it.

## 2 users 1 bug

Two people encountered the same bug, that had a simple 15 line fix. One had to wait 2,467 days for a fix, while the other one only had to wait 5.

* [Issue raised](https://github.com/brentp/cyvcf2/issues/17) Aug 3, 2016
* 1st developer reply Aug 3, 2016 - 0 days
* [Another user mentioned bug is still active](https://github.com/brentp/cyvcf2/issues/17#issuecomment-1529468360) - May 1, 2023 - 2462 days later
* [Pull request](https://github.com/brentp/cyvcf2/pull/262) May 4, 2023 - 3 days later
* Merged into code base May 6, 2023 - 2 days later

What's the lesson here? If you want a fix fast, do it yourself!

Who fixes open source code? People like me, and people like you!

## Forking

You generally don't have permission to modify other people's projects, so how do we modify it, or send them changes? Start by forking.

Forking is a way to create a personal copy of someone else's repository on GitHub. When you fork a repo, you're essentially duplicating it under your GitHub account. This allows you to freely experiment with changes without affecting the original project.

Steps for Forking on GitHub:

    Navigate to the repository you want to fork.
    Click the "Fork" button, usually located at the top right corner of the repo's page.
    Choose where you want to fork the repo. Usually, it's your personal GitHub account.

What Happens After Forking:

    You get a new repository under your account that's identical to the original repo at the time of forking.
    Your fork exists independently, meaning changes to the fork won't affect the original repo, and vice versa.

Common Use Cases:

    Contributing to a Project: Fork the repo, clone your fork locally, make changes, push back to your fork on GitHub, and then create a pull request to the original repo.
    Personal Experiments: You can freely modify your fork without worrying about the original project.

**Task:** Fork a bioinformatics project of your choice, then clone it to your local VM

## Modfiying tools or libraries

I always raise an issue before starting a fix, to confirm it is real, and 

When making bug fixes, I like to write automated unit tests that reproduce the bug, and would fail with the current code.

Then I make the fix, and verify that the new tests pass. I make sure all of the existing tests still pass, to show I haven't accidentally broken anything.

I like to use my code for a while to thoroughly use it myself for a few days under real world conditions, this means installing it to the system, so I can use the library from my code, or my modified tool in my pipelines.

Then, when I'm happy, I'll try to get it merged into the original project

**Question:** I have selfish reasons to get my bug fixes merged into the project. What are they?

## Pull requests

Since you don't have permission to write to other people's projects (ie you can't push), you need them to pull your changes. This is called a pull request.

If a pull request is still open, and you push more commits to your fork, it adds onto the currently open pull request.

I like to create a branch in my repository for just that fix, then make the pull request from there.

Unit tests drastically increase the chance of your merge being accepted.

This allows me to make multiple pull requests, separated out per issue. This makes it easier to validate (rather than a giant pull request that fixes 4 independent issues)

I put all of the fixes in the main/master branch - so if I need to replace the library/tool - this is what I use until the original project merges the fixes and makes a release.

## More examples

Developers have busy lives, sometimes it can take a while for the merges to get in. I try to gently remind people, and fulfill any requests they have.

Be patient, don't be rude, and just keep using your forks for personal use.

**Cutadapt:**

* [Issue raised](https://github.com/marcelm/cutadapt/issues/64) October 30, 2012
* 1st developer reply December 04, 2012 - 35 days
* [Pull request](https://github.com/marcelm/cutadapt/pull/3) December 14, 2012 - 10 days
* Merged into code base December 17, 2012 - 3 days

**HTSeq:**

* [Issue raised](https://github.com/htseq/htseq/issues/63) Jun 9 2023
* 1st developer reply Jun 9 2023 - 0 days
* [Pull request](https://github.com/htseq/htseq/pull/65) Jun 13 2023 - 4 days
* Merged into code base Aug 1 2023 - 49 days (had to fight the automated build system)

**Pronto**

* [Issue raised](https://github.com/althonos/pronto/issues/176) - Jun 21, 2022
* 1st developer reply - n/a
* [Pull request](https://github.com/althonos/pronto/pull/179) Jun 23, 2022 - 2 days
* Merged into code base Jun 29, 2022 - 6 days

**HGVS**

* [Issue raised](https://github.com/biocommons/hgvs/issues/447) - May 30, 2017
* [Pull request](https://github.com/biocommons/hgvs/pull/660) - June 30, 2023 - 2,222 days later
* Merged into code base - Sep 11, 2023 - 73 days


