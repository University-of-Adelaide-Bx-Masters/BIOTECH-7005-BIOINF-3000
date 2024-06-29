# Assignment 5 [*30 marks*] Due Friday 18 October 2024

## Assignment overview [*4 marks for general format*]

This assignment is based on the Local version control and GitHub repository practicals. 

The repository at https://github.com/davmlaw/2023_assignment_5_adelaide_uni_bioinfo is a simulation of an open source project with bugs.

This assignment involves performing actions on GitHub [*4 marks total for general format*]:

* Issue for high severity bug [*1 mark*]
* Issue for low severity bug [*1 mark*]
* Fork repo [*1 mark*]
* Pull request for fixing the issues [*1 mark*] 

These have timestamps on them, and these will be checked against the submission date.

Answers to theory questions should be submitted to MyUni as a file (plain text or PDF)

## Practical questions [*21 marks*] 

Please read the instructions in the initial README of the repository, and only modify the files assigned to you.

There are two kinds of changes/bugs:

* insertion of XXXXXXX
* changing some lower case words to UPPERCASE

These should be treated as two separate incidents, with different severities, and will be handled separately.

1. Raise an issue in the [repo](https://github.com/davmlaw/2023_assignment_5_adelaide_uni_bioinfo) for the high severity bug (XXXXXXX)
    + Giving a good title succinctly stating the problem and identifying the file (such that it is clear why it is distinct from other issues) [*1 point*]
    + Thank the maintainer and use a polite, kind and non-accusatory tone  [*1 point*]
    + Provide the person's name mentioned in the comment of the commit that introduced the error to your file [*1 point*]
    + Write a expected / actual for the lines affected, only mentioning the high severity bug (XXXXXXX) [*2 points*]
    + Add a comment that you plan on submitting a pull request [*1 point*]

(total: 6 points)

2. Raise an issue in the [repo](https://github.com/davmlaw/2023_assignment_5_adelaide_uni_bioinfo) for the low severity bug (UPPERCASE)
    + Giving a good title succinctly stating the problem and identifying the file (such that it is clear why it is distinct from other issues) [*1 point*]
    + Provide the date and Git hash of the commit that introduced the error to your file [*1 point*]
    + Write a expected / actual for the lines affected, only mentioning the low severity bug (UPPERCASE) [*2 points*]

(total: 4 points)

3. Fork the repository to your personal account.
    + Make sure you untick the "Copy the main branch only" - you want the other "stable" branch to also be copied

4. Fix the high severity issue (XXXXXXX), in the "main"" and "stable" branches
    + You will need to clone your private fork of the repo to your local machine.
    + Edit your data/ file, and remove the high severity issue (XXXXXXX) ONLY (ie DO NOT change the low severity (UPPERCASE) at the same time). There should only be 1 space between the words after removing the XXXXX [*1 point*]
    + Edit your doc/ file, which is in the format of a [changelog](https://keepachangelog.com/en/1.0.0/) adding a CHANGED entry for your fix, that links to the high severity issue above [*1 points*]
    + Write a descriptive commit message, which references the high severity issue above [*1 point*]
    + This bug needs to be fixed in the "stable" branch as well, so switch to that branch and cherry-pick the commit you just made [*3 points*]  
    + Push both branches to your repo

(total: 5 points)

5. Fix the low severity issue (UPPERCASE), in the "main"" branch ONLY
    + Edit your data/ file, and remove the low severity issue. [*1 points*]
    + Edit your doc/ file, which is in the format of a [changelog](https://keepachangelog.com/en/1.0.0/) adding a CHANGED entry for your fix, that links to the low severity issue above [*1 points*]
    + Write a descriptive commit message, which references the low severity issue above [*1 point*]
    + Push both branches to your repo

(total: 4 points)

6. Submit a pull request from your main/master branch
    + Comments should link to both issues and give a brief description of the change [*2 points*] 

(total: 2 points)
    
## Theory questions [*5 marks*]

1. How do you tell Git to stop tracking / ignore some files. Name 2 types of files you might want to ignore and why. [*1 point*]

2. What is the staging area in Git, and why is it used? [*1 point*]

3. Why would you want a long lasting branch (eg "stable" above). [*1 point*]
 
4. Why would you want a short lived one (ie a "dev" or "feature" branch) [*1 point*]

5. Why is it called a pull request? [*1 point*]

### Resources

- [How to write a good bug report](https://musescore.org/en/node/309537)
- [Issues - Actual vs Expected](https://medium.com/we-are-testers/chapter-2-how-to-write-useful-actual-and-expected-results-details-in-your-bug-report-10b83e5aaa75)
- [Git - branches in a nutshell](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell)
- [GitHub - about branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches)
- [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
- [GitHub - creating a pull request from a fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork)
