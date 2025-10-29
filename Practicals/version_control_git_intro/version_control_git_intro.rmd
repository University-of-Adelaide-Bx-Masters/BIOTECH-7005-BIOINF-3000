---
title: "Version Control - Git Intro"
author: "David Lawrence"
date: "2025-10-19"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Vocab review

Make sure you are familiar with these terms!

* Diff
* Hashing

## Git help

Like most command line programs, git has command line help. Run:

```
git --help
```

This shows an overview of sub commands. It's short, read it all!

```
start a working area (see also: git help tutorial)
   clone     Clone a repository into a new directory
   init      Create an empty Git repository or reinitialize an existing one

work on the current change (see also: git help everyday)
   add       Add file contents to the index
   mv        Move or rename a file, a directory, or a symlink
   restore   Restore working tree files
   rm        Remove files from the working tree and from the index

examine the history and state (see also: git help revisions)
   bisect    Use binary search to find the commit that introduced a bug
   diff      Show changes between commits, commit and working tree, etc
   grep      Print lines matching a pattern
   log       Show commit logs
   show      Show various types of objects
   status    Show the working tree status

grow, mark and tweak your common history
   branch    List, create, or delete branches
   commit    Record changes to the repository
   merge     Join two or more development histories together
   rebase    Reapply commits on top of another base tip
   reset     Reset current HEAD to the specified state
   switch    Switch branches
   tag       Create, list, delete or verify a tag object signed with GPG

collaborate (see also: git help workflows)
   fetch     Download objects and refs from another repository
   pull      Fetch from and integrate with another repository or a local branch
   push      Update remote refs along with associated objects
```

You can get help for a subcommand via:

```
git clone --help
```

This opens a man (manual) page, and is the equivalent of ```man git-clone```

Scroll through the man page with up/down keys, then page up/page down keys.

Press q to quit the man page.

## Git setup / config

Git tracks which user did what, so you need to identify yourself.

Copy/paste this into the terminal, EDIT IT WITH YOUR ACTUAL INFO then run:

```
# You need to modify things below to your actual info!
git config --global user.name "Your name"
git config --global user.email "your.name@student.adelaide.edu.au"
git config --global core.editor "nano -w"
```

This is stored in your home directory and so works on all git repositories (global)

To see it:

```
cat ~/.gitconfig
```

## Local repository

If you have already have the repo from a previous tutorial, you need to update it:

```
cd BIOTECH-7005-BIOINF-3000 # Wherever you cloned it before
git pull  # fetch changes from GitHub and update your local repository
# If you get errors due to changes from the last tutorial, perhaps it is easier to just delete it and clone it again.
```

If you haven't cloned the repo before, do so via:

```
# Clone a repo
git clone https://github.com/University-of-Adelaide-Bx-Masters/BIOTECH-7005-BIOINF-3000

# Go into the repo
cd BIOTECH-7005-BIOINF-3000
```

## Your repo

When you clone a repo - you make a local copy of the whole version control database (history, deleted files etc)

Remember this is your own personal repository, you can do whatever you like! You are not modifying the GitHub repo unless you send your changes up, and you don't have permission to do that, so don't worry!

The files you see in this directory, the "normal" files, are not Git’s real copy of your project. To see them look in the .git directory.

```
find .git | less
```

Compare the size with the total repo:

```
du -sh .git
du -sh
```

The .git directory contains all of the changes ever made to files, git can move files through their history depending on what changes it applies

The "normal files" are just temporary ones made by git from that database (or edits that you have yet to tell Git about)

You're looking at the latest copy - ie files made by applying all of the changes (diffs)

## Changelog / History

A git repository is a sequence of changes, you can see a history of these via:

```
git log
```

This opens an interactive log (shows a ":" at the bottom) of commits in reverse chronological order (ie latest commits first). Press the up/down arrow keys to scroll, or page up and down. Press q to quit.

Notice that at the top of every commit is a hash - which is the unique label you can use to refer to that change. An example:

```
commit 97aba47a9e8a5db9fd0960df025fb11eeb20138f
Author: Dave Lawrence <davmlaw@gmail.com>
Date:   Tue Sep 6 03:10:15 2022 +0930

    Fix main README link to coding comments
```

At the top is a commit hash. This is the label for this change.
There's also the user, date and a change comment.

There are lots of options, eg to show a summary of changes.

```
git log --stat --summary
```

Or to see a visual of branching/merging:

```
git log --graph
```

Sample output:

```
| * |   commit 4ac3e07643488adad8a5a4462c756503cc5e975a
| |\ \  Merge: 6dea320 bbef0f4
| | | | Author: Dave Adelson <david.adelson@adelaide.edu.au>
| | | | Date:   Wed Sep 22 09:41:55 2021 +0930
| | | |
| | | |     Merge pull request #27 from University-of-Adelaide-Bx-Masters/NewTranscriptome
| | | |
| | | |     update link to merged major project document
| |\ \  Merge: 6dea320 bbef0f4
| | | | Author: Dave Adelson <david.adelson@adelaide.edu.au>
| | | | Date:   Wed Sep 22 09:41:55 2021 +0930
| | | |
| | | |     Merge pull request #27 from University-of-Adelaide-Bx-Masters/NewTranscriptome
| | | |
| | | |     update link to merged major project document
| | | |
| * | |   commit 6dea320691f4e3d5d5ec1e3b02c1c84db46f5ec6
| |\ \ \  Merge: 8a40db2 5b9a8b1
| | | | | Author: Dave Adelson <david.adelson@adelaide.edu.au>
| | | | | Date:   Wed Sep 22 09:33:33 2021 +0930
| | | | |
| | | | |     Merge pull request #26 from University-of-Adelaide-Bx-Masters/NewTranscriptome
| | | | |
| | | | |     New transcriptome
| | | | |
* | | | |   commit cff75d103c70161e44decc9bdf59f4cd8274cfeb
|\ \ \ \ \  Merge: 8a40db2 9d41f5e
| |/ / / /  Author: David Adelson <david.adelson@adelaide.edu.au>
|/| | | /   Date:   Wed Sep 22 09:48:38 2021 +0930
| | |_|/
| |/| |         Merge branch 'NewTranscriptome'
| | | |
| * | | commit 9d41f5e071da6e179660adb8a124beb2b190c2c6
| | |/  Author: David Adelson <david.adelson@adelaide.edu.au>
| |/|   Date:   Wed Sep 22 09:44:34 2021 +0930
| | |
| | |       deleted html
| | |
| * | commit bbef0f42d9e7efd2be26cd3d281f55dfcfcfef01
| |/  Author: David Adelson <david.adelson@adelaide.edu.au>
| |   Date:   Wed Sep 22 09:36:20 2021 +0930
| |
| |       update link to merged major project document
```

Or a more concise graph:

```
git log --all --decorate --oneline --graph
```

You're looking at a change log, a log of changes (diffs) applied in a sequence, where each is labelled with a hash. Lets look at the changes:

## Show / Diffs

To see the details of a particular commit, copy/paste a hash, eg:

```
git show 23d004b56503d2aaccc5f6e39993a8d10aee334c
```

To see the differences between two commits, you can use diff:

Show is the same as the difference between a commit and the previous commit:

```
git diff 11a25fea8e15523049896381de9312501053843d 23d004b56503d2aaccc5f6e39993a8d10aee334c
```

but you can eg look across 3 commits:

```
git diff 141d48db593f24783169062b1e94c4259bacfaeb f6cec9c27eea6ef854f754f1adbe3c3310cb708d
```

Exercise:

Use git log to find 2 commits, and run git diff between them.

Now add a few commits before/after it to see how the size of the diff grows

## HEAD - the latest commit

It's changed now, but when writing this the top commit was: 

```
commit f6cec9c27eea6ef854f754f1adbe3c3310cb708d (HEAD -> master, origin/master, origin/HEAD)
Author: dadelson <david.adelson@adelaide.edu.au>
Date:   Mon Sep 12 12:33:39 2022 +0930

    updated assignment 4
```

This commit has some extra stuff after it. These are special labels that are automatically created.

**HEAD** - is the end of a branch, ie the latest commit
**master (or main)** - is the default name of a branch (we'll come back to this)

To see find latest commit:

```
git log -1
```

Copy/paste the hash

```
git show <HASHCODE>
```

HEAD is simply a shortcut to the latest hash:

```
# Should show the same output as above
git show HEAD
```

Now get the 2nd last commit:

```
git log -2
```

copy 2nd to last hash

```
git diff paste_second_to_last_hash_here HEAD
```

This compares that commit with head

You use HEAD so often that it is the default argument in Git. Ie if you don't specify a commit, you are assumed to be comparing it to Git. Ie:

```
git diff paste_second_to_last_hash_here # Implicitly comparing to HEAD
```

**Shortcuts**

Every commit usually has one "parent" commit which points to the previous state, some shortcuts:

```
# These also work on normal hashes not just HEAD
git show HEAD^  # to see the parent of HEAD
git show HEAD^^ # to see the grandparent of HEAD
git show HEAD~4 # to see the great-great grandparent of HEAD
```

```
# See last 2 commits:
git diff HEAD^^  # Implicitly comparin to HEAD
```


## Moving through history

Git stores (in .git) all of the changes ever applied to the files in the repo. The "normal files" (under management of Git) are created by git by replying changes.

To move through history - you want to replace the current files (generated by running all diffs up till the latest, ie HEAD) with ones that only have the changes done up until a certain time (commit)

So - to see a snapshot of how the repository looked at a commit - you can "check out" that commit:

```
git checkout 141d48db593f24783169062b1e94c4259bacfaeb
```

All files (under management in the repo) are as they were at that time.

```
# You can copy this version of the file somewhere else
cp README.md /tmp
```

If you run "git log" now - you can see that the final commit is not pointing to "master" anymore.

To get back to the latest commit, run:

```
git checkout master
```

You now have 2 copies of the file - and can compare them

```
# Normal diff, not git diff
diff /tmp/README.md README.md
```

## Restoring deleted files

```
git checkout master  # Just in case we forgot where we are
# Delete the README file
rm README.md
```

Now run:

```
git status  # Shows file was deleted
git diff  # Shows patch
```

The file has been deleted (all rows with "-" in front)

Remember you can look for help - Google [git restore deleted file](https://stackoverflow.com/a/21307473)

So you can restore the file via:

```
git checkout README.md
```

## Something scary

Feel free to sit these ones out... I don't have any assignments to accidentally delete!

```
# Make sure you are in the right directory!
# Then run:
rm -rf *
```

Whoops! Now Google for "git restore all deleted files" - the top result is from [stack overflow](https://stackoverflow.com/a/26892936)

You should be able to restore everything via:

```
git ls-files -z -d | xargs -0 git checkout --
```

Note: I didn't know how to do this before looking on Stack Overflow.

Q: What would have happened if I deleted the .git directory?


## Modifying a file

* Go and edit the README.md file and delete lines, add lines, just totally trash it
* To see what you changed (can you remember?)

```git diff``

* Imagine if you had to restore this by manually undoing those changes?
* Revert all changes via:

```git checkout README.md```

## Committing a file

* Modify the README.md and add a section called "students" and add your name

Or, modify the file in any other way you want.

Run

```
git diff
```

To see the lines you added.

Commit the changes - the "-m" adds a message. You should write good messages!

```
git add README.md
git commit -m "Added student section"
```

Now run ```git log```

You can see we're out of sync with "origin/master" (which is the remote - ie GitHub)

## Blame Game

Run

```git blame README.md```

Scroll through the file until you get to the student section. You should see your own name, indicating you were responsible for the last changes to this line. They're onto you!

## Destroying commits

The latest change was put into our repository. We haven't "pushed" it up to GitHub yet, so it's only here on our local repo

If it's not secret information (or embarassing) you should probably just fix things like this via a "soft" restore:

* Restore the file from before the change
* Check in that file over the top of the change

This will leave 2 commits, the change and undoing the change.

But - you can delete the commits - Googling for:

[git delete unpushed commit](https://stackoverflow.com/a/3197432)

Says you can do a hard reset via:

```
# This is ok on this test repo but
# DANGER: things like "--force" and "--hard" are some of the few ways to permanently lose code in Git
git reset --hard HEAD~1
```

Note: As this is so dangerous I look it up on Stack overflow to make sure I got the command exactly right

## Emergency maneuvers

You could have also done it this way:

![](https://imgs.xkcd.com/comics/git_2x.png)
