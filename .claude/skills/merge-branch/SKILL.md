---
name: merge-branch
description: Land the current feature branch the way this repository does - a --no-ff merge commit into master, delete the branch, push master to origin. Use only when the user asks to merge a branch and/or push it.
---

# Merge a branch into master and push

Stop and report, rather than work around it, whenever a check below fails.

1. **Preconditions.** The current branch is not `master`, and `git status --short` is empty.
   Do not stash to make it empty (and never `git stash --all` or `git clean -x`: they reach the
   gitignored `data/` - see CLAUDE.md).
2. **Origin has not moved.** `git fetch -q origin`, then
   `git rev-list --left-right --count origin/master...master` must start with `0`. If origin is
   ahead, stop and ask; never rebase, reset or force-push to get past it.
3. **What lands.** `git log --oneline master..<branch>` - summarise it in one line.
4. **Merge.** `git checkout master`, then
   `git merge --no-ff <branch> -m "Merge <branch>: <one-line summary>"` with a blank line and the
   commit attribution trailer this session's instructions give. Match the style of
   `git log --merges --format=%s -3`.
5. **Delete the branch.** `git branch -d <branch>` (lowercase `-d`: it refuses an unmerged branch).
6. **Push.** `git push origin master`. Never `--force`.
7. **Report** the merge commit, the pushed range (`old..new`), and that `git status -sb` shows
   `master...origin/master` with nothing ahead or behind.
