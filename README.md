# MonkeyPipeline

## Live version  
Sean Thomas and Alex WIlliams 2013-2016  
**version 2** - parallelized for Rigel / works with TORQUE PBS scheudler (should work with other PBS schedulers). Can operate without a scheduler if you add "--notorque" to the command line.

----------

### Basic Rules:

- Don't implement changes directly on the `master` branch.
- Implement fixes and features in a feature branch in your personal/local copy of the repo.
    - Push feature branches to a counterpart branch in remote/origin (i.e. GitHub) with the same name.
    - When done implementing fix or feature open a Pull Request (PR) to merge the changes into `master`.
    - Ideally have someone else review your code and approve those changes, i.e. approve the PR.
    - **feature branch naming convention:** *initials*-feature-branch-name, e.g. df-add_new_feature.
- Generally use `git pull --rebase` when pulling changes (say from `master`) from GitHub, and `git merge` when pushing changes up the hierarchy/GitHub. See these:  
    - https://www.derekgourlay.com/blog/git-when-to-merge-vs-when-to-rebase/  
    - https://www.atlassian.com/git/tutorials/merging-vs-rebasing
- **TYPICAL ORDER OF STEPS** (when you have local changes and you want to come up-to-date with remote)**:** `commit` local changes, `pull` from remote, `push` to remote.
