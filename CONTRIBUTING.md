# BifrostTools.jl contributions

BifrostTools.jl includes the most basic functions to read Bifrost simulation input, and do interpolations and derivatives on the staggered mesh. If you want more functionality, feel free to open a pull request or raise an issue in the repository.

### Contributing code

If you want to contribute code, fork the repository and open a pull request. It best to open a PR directly into the develop branch.

We are not strict with coding style, but please document new functions you add with docstrings. It is also a huge bonus if you can add testing for new functions, but this is not necessary.

### Bugs

If you find bugs in the code, you have three options:

1. raise an issue in the repository
2. message the developers (Elias and Eilif)
3. fix it yourself in a pull request (if you are familiar with Julia)

### Git commit message manners
The commit message is mainly for other people, so they should be able to understand it now and six months later. Commit messages cannot be longer than one sentence (line) and should start with a tag identifier (see the end of this section).

Use the imperative form of verbs rather than past tense when referring to changes introduced by the commit in question. For example, "Remove property X", not "Removed property X" and not "I removed...". This tense makes picking, reviewing or reverting commits more readable.

Use following tags for commit messages:

       [DEV] : Code development (including additions and deletions)
       [ADD] : Adding new feature
       [DEL] : Removing files, routines
       [FIX] : Fixes that occur during development, but which have essentially no impact on previous work
       [BUG] : Bug with significant impact on previous work -- `grep`-ing should give restricted list
       [OPT] : Optimisation
       [DBG] : Debugging
       [ORG] : Organisational, no changes to functionality
       [SYN] : Typos and misspellings (including simple syntax error fixes)
       [DOC] : Documentation only
       [REP] : Repository related changes (e.g., changes in the ignore list, remove files)
       [UTL] : Changes in utils

Commit message examples:

* "[BUG] Add missing initialisation to tg array"
* "[FIX] Add lowercase castig to params"
* "[CLN] Remove unnecessary allocation for dynamic arrays."
