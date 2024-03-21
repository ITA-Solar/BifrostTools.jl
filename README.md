# BifrostTools.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ITA-Solar.github.io/BifrostTools.jl/dev/)
[![Build Status](https://github.com/ITA-Solar/BifrostTools.jl/actions/workflows/CI.yml/badge.svg?branch=develop)](https://github.com/ITA-Solar/BifrostTools.jl/actions/workflows/CI.yml?query=branch%3Adevelop)
[![Coverage](https://codecov.io/gh/ITA-Solar/BifrostTools.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/ITA-Solar/BifrostTools.jl)

Tools for working with Bifrost in Julia.

This Julia package is created for working with *Bifrost* data **faster**. 
Load single or multiple simulation snapshots, and analyse data efficiently. 

This package is an extension of `Bifrost.jl`, a script written by Mikolaj Szydlarski.

## Documentation
The documentation is currently available at https://ita-solar.github.io/BifrostTools.jl/build/

___
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
