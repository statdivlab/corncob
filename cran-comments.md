## Test environments
* local OS X install, R 4.0.3
* ubuntu 20.04 (via GitHub Actions)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Bryan D Martin <bmartin6@uw.edu>'

New submission

This is a new submission with myself as a maintainer.

## Version 0.3.0

* This version fixes a dependency issue with differentialTest
* Fixes all notes with if() conditions for comparing class() to string
* Fixes dead link in README
* Remove codecov link for CRAN invalid URL issue
* Other minor revisions to fix notes and improve performance

## Version 0.3.1

* Add --compact-vignettes=both to build options to reduce PDF size for WARN
