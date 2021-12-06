
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntsIUTA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The objective of `ntsIUTA` is to facilitate the development and routine
use of application driven workflows for non-target screening, using open
source tools. The packages
[patRoon](https://github.com/rickhelmus/patRoon),
[MSnbase](https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html)
and
[xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html)
play a fundamental role to `ntsIUTA`. Therefore, installation and
experience with these packages is recommended.

## Installation

`ntsIUTA` depends on external open source software. Therefore,
installation of these dependencies should occur beforehand as described
in following Prerequisites section.

## Prerequisites

Confirm that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is
installed and the PATH is set up correctly in the environmental
variables. To set the PATH the following line can be used where PATH is
the location of the RTools bin folder.

``` r
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
```

`ntsIUTA` largely depends on
[patRoon](https://github.com/rickhelmus/patRoon) package and both have a
large number of common dependencies. Therefore, installation of
[patRoon](https://github.com/rickhelmus/patRoon) is recommended before
installing `ntsIUTA`. An installation guide for the
[patRoon](https://github.com/rickhelmus/patRoon) package can be found
[here](https://rickhelmus.github.io/patRoon/handbook_bd/installation.html).
[patRoon](https://github.com/rickhelmus/patRoon) also provides an
[automatic installation
script](https://rickhelmus.github.io/patRoon/handbook_bd/automatic-installation-windows-only.html)
for Windows users. Alternatively, a [manual installation
guide](https://rickhelmus.github.io/patRoon/handbook_bd/manual-installation.html)
is also provided. The most basic steps of the
[patRoon](https://github.com/rickhelmus/patRoon) installation are as
follows.

``` r
install.packages("remotes")
remotes::install_github("rickhelmus/patRoon", upgrade = "never")

#Direct installation of packages from Bioconductor
install.packages("BiocManager")
BiocManager::install(c("CAMERA", "mzR", "xcms", "Biobase", "MSnbase"))
```

The code above installs a basic configuration of
[patRoon](https://github.com/rickhelmus/patRoon). External software,
such as
[proteowizard](http://proteowizard.sourceforge.net/download.html),
[MetFrag](https://ipb-halle.github.io/MetFrag/projects/metfragcl/),
needs to be installed and referenced in the `.Rprofile` manually as
suggested
[here](https://rickhelmus.github.io/patRoon/handbook_bd/installation.html).
Alternatively, the automatic installation of
[patRoon](https://github.com/rickhelmus/patRoon) can be used which
creates the `.Rprofile-patRoon.R` with all required locations of the
external software used in
[patRoon](https://github.com/rickhelmus/patRoon) and consequently
`ntsIUTA`.

## Installation of the ntsIUTA package

After installing the prerequisites, the `ntsIUTA`package can then be
installed from the [ntsIUTA GitHub
repository](https://github.com/ricardobachertdacunha/ntsIUTA) using the
code line below.

``` r
remotes::install_github("ricardobachertdacunha/ntsIUTA")


#If the installation fails with an error stating patRoon was not installed for i386, please use this command instead:
remotes::install_github("ricardobachertdacunha/ntsIUTA", INSTALL_opts = "--no-multiarch")

```

<!-- ## Information -->
<!-- Loading `ntsIUTA`: -->
<!-- See vignettes: -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
