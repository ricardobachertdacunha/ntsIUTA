
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntsIUTA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The objective of `ntsIUTA` is to simplify and facilitate application of
highly efficient open source tools for processing of chromatographically
separated mass spectral (MS) data, using simple
[R](https://cran.r-project.org/bin/windows/base/) objects and
comprehensive visualization for both optimization of data extraction and
production of results. `ntsIUTA` can be seen as an oriented guide for
NTS with a generic but flixible and application-driven functionality.
The package `ntsIUTA` relies heavily on open source tools for non-target
screening (NTS) using MS data. It does not invent the wheel, but eases
its usability. In particular, the packages `patRoon` and `Bioconductor`
play a fundamental role on the functionality of `ntsIUTA`. While
`MSnbase` and `xcms` from the `Bioconductor` are used to extract and
visualize raw spectral data, `patRoon` is crucial for downstream
integration of external software. Although, simple objects, such as
[list](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/list),
[data.frame](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/data.frame)
and
[vector](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/vector),
are used as the main output, native objects from `patRoon` and
`Bioconductor` are preserved to facilitate the use of respective native
functions.

## Installation

To add.
<!--You can install the released version of ntsIUTA from [CRAN](https://CRAN.R-project.org) with:-->
<!-- ``` r --> <!-- install.packages("ntsIUTA") --> <!-- ``` -->

## Information

Loading `ntsIUTA`:

``` r
library(ntsIUTA)
```

See vignettes:

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
