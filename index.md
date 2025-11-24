## Automation and visualization of flow cytometry data analysis pipelines

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/CytoPipelineUtils/coverage.svg?branch=main)](https://codecov.io/github/UCLouvain-CBIO/CytoPipelineUtils?branch=main)
[![license](https://img.shields.io/badge/license-GPL3.0-blue)](https://opensource.org/licenses/GPL-3.0)

### What is CytoPipelineUtils?

`CytoPipelineUtils` is meant to be used in conjunction with
`CytoPipeline` package (see
[here](https://github.com/UCLouvain-CBIO/CytoPipeline)), which provides
support for automation and visualization of flow cytometry data analysis
pipelines. `CytoPipelineUtils` provides a series of wrapper
implementations that can in turn be defined as CytoProcessingSteps in
CytoPipeline. It is therefore a helper package, which is aimed at
hosting wrapper implementations of various functions of various
packages.

`CytoPipelineUtils` is open to contributions. If you want to implement
your own wrapper of your favourite pre-processing function and use it in
a `CytoPipeline` object, this is the place to do it! If you do so,
please issue a [pull
request](https://github.com/UCLouvain-CBIO/CytoPipelineUtils/pulls).

### License

The `CytoPipelineUtils` code is provided under [GPL license version 3.0
or higher](https://opensource.org/licenses/GPL-3.0). The documentation,
including the manual pages and the vignettes, are distributed under a
[CC BY-SA 4.0 license](https://creativecommons.org/licenses/by-sa/4.0/).
