# Multivariate Longitudinal Factor Analysis

The *mlfa* R package implements a three-stage multivariate longitudinal factor analysis model, as well as helper functions for data perprocessing and plotting functionality. The model, which is detailed in the [accompanying manuscript](TODO), allows to characterise latent longitudinal trajectories and relate them linearly with the measurements underlying the data. On a high-level, the model consists of a sparsified Probabilistic Principal Component Analysis based dimensionality reduction (stages one and two) followed by a longitudinal smoothing with Generalized Additive Mixed models (stage three).

## Installation

You can install *mlfa* from CRAN: 

```
install.packages("mlfa")
```

Alternatively, you can install the latest development version from Github: 

```
# install.packages("devtools")
devtools::install_github("FabianFalck/mlfa")
```

If you want to load the package from source without installing it (e.g. to develop the package), you can do:

```
# install.packages("devtools")
devtools::load_all("<path-to-mlfa-package-code>/mlfa")
```

## Overview

This package has 2 main functions: 
* `mfd`: This function performs stage 1 and 2 of the model. 
* `mlfa`: This function performs stage 3 of the model.
Each of these function calls several package functions and returns an S3 object of the respective type. An object from the `mfd` class is input to a subsequent `mlfa` function call.

`mfd` and `mlfa` are accompanied by several helper functions for plotting and other diagnostics. Both classes implement the generic functions `plot` and `summary`. In particular `plot.mfd` and `plot.mlfa` call several helper functions at once. This provides a convenient way to use the package, as its interface and key functionality boil down to four function calls: Two for creating `mfd` and `mlfa` objects, and two for plotting their results.

## How to use

TODO short excerpt from Vignette

## Supplementary ressources

This package provides the following material for further reference: 

* [Vignette](TODO): The package Vignette applies the package on simulated data and demonstrates the core functions and their usage. It is understood as an ideal starting point for someone using the package for the first time.
* [Paper](TODO): Our Paper explains in detail the methodology that this package implements. 
* [Reference Manual](TODO): The Reference Manual documents all package functions. 
* [UML class chart](TODO): A UML class chart of the two S3 classes and its used functions.

## Cite

If you find this R package useful, please cite our paper: 

TODO 