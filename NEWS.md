# pedmut 0.9.0

## New features

* New function `lumpMutSpecial()` implementing special lumping of mutation models. Only some cases are implemented so far, more may be added in the future.

* `mutationMatrix()` gains argument `validate`.

* New function `makeStationary()` replaces the previous `stabilize()` for PM-transformation to stationarity.

* As a shortcut for `makeStationary()`, the PM transformation may also be applied directly with `mutationMatrix(..., transform = "PM")` or `mutationModel(..., transform = "PM")`.

## Bug fixes

* Fixed a bug in `makeReversible(..., method = "PR")` failing to catch non-reversible inputs.

* Various minor fixes and code simplifications.


# pedmut 0.8.0

## New features

* New function `makeReversible()` implementing three transformations to reversibility.

* `mutationMatrix()` and `mutationModel()` gain a `transform` argument for applying transformations to reversibility on the fly.

* New function `adjustRate()` for adjusting the overall mutation rate of a mutation matrix.

* `mutationModel(model = "dawid", ...)` replaces the function `stepwiseReversible()`.

* Random mutation models can now condition on a fixed overall mutation rate.

* Include info on boundedness in the print method for mutation models.

## Other

* Refactored `mutationMatrix()` for better code maintainability.

* Bug fix: Don't trivialise model name when there is only 1 allele.

* Several minor improvements in documentation and examples.

* Use rhub v2.


# pedmut 0.7.0

## New features

* New function `mutRate()` computing the overall mutation rate of a mutation model.

* New function `isBounded()` checking if a mutation matrix M is bounded by the allele frequencies, i.e., that `M[i,j] <= afreq[j]` for all `i != j`.

* New function `stepwiseReversible()`.

* Added Thore Egeland as package contributor.

## Bug fixes

* Fixed a floating point bug in the creation of stepwise models.

* Fix broken URLs.


# pedmut 0.6.0

## New features

* Implement the PM stabilisation method of Simonsson & Mostad; liftover from `Familias:::stabilize`.

* Implement more general lumping (multiple lumps!) under strong lumpability. 

* New utility `getParams()` for extracting model parameters.

* Extend `isStationary()`, `isReversible()` and `alwaysLumpable()` to full models.

* New functions `isMutationMatrix()` and `isMutationModel()`.

* New S3 method `as.matrix()` for classes `mutationMatrix` and `mutationModel`.

* New function `findStationary()` for obtaining the stationary distribution of a mutation matrix.


# pedmut 0.5.0

## New features

*  Various speedups, mainly concerning lumped models. These improvements should give better performance in other parts of the `pedsuite`, especially in likelihood calculations with many markers.

* New function `lumpedModel()` which is a convenient wrapper of `lumpedMatrix()`.

* New function `sexEqual()` for reading the `sexEqual` attribute of a mutation model.



# pedmut 0.4.0

## New features

* `mutationModel()` gains a new parameter `validate`.
* Various code improvements in `validateMutationModel()` and `validateMutationMatrix()`, making them significantly faster.

# pedmut 0.3.0

This is maintenance release with mostly minor changes.

# pedmut 0.2.0

## New features
* Added `onestep` model.
* Added `toString()` method producing a short description of a mutation model.
This is used in `pedtools::print.marker()`.

## Other changes

* Updated README.
* Minor bug fixes.


# pedmut 0.1.0

* Initial CRAN release
