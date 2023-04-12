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

*  Various speedups, mainly concerning lumped models. These improvements should give better performance in other parts of the `ped suite`, especially in likelihood calculations with many markers.

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
