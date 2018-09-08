context("properties of mutation matrices")

test_that("equal model is station/revers only for equal freqs", {
  mutmat = mutationMatrix(alleles = 1:3, model = "equal", rate = 0.1)
  expect_false(isStationary(mutmat, c(0.2, 0.3, 0.5)))
  expect_false(isReversible(mutmat, c(0.2, 0.3, 0.5)))

  expect_true(isStationary(mutmat, c(1,1,1)/3))
  expect_true(isReversible(mutmat, c(1,1,1)/3))
})

test_that("proportial model is stationary and reversible", {
  afreq = c(0.2, 0.3, 0.5)
  mutmat = mutationMatrix(alleles = 1:3, model = "prop", rate = 0.1, afreq = afreq)
  expect_true(isStationary(mutmat, afreq))
  expect_true(isReversible(mutmat, afreq))
})

test_that("random model 3*3 is not stationary (nor reversible)", {
  afreq = c(0.2, 0.3, 0.5)
  mutmat = mutationMatrix(alleles = 1:3, model = "random", rate = 0.1)
  expect_false(isStationary(mutmat, afreq))
  expect_false(isReversible(mutmat, afreq))
})

test_that("equal model is always lumpable", {
  mutmat = mutationMatrix(alleles = 1:4, model = "equal", rate = 0.1)
  expect_true(isLumpable(mutmat, lump = 1:2))
  expect_true(isLumpable(mutmat, lump = c(2,4)))
  expect_true(isLumpable(mutmat, lump = c(1:2,4)))
})

test_that("proportional model is always lumpable", {
  afreq = c(0.2, 0.3, 0.5)
  mutmat = mutationMatrix(alleles = 1:3, model = "prop", rate = 0.1, afreq = afreq)
  expect_true(isLumpable(mutmat, lump = 1:2))
  expect_true(isLumpable(mutmat, lump = c(1,3)))
})
