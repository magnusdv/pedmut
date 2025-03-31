
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

test_that("random model 3*3 has correct rate", {
  afreq = c(a=0.2, b=0.3, c=0.5)
  m1 = mutationMatrix(afreq = afreq, model = "random", rate = 0.001, seed = 1)
  expect_equal(mutRate(m1), 0.001)
  m2 = mutationMatrix(afreq = afreq, model = "random", rate = 0.999, seed = 1)
  expect_equal(mutRate(m2), 0.999)
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

test_that("isStationary and isReversible works with full models", {
  afreq = c(0.5, 0.5)
  m1 = mutationMatrix("eq", alleles = 1:2, rate = 0.1, afreq = afreq)
  m2 = mutationMatrix("random", alleles = 1:2, afreq = afreq)
  expect_true(isStationary(mutationModel(list(female = m1, male = m1))))
  expect_false(isStationary(mutationModel(list(female = m1, male = m2))))
  expect_false(isStationary(mutationModel(list(female = m2, male = m2))))

  expect_true(isReversible(mutationModel(list(female = m1, male = m1))))
  expect_false(isReversible(mutationModel(list(female = m1, male = m2))))
  expect_false(isReversible(mutationModel(list(female = m2, male = m2))))

  stb = stabilize(mutationModel(list(female = m1, male = m2)))
  expect_true(isStationary(stb))
})
