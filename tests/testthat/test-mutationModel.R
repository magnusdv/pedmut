context("creation of mutation models")

test_that("mutation models are validated", {
  alleles = 1:3
  afr = c(.2,.3,.5)
  rate = 0.1
  M = mutationModel("prop", alleles=alleles, afreq=afr, rate=rate)
  expect_silent(validateMutationModel(M))
})


test_that("male and female models are equal by default", {
  alleles = 1:3
  afr = c(.2,.3,.5)
  rate = 0.1
  m = mutationMatrix("prop", alleles=alleles, afreq=afr, rate=rate)

  M = mutationModel("prop", alleles=alleles, afreq=afr, rate=rate)
  expect_identical(M$female, m)
  expect_identical(M$male, m)
})
