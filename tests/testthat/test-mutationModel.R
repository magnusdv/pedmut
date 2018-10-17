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

test_that("stepwise mutation model splits correctly in male/female", {
  alleles = c(1, 1.5, 2,3)
  M = mutationModel("step", alleles = alleles,
                    rate = list(male=0, female=0.1),
                    rate2 = list(male=0.2, female=0),
                    range = list(male=2, female=3))
  matM = mutationMatrix("step", alleles = alleles,
                        rate = 0, rate2 = 0.2, range = 2)
  matF = mutationMatrix("step", alleles = alleles,
                         rate = 0.1, rate2 = 0, range = 3)
  expect_identical(M, mutationModel(list(male=matM, female=matF)))
})
