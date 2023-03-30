test_that("stabilization works", {
  m = mutationMatrix("random", alleles = 1:2, afreq = c(.5,.5), seed = 123)
  m2 = stabilize(m)
  expect_true(isStationary(m2))
})
