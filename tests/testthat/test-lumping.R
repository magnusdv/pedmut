
test_that("lumping works for custum matrix", {
  m1 = mutationMatrix("custom", matrix = matrix(1/3, ncol=3, nrow=3), alleles=1:3)
  lumped1 = lumpedMatrix(m1, lump = 2:3)
  expect_equal(dim(lumped1), c(2,2))
  expect_equal(colnames(lumped1), c("1", "lump"))
  expect_equal(attr(lumped1, 'lumpedAlleles'), c("2", "3"))
})

test_that("wrong input is caught when lumping matrix", {
  mat = mutationMatrix("eq", alleles = 1:4, rate = 0.1)
  expect_error(lumpedMatrix(mat, lump = 2:3, afreq = 1),
               "Frequency vector does not match the number of alleles")
  expect_error(lumpedMatrix(mat, lump = 2:3, afreq = 1:4),
               "Allele frequencies do not sum to 1")

  mod = mutationModel("eq", alleles = 1:4, rate = 0.1)

  expect_error(lumpedModel(mod, lump = 2:3, afreq = 1),
               "Frequency vector does not match the number of alleles")
  expect_error(lumpedModel(mod, lump = 2:3, afreq = 1:4),
               "Allele frequencies do not sum to 1")
})
