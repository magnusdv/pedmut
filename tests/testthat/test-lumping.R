
test_that("lumping works for custum matrix", {
  m1 = mutationMatrix("custom", matrix = matrix(1/3, ncol=3, nrow=3), alleles=1:3)
  lumped1 = lumpedMatrix(m1, lump = 2:3)
  expect_equal(dim(lumped1), c(2,2))
  expect_equal(colnames(lumped1), c("1", "lump"))
  expect_equal(attr(lumped1, 'lumpedAlleles'), 2:3)
})
