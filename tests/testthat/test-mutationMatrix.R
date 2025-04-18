
test_that("bad custom matrices are caught", {
  expect_error(mutationMatrix(),
               "`matrix` cannot be NULL with the `custom` model")
  expect_error(mutationMatrix(matrix=data.frame(a=1)),
               "Custom matrix must be a matrix, not a data.frame")
  expect_error(mutationMatrix(matrix=list(a=1)),
               "Custom matrix must be a matrix, not a list")
  expect_error(mutationMatrix(matrix=1),
               "Custom matrix must be a matrix, not a numeric")
  expect_error(mutationMatrix(matrix=matrix("a", ncol=1)),
               "Custom matrix must be numeric, not character")
  expect_error(mutationMatrix(matrix=matrix(0, nrow=2, ncol=1)),
               "Custom matrix must be square, not 2*1", fixed=T)
  expect_error(mutationMatrix(matrix=matrix(1, ncol=1), alleles=1:2),
               "Length of `alleles` must equal the dimension of `matrix`")
  expect_error(mutationMatrix(matrix=matrix(1, ncol=1)),
               "When custom matrix lacks names, the argument `alleles` cannot be NULL")
  m = matrix(1/2, 2, 2, dimnames = list(c("a", "b"), NULL))
  expect_error(mutationMatrix(matrix = m),
               "Mutation matrix has unequal row and column names")
})

test_that("bad inputs are caught", {
  expect_error(mutationMatrix(model="eq", matrix=matrix(1, ncol=1)),
               "`matrix` cannot be used with the `equal` model")
  expect_error(mutationMatrix(model="eq"),
               "`alleles` cannot be NULL with the `equal` model")
  expect_error(mutationMatrix(model="eq", alleles=1:2),
               "`rate` cannot be NULL with the `equal` model")
  expect_error(mutationMatrix(model="prop", alleles=1:2, rate=0),
               "`afreq` cannot be NULL with the `proportional` model")
  expect_error(mutationMatrix(model="prop", alleles=1:2, afreq=1),
               "Frequency vector does not match the number of alleles")
  expect_error(mutationMatrix(model="prop", alleles=1:2, afreq=c(0.5, 0.501)),
               "Allele frequencies do not sum to 1")
  expect_error(mutationMatrix(model="prop", alleles=1:2, afreq=c(.5,.5), rate=2),
               "`rate` must be a number in the interval")
  expect_error(mutationMatrix(model="prop", alleles=1:2, afreq=c(.3,.7), rate=1),
               "Model undefined: max `rate` for the given input is: 0.6")
  expect_error(mutationMatrix(model="step", alleles=1:2, rate=1, rate2=0.5, range=0.1),
               "The total mutation rate")
  expect_error(mutationMatrix(model="step", alleles=1:2, rate=0, rate2=0),
               "`range` cannot be NULL with the `stepwise` model")
  expect_error(mutationMatrix(model="dawid", alleles=1:2, afreq = c(.3,.7), rate=0, rate2=0, range=1),
               "For the `dawid` model, `range` must be in the interval")
})


test_that("equal model works", {
  expect_equivalent(
    mutationMatrix(model = "eq", alleles = 1:3, rate = 0),
    mutationMatrix(matrix = diag(3), alleles = 1:3))

  m = mutationMatrix(alleles = 1:3, model = "equal", rate = 0.1)
  expect_silent(validateMutationMatrix(m, alleles = 1:3))
  expect_equivalent(diag(m), rep(0.9, 3))
  expect_equivalent(rowSums(m), rep(1,3))

  expect_equivalent(m, mutationMatrix(matrix = m))
})

test_that("proportional model works", {
  afr = c(0.2, 0.3, 0.5)
  expect_equivalent(
    mutationMatrix(model = "prop", alleles = 1:3, rate = 0, afreq=afr),
    mutationMatrix(matrix = diag(3), alleles = 1:3))

  m = mutationMatrix(model = "prop", alleles = 1:3, rate = 0.1, afreq=afr)
  expect_silent(validateMutationMatrix(m, alleles = 1:3))
  expect_equivalent(rowSums(m), rep(1,3))

  expect_equal(c(m[2,1], m[1,2], m[1,3]), c(m[3,1], m[3,2], m[2,3]))

  expect_equivalent(m, mutationMatrix(matrix = m))
})

test_that("random model works", {
  m = mutationMatrix(alleles = 1:3, model = "random")
  expect_silent(validateMutationMatrix(m, alleles = 1:3))
  expect_equivalent(rowSums(m), rep(1,3))

  expect_equivalent(m, mutationMatrix(matrix = m))
})

test_that("trivial model works", {
  m = mutationMatrix(alleles = 1:3, model = "triv")
  expect_silent(validateMutationMatrix(m, alleles = 1:3))
  expect_equivalent(m, mutationMatrix(matrix = diag(3), alleles = 1:3))

  expect_equivalent(m, mutationMatrix(matrix = m))
})

test_that("trivial stepwise model is diagonal", {
  m = mutationMatrix(alleles = 1:3, model = "step",
                     rate=0, rate2=0, range=0.5)

  expect_silent(validateMutationMatrix(m, alleles = 1:3))
  expect_equivalent(m, mutationMatrix(matrix = diag(3), alleles = 1:3))
})

test_that("non-trivial stepwise mutation matrix is correct", {
  alleles = c(1, 1.5, 2,3)
  rate = 0.6
  rate2 = 0.3
  range = 0.5
  M = mutationMatrix("step", alleles=alleles, rate=rate, rate2=rate2, range = range)
  expect_equivalent(unclass(M),
                    matrix(c(.1,.1,.3,.2,.3,.7,.3,.3,.4,.1,.1,.4,.2,.1,.3,.1),
                           ncol=4, dimnames = list(alleles, alleles)))
})
