context("test-general")

#RUN CODE
output_string <- rMSIcleanup::hello()
t=0
#elapsed_t=system.time(for(i in 1:10000) x <- mean(rt(1000, df = 4)))
#t=as.numeric(elapsed_t)[1]

#NUMERIC
a=2*2
test_that("equal", {
  expect_equal(a, 4)
})
test_that("less", {
  expect_lt(a, 5)
})
test_that("less or equal", {
  expect_lte(a, 4)
})
test_that("greater", {
  expect_gt(a, 3)
})
test_that("greater or equal", {
  expect_gte(a, 4)
})

test_that("multiplication works 2", {
  expect_equal(a, a)
})

#CHECK OUTPUT STRING
test_that("Not empty", {
  expect_match(output_string, ".")
})

test_that("Includes world", {
  expect_match(output_string, "world")
})

#CHECK RUNTIME
test_that("Fast code", {
  expect_lt(t, 3)
})
