context("Generated examples")
test_that("getWindEstimates",{
data("storks")
  expect_equivalent(class(tmp<-getWindEstimates(storks[[1:2]])),"MoveStack")
})