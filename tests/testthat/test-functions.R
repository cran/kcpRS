test_that("kcpRS finds change point", {
  set.seed(15)
  S=100
  X1=cbind(rnorm(S,0,1),rnorm(S,0,1)) #phase1: Means=0
  X2=cbind(rnorm(S,1,1),rnorm(S,1,1)) #phase2: Means=1
  X=rbind(X1,X2)
  X=scale(X)
  
  result = kcpRS(X, RS_fun = runMean, RS_name = "Mean")

  expect_equal(result$changePoints, 101, tolerance = 5)
  expect_output(str(result), "List of 16")
  expect_s3_class(result, "kcpRS")
  expect_output(summary(result))
  expect_output(print(result))
})

test_that("kcpRS finds change point in workflow", {
  set.seed(15)
  S=100
  X1=cbind(rnorm(S,0,1),rnorm(S,0,1)) #phase1: Means=0
  X2=cbind(rnorm(S,1,1),rnorm(S,1,1)) #phase2: Means=1
  X=rbind(X1,X2)
  X=scale(X)
  
  result = kcpRS_workflow(X)
  
  expect_equal(result$kcpMean$changePoints, 101, tolerance = 5)
  expect_output(str(result), "List of 16")
  expect_s3_class(result, "kcpRS_workflow")
  expect_output(summary(result))
  expect_output(print(result))
})
