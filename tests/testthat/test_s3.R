library(corncob)
context("Test S3")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out <- bbdml(formula = cbind(W, M - W) ~ X1,
        phi.formula = ~ X1,
        data = test_data,
        link = "logit",
        phi.link = "logit",
        nstart = 1)

out_fish <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "fishZ",
             nstart = 1)

test_that("bbdml is S3", {
  expect_is(out, "bbdml")
  expect_is(out_fish, "bbdml")
})

test_that("bbdml S3 plotting works", {
  expect_is(plot(out, facet = "X1"), "ggplot")
  expect_is(plot(out, color = "X1", shape = "X1", facet = "X1", AA = TRUE), "ggplot")
  expect_is(plot(out, color = as.character(1:20), shape = 1:20), "ggplot")
  expect_error(plot(out, color = c(1,2)))
  expect_error(plot(out, shape = c(1,2)))
})

test_that("bbdml print generic works", {
  expect_is(print(out), "NULL")
})

test_that("summary function works", {
  expect_is(summary(out), "summary.bbdml")
})

test_that("summary printing works", {
  expect_null(print(summary(out)))
})
