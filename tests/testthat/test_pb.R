
set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

my_counts_bad <- my_counts
my_counts_bad[1:10] <- 0
test_data_bad <- data.frame("W" = my_counts_bad, "M" = seq_depth, my_covariate)

out <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             nstart = 1)
out_nullmu <- bbdml(formula = cbind(W, M - W) ~ 1,
                    phi.formula = ~ X1,
                    data = test_data,
                    link = "logit",
                    phi.link = "logit",
                    nstart = 1)
out_nullphi <- bbdml(formula = cbind(W, M - W) ~ X1,
                     phi.formula = ~ 1,
                     data = test_data,
                     link = "logit",
                     phi.link = "logit",
                     nstart = 1)

out_bad <- bbdml(formula = cbind(W, M - W) ~ X1,
                 phi.formula = ~ X1,
                 data = test_data_bad,
                 link = "logit",
                 phi.link = "logit",
                 nstart = 1)


test_that("pbWald works", {
  expect_is(pbWald(out, out_nullmu, B = 10), "numeric")
  expect_true(is.na(pbWald(out_bad, out_nullmu, B = 10)))
})

test_that("pbLRT works", {
  expect_is(pbLRT(out, out_nullmu, B = 10), "numeric")
})


out_error <- out
out_error$phi.link <- "break"

test_that("doBoot breaks properly", {
  expect_true(is.na(doBoot(out_error, out_nullmu, test = "LRT")))
  expect_true(is.na(doBoot(out, out, test = "Wald")))
})
