library(corncob)
context("Test noninteger Ws work")

set.seed(2)
seq_depth <- rpois(20, lambda = 10000) + rnorm(20)
my_counts <- rbinom(20, size = round(seq_depth), prob = 0.001) * 10 + rnorm(20)
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out_integer <- bbdml(formula = round(cbind(W, M - W)) ~ X1,
                     phi.formula = ~ X1,
                     data = test_data,
                     link = "logit",
                     phi.link = "logit",
                     nstart = 1)


test_that("dbetabinom_cts works as intended", {

  ws <- c(149.102302643938, 80.6917805312822, 131.116260412284, 120.183214155073,
          180.445439532545, 109.74950177694, 121.064836349872, 119.542685306637,
          139.587723111623, 109.73892259761, 80.1519520682902, 129.146809463979,
          88.6776462197223, 89.8491110259799, 89.8176118808277, 80.9499143980439,
          50.242792420626, 79.8828661542037, 78.8673326216628, 50.6320147222528)
  ms <- c(9908.88040089543, 9903.33561720997, 9886.49479576711, 10012.1380527087,
          9928.88120797422, 9976.19768426235, 9922.93130728875, 10040.1967867826,
          10096.8862348637, 9961.5800916837, 9964.49781876104, 9768.26264545866,
          10110.7670988004, 10031.9962764662, 9944.51167228281, 10208.5243017156,
          9880.79791643753, 10081.0259974389, 9961.68937269777, 9753.04416089672)

  expect_silent(dbetabinom_cts(x=ws, size=ms, prob = .5, rho = .1, log = TRUE))
  expect_equal(dbetabinom_cts(x=round(ws), size=ms, prob = .5, rho = .1, log = TRUE),
               VGAM::dbetabinom(x=round(ws), size=ms, prob = .5, rho = .1, log = TRUE))

  expect_equal(VGAM::dbetabinom(x=round(ws), size=round(ms), prob = .5, rho = .1, log = TRUE),
               dbetabinom_cts(x=round(ws), size=round(ms), prob = .5, rho = .1, log = TRUE))

})

test_that("Noninteger W's throw errors if not forced", {

  expect_error({
    out_noninteger_error <- bbdml(formula = cbind(W, M - W) ~ X1,
                                  phi.formula = ~ X1,
                                  data = test_data,
                                  link = "logit",
                                  phi.link = "logit",
                                  nstart = 1)
  })
})

out_noninteger <- bbdml(formula = cbind(W, M - W) ~ X1,
                        phi.formula = ~ X1,
                        data = test_data,
                        link = "logit",
                        phi.link = "logit",
                        nstart = 1,
                        allow_noninteger=TRUE)

test_that("Noninteger W's work", {

  expect_is(summary(out_noninteger), "summary.bbdml")

  # Noninteger W's are close to integer W's
  expect_equal(out_noninteger$param, out_integer$param, tolerance=0.05)

  # Sandwich standard errors work for noninteger W's

  expect_is(sand_vcov(out_integer), "matrix")
  expect_is(sand_vcov(out_noninteger), "matrix")
})

sand_vcov(out_noninteger) %>% diag %>% sqrt
out_noninteger %>% summary %>% coef

