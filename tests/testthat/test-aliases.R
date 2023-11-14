# ======================================================== #
# rtrunc methods vs aliases                                #
# ======================================================== #

context("Matching output of rtrunc aliases")

test_that("rtrunc works the same from generic and alias", {
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 32, 22, .4, .6, family = "beta", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncbeta(1000, 32, 22, .4, .6)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 4, .2, 2, family = "binomial", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncbinom(1000, 4, .2, 2)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 45, 40, family = "chisq", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncchisq(1000, 45, 40)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, .2, .2, .7, family = "contbern", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtrunccontbern(1000, .2, .2, .7)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 0.584, 0.5, 1, family = "exp", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncexp(1000, 0.584, 0.5, 1)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, shape = 45, rate = 12, a = 1, b = 8, family = "gamma", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncgamma(1000, shape = 45, rate = 12, a = 1, b = 8)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, shape = 45, rate = 12, a = 0.2, b = 0.3, family = "invgamma", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncinvgamma(1000, shape = 45, rate = 12, a = 0.2, b = 0.3)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 38, 38, 0, 100, family = "invgauss", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncinvgauss(1000, 38, 38, 0, 100)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, -32, 2, family = "lognormal", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtrunclnorm(1000, -32, 2)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 152, .37, family = "nbinom", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncnbinom(1000, 152, .37)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 100, 200, 0, family = "normal", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncnorm(1000, 100, 200, 0)
    },
  )
  expect_identical(
    object = {
      set.seed(8)
      rtrunc(1000, 7521, 7500, family = "poisson", faster = FALSE)
    },
    expected = {
      set.seed(8)
      rtruncpois(1000, 7521, 7500)
    },
  )
})

# ======================================================== #
# rtrunc functions vs stats functions                      #
# ======================================================== #

context("Matching output of stats::r*")

test_that("Output of rtrunc matches stats::r*", {
  expect_setequal(
    object = {
      set.seed(1)
      rbeta(500, shape1 = 8, shape2 = 86)
    },
    expected = {
      set.seed(1)
      rtrunc(500, shape1 = 8, shape2 = 86, family = "beta")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rbinom(50, size = 10, prob = .3)
    },
    expected = {
      set.seed(1)
      rtrunc(50, size = 10, prob = .3, family = "binomial")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rchisq(50, df = 23)
    },
    expected = {
      set.seed(1)
      rtrunc(50, df = 23, family = "chisq")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rcontbern(50, lambda = 0.5)
    },
    expected = {
      set.seed(1)
      rtrunc(50, lambda = 0.5, family = "contbern")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rexp(50, rate = 26)
    },
    expected = {
      set.seed(1)
      rtrunc(50, rate = 26, family = "exp")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rgamma(50, shape = 4, rate = 5)
    },
    expected = {
      set.seed(1)
      rtrunc(50, shape = 4, rate = 5, family = "gamma")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rgamma(50, shape = 4, scale = 6)
    },
    expected = {
      set.seed(1)
      rtrunc(50, shape = 4, scale = 6, family = "gamma")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rinvgamma(50, shape = 6, scale = 9)
    },
    expected = {
      set.seed(1)
      rtrunc(50, shape = 6, scale = 9, family = "invgamma")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rinvgamma(50, shape = 52, scale = .21)
    },
    expected = {
      set.seed(1)
      rtrunc(50, shape = 52, scale = .21, family = "invgamma")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rinvgauss(500, m = 1, s = 3)
    },
    expected = {
      set.seed(1)
      rtrunc(500, m = 1, s = 3, family = "invgauss")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rlnorm(50, meanlog = 7, sdlog = 2)
    },
    expected = {
      set.seed(1)
      rtrunc(50, meanlog = 7, sdlog = 2, family = "lognormal")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rnbinom(500, size = 55, prob = .4)
    },
    expected = {
      set.seed(1)
      rtrunc(500, size = 55, prob = .4, family = "nbinom")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rnbinom(500, size = 55, mu = 4)
    },
    expected = {
      set.seed(1)
      rtrunc(500, size = 55, mu = 4, family = "nbinom")
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rnorm(500, mean = 1, sd = 3)
    },
    expected = {
      set.seed(1)
      rtrunc(500, mean = 1, sd = 3)
    }
  )
  expect_setequal(
    object = {
      set.seed(1)
      rpois(500, lambda = 72)
    },
    expected = {
      set.seed(1)
      rtrunc(500, lambda = 72, family = "poisson")
    }
  )
})
