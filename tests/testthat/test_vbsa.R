
context("Function vbsa()")

# Example: Ishigami test function 
ishigami <- function(x, a=7, b=0.1) {
  sin(x[1]) + a*(sin(x[2]))^2 + b*(x[3]^4)*sin(x[1]) 
}
nparam <- 3
pars_lower <- rep(-pi, nparam)
pars_upper <- rep(pi, nparam)

# TESTS #

test_that("default output with default arguments is returned as expected", {
  res_default <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper)
  
  expect_type(res_default, "list")
  expect_length(res_default, 5)
  expect_named(res_default, c("N", "fn.counts", "Si", "St", "ranking"))
  expect_equivalent(res_default$N, 1000)
  expect_equivalent(res_default$fn.counts, 1000*(nparam+2))
  invisible(lapply(res_default[c("Si", "St", "ranking")], function(x) expect_length(x, n=nparam)))
  expect_type(res_default$Si, "double")
  expect_type(res_default$St, "double")
  expect_type(res_default$ranking, "integer")
  expect_equivalent(res_default$ranking, order(res_default$Si, decreasing = T))
})

test_that("extended default output with argument 'full.output' is returned as expected", {
  res_full <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper,
                   full.output = T)
  
  expect_type(res_full, "list")
  expect_length(res_full, 11)
  expect_named(res_full, c("N", "fn.counts", "Si", "St", "ranking",
                              "Matrix.A", "Matrix.B", "Matrix.Ab", "fn.A", "fn.B", "fn.Ab"))
  invisible(lapply(res_full[c("Matrix.A", "Matrix.B", "Matrix.Ab")], 
                   function(x) expect_is(x, class="matrix")))
  expect_true(all(dim(res_full$Matrix.A) == c(1000, nparam)))
  expect_true(all(dim(res_full$Matrix.B) == c(1000, nparam)))
  expect_true(all(dim(res_full$Matrix.Ab) == c(nparam*1000, nparam)))
  expect_length(res_full$fn.A, 1000)
  expect_length(res_full$fn.B, 1000)
  expect_true(all(dim(res_full$fn.Ab) == c(1000, nparam)))
})

test_that("output is reasonable, i.e. correctly calculated from parameters and returned", {
  res_full <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, full.output = T)
  res_check_A <- apply(res_full$Matrix.A, 1, ishigami)
  res_check_B <- apply(res_full$Matrix.B, 1, ishigami)
  res_check_Ab <- matrix(apply(res_full$Matrix.Ab, 1, ishigami), ncol=nparam, byrow = T)
  
  expect_equivalent(res_check_A, res_full$fn.A)
  expect_equivalent(res_check_B, res_full$fn.B)
  expect_equivalent(res_check_Ab, res_full$fn.Ab)
})

test_that("bootstrapping works (returns expected output)", {
  res_boot <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper,
                   Nb = 500, Nb.sig = 0.05, full.output = T)
  
  expect_type(res_boot, "list")
  expect_length(res_boot, 21)
  expect_true(all(c("Si.sd", "Si.lo", "Si.up", "St.sd", "St.lo", "St.up",
                    "ranking.boot", "Si.boot", "St.boot", "B") %in% names(res_boot)))
  invisible(lapply(res_boot[c("Si.sd", "Si.lo", "Si.up", "St.sd", "St.lo", "St.up")], function(x) expect_length(x, n=nparam)))
  invisible(lapply(res_boot[c("ranking.boot", "Si.boot", "St.boot")], function(x) expect_true(all(dim(x) == c(500, nparam)))))
  expect_true(all(dim(res_boot$B[[1]]) == c(1000, 500)))
})

test_that("subsampling works (returns expected output)", {
  res_sub <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper,
                   subsamp = 100)
  
  expect_type(res_sub, "list")
  expect_length(res_sub, 5)
  
  invisible(lapply(res_sub[c("Si", "St", "ranking")], function(x) expect_equivalent(dim(x), c(10, nparam))))
  invisible(lapply(res_sub[c("Si", "St", "ranking")], function(x) expect_equivalent(rownames(x), paste(seq(100,1000,100)))))
  
  expect_equivalent(nrow(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = 100)$Si), 10)
  expect_equivalent(nrow(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = c(10, 100, 1000))$Si), 3)
  expect_equivalent(nrow(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = c(10, 100, 999))$Si), 4)
  
  expect_warning(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = 1), "Argument 'subsamp' should only contain values > 1.")
  expect_equivalent(nrow(suppressWarnings(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = c(10,1, 100, 900)))$Si), 4)
  expect_length(suppressWarnings(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, subsamp = 1))$Si, nparam)
})

test_that("bootstrapping and subsampling combined works as expected", {
  res_boot_sub <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper,
                       subsamp = 100, Nb = 500, full.output = T)
  
  expect_type(res_boot_sub, "list")
  expect_length(res_boot_sub, 21)
  invisible(lapply(res_boot_sub[c("ranking.boot", "Si.boot", "St.boot")], function(x) expect_true(all(dim(x) == c(10, 500, nparam)))))
  
  expect_type(res_boot_sub$B, "list")
  expect_length(res_boot_sub$B, 10)
  expect_true(all(dim(res_boot_sub$B$`500`) == c(500,500)))
  expect_true(all(dim(res_boot_sub$B$`700`) == c(700,500)))
})

test_that("argument 'verbose' works correctly", {
  expect_message(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, verbose = T),
                 "[ START SENSITIVITY ANALYSIS! ]")
  expect_silent(vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, verbose = F))
})
 
test_that("argument 'debug' produces files for debugging", {
  res_dbg <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, debug = T)
  expect_true(all(file.exists(paste(getwd(), c("vbsa_backup1.RData", "vbsa_backup2.RData", "vbsa_backup3.RData", "vbsa_backup4.RData"), sep="/"))))
  
  file.remove(paste(getwd(), c("vbsa_backup1.RData", "vbsa_backup2.RData", "vbsa_backup3.RData", "vbsa_backup4.RData"), sep="/"))
  res_dbg_off <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper, debug = F)
  expect_true(all(!file.exists(paste(getwd(), c("vbsa_backup1.RData", "vbsa_backup2.RData", "vbsa_backup3.RData", "vbsa_backup4.RData"), sep="/"))))
})

test_that("specified parameter names occur in the output", {
  names(pars_lower) <- c("p1", "p2", "p3")
  names(pars_upper) <- c("p1", "p2", "p3")
  res_pars <- vbsa(fn = "ishigami", lower = pars_lower, upper = pars_upper)
  invisible(lapply(res_pars[c("Si", "St", "ranking")], function(x) expect_named(x, c("p1", "p2", "p3"))))
})
