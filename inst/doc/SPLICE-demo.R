## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SPLICE)
set.seed(20201006)
ref_claim <- return_parameters()[1] # 200,000
time_unit <- return_parameters()[2] # 0.25

## -----------------------------------------------------------------------------
test_claims <- SynthETIC::test_claims_object

## -----------------------------------------------------------------------------
# major revisions
major <- claim_majRev_freq(test_claims)
major <- claim_majRev_time(test_claims, major)
major <- claim_majRev_size(major)

# minor revisions
minor <- claim_minRev_freq(test_claims)
minor <- claim_minRev_time(test_claims, minor)
minor <- claim_minRev_size(test_claims, major, minor)

# development of case estimates
test <- claim_history(test_claims, major, minor)
test_inflated <- claim_history(test_claims, major, minor,
                               base_inflation_vector = rep((1 + 0.02)^(1/4) - 1, times = 80))

# transactional data
test_incurred_dataset_noInf <- generate_incurred_dataset(test_claims, test)
test_incurred_dataset_inflated <- generate_incurred_dataset(test_claims, test_inflated)

# incurred cumulative triangles
incurred_inflated <- output_incurred(test_inflated, incremental = FALSE)

## -----------------------------------------------------------------------------
## paramfun input
# lambda as a function of claim size
no_majRev_param <- function(claim_size) {
  majRevNo_mean <- pmax(1, log(claim_size / 15000) - 2)
  c(lambda = majRevNo_mean)
}

## implementation and output
major_test <- claim_majRev_freq(
  test_claims, rfun = actuar::rztpois, paramfun = no_majRev_param)
# show the distribution of number of major revisions
table(unlist(major_test))

## -----------------------------------------------------------------------------
## paramfun input
# an extended parameter function
majRevNo_param <- function(claim_size, no_payment) {
  majRevNo_mean <- pmax(0, log(claim_size / 1500000)) + no_payment / 10
  c(lambda = majRevNo_mean)
}

## -----------------------------------------------------------------------------
## implementation and output
no_payments_vect <- unlist(test_claims$no_payments_list)
# sample the frequency of major revisions from zero-truncated Poisson
# with parameters above
major_test <- claim_majRev_freq(
  test_claims, rfun = actuar::rztpois, paramfun = majRevNo_param,
  no_payment = no_payments_vect)
# show the distribution of number of major revisions
table(unlist(major_test))

## -----------------------------------------------------------------------------
## input
# package default function for frequency of major revisions
dflt.majRev_freq_function <- function(
  n, claim_size, claim_size_benchmark = 0.075 * ref_claim) {
  
  # construct the range indicator
  test <- (claim_size > claim_size_benchmark)

  # if claim_size <= claim_size_benchmark
  # "small" claims assumed to have no major revisions except at notification
  no_majRev <- rep(1, n)
  # if claim_size is above the benchmark
  # probability of 2 major revisions, increases with claim size
  Pr2 <- 0.1 + 0.3 * 
    min(1, (claim_size[test] - 0.075 * ref_claim)/(0.925 * ref_claim))
  # probability of 3 major revisions, increases with claim size
  Pr3 <- 0.5 *
    min(1, max(0, claim_size[test] - 0.25 * ref_claim)/(0.75 * ref_claim))
  # probability of 1 major revision i.e. only one at claim notification
  Pr1 <- 1 - Pr2 - Pr3
  no_majRev[test] <- sample(
    c(1, 2, 3), size = sum(test), replace = T, prob = c(Pr1, Pr2, Pr3))
  
  no_majRev
}

## -----------------------------------------------------------------------------
## implementation and output
# simulate the number of major revisions
major <- claim_majRev_freq(
  claims = test_claims,
  rfun = dflt.majRev_freq_function
)

# show the distribution of number of major revisions
table(unlist(major))

# view the major revision history of the first claim in the 1st occurrence period
# note that the time and size of the major revisions are yet to be generated
major[[1]][[1]]

## -----------------------------------------------------------------------------
major_test <- claim_majRev_freq(
  claims = test_claims,
  claim_size_benchmark = 30000
)

## -----------------------------------------------------------------------------
## input
majRev_time_rfun <- function(n, min, max) {
  # n = number of major revisions of an individual claim
  majRev_time <- vector(length = n)
  majRev_time[1] <- 0 # first major revision at notification
  if (n > 1) {
    majRev_time[2:n] <- sort(stats::runif(n - 1, min, max))
  }
  
  return(majRev_time)
}
majRev_time_paramfun <- function(setldel, ...) {
  # setldel = settlement delay
  c(min = setldel/3, max = setldel)
}

## implementation and output
major_test <- claim_majRev_time(
  test_claims, major, rfun = majRev_time_rfun, paramfun = majRev_time_paramfun
)
major_test[[1]][[1]]

## -----------------------------------------------------------------------------
## package default function for time of major revisions
dflt.majRev_time_function <- function(
  # n = number of major revisions
  # setldel = settlement delay
  # penultimate_delay = time from claim notification to second last payment
  n, claim_size, setldel, penultimate_delay) {

  majRev_time <- rep(NA, times = n)
  
  # first revision at notification
  majRev_time[1] <- 0
  if (n > 1) {
    # if the claim has multiple major revisions
    # the probability of having the last revision exactly at the second last partial payment
    p <- 0.2 *
      min(1, max(0, (claim_size - ref_claim) / (14 * ref_claim)))
    at_second_last_pmt <- sample(c(0, 1), size = 1, replace = TRUE, prob = c(1-p, p))
    
    # does the last revision occur at the second last partial payment?
    if (at_second_last_pmt == 0) {
      # -> no revision at second last payment
      majRev_time[2:n] <- sort(rtri(n - 1, min = setldel/3, max = setldel, mode = setldel/3))
    } else {
      # -> yes, revision at second last payment
      majRev_time[n] <- penultimate_delay
      if (n > 2) {
        majRev_time[2:(n-1)] <- sort(
          rtri(n - 2, min = majRev_time[n]/3, max = majRev_time[n], mode = majRev_time[n]/3))
      }
    }
  }
  majRev_time
}

## -----------------------------------------------------------------------------
dflt.majRev_time_paramfun <- function(payment_delays, ...) {
  c(penultimate_delay = sum(payment_delays[1:length(payment_delays) - 1]),
    ...)
}

## -----------------------------------------------------------------------------
## implementation and output
major <- claim_majRev_time(
  claims = test_claims,
  majRev_list = major, # we will update the previous major list
  rfun = dflt.majRev_time_function,
  paramfun = dflt.majRev_time_paramfun
)

# view the major revision history of the first claim in the 1st occurrence period
# observe that we have now updated the time of major revisions
major[[1]][[1]]

## ---- eval=FALSE--------------------------------------------------------------
#  major <- claim_majRev_time(claims = test_claims, majRev_list = major)

## -----------------------------------------------------------------------------
## input
majRev_size_rfun <- function(n, shape, rate) {
  # n = number of major revisions of an individual claim
  majRev_size <- vector(length = n)
  majRev_size[1] <- 1 # first major revision at notification
  if (n > 1) {
    majRev_size[2:n] <- stats::rgamma(n - 1, shape, rate)
  }
  
  majRev_size
}

majRev_size_paramfun <- function(claim_size) {
  shape <- max(log(claim_size / 5000), 1)
  rate <- 10 / shape
  c(shape = shape, rate = rate)
}

## -----------------------------------------------------------------------------
## implementation and output
claim_size_vect <- unlist(test_claims$claim_size_list)
major_test <- claim_majRev_size(
  majRev_list = major,
  rfun = majRev_size_rfun,
  paramfun = majRev_size_paramfun,
  claim_size = claim_size_vect
)

# view the major revision history of the first claim in the 1st occurrence period
# observe that we have now updated the size of major revisions
major_test[[1]][[1]]

## -----------------------------------------------------------------------------
## input
# package default function for sizes of major revisions
dflt.majRev_size_function <- function(n) {
  majRev_factor <- rep(NA, times = n)
  # set revision size = 1 for first revision (i.e. the one at notification)
  majRev_factor[1] <- 1
  if (n > 1) {
    # if the claim has multiple major revisions
    majRev_factor[2] <- stats::rlnorm(n = 1, meanlog = 1.8, sdlog = 0.2)
    if (n > 2) {
      # the last revision factor depends on what happened at the second major revision
      mu <- 1 + 0.07 * (6 - majRev_factor[2])
      majRev_factor[3] <- stats::rlnorm(n = 1, meanlog = mu, sdlog = 0.1)
    }
  }

  majRev_factor
}

## implementation and output
major <- claim_majRev_size(
  majRev_list = major,
  rfun = dflt.majRev_size_function
)

# view the major revision history of the first claim in the 1st occurrence period
# observe that we have now updated the size of major revisions
major[[1]][[1]]

## -----------------------------------------------------------------------------
## input
# package default function for frequency of minor revisions at partial payments
dflt.minRev_freq_notatP_function <- function(n, setldel) {
  # setldel = settlement delay
  minRev_freq_notatP <- stats::rgeom(n, prob = 1 / (min(3, setldel/4) + 1))
  minRev_freq_notatP
}

## implementation and output
minor <- claim_minRev_freq(
  test_claims,
  prob_atP = 0.5,
  rfun_notatP = dflt.minRev_freq_notatP_function)

# view the minor revision history of the 10th claim in the 1st occurrence period
minor[[1]][[10]]

## ----eval=FALSE---------------------------------------------------------------
#  minRev_freq_notatP_paramfun <- function(setldel) {
#    c(prob = 1 / (min(3, setldel/4) + 1))
#  }
#  
#  minor <- claim_minRev_freq(
#    test_claims,
#    prob_atP = 0.5,
#    rfun_notatP = stats::rgeom,
#    paramfun_notatP = minRev_freq_notatP_paramfun)

## ----eval=FALSE---------------------------------------------------------------
#  minor <- claim_minRev_freq(claims = test_claims)

## -----------------------------------------------------------------------------
minRev_freq_notatP_paramfun <- function(setldel) {
  c(prob = 1 / (min(3, setldel/4) + 2))
}

minor_test <- claim_minRev_freq(
  test_claims,
  prob_atP = 0,
  rfun_notatP = stats::rgeom,
  paramfun_notatP = minRev_freq_notatP_paramfun)

minor_test[[1]][[10]]

## -----------------------------------------------------------------------------
## input
# package default function for time of minor revisions that do not coincide with a payment
dflt.minRev_time_notatP <- function(n, setldel) {
  sort(stats::runif(n, min = setldel/6, max = setldel))
}

## implementation and output
minor <- claim_minRev_time(
  claims = test_claims,
  minRev_list = minor, # we will update the previous minor list
  rfun_notatP = dflt.minRev_time_notatP
)

# view the minor revision history of the 10th claim in the 1st occurrence period
# observe that we have now updated the time of minor revisions
minor[[1]][[10]]

## -----------------------------------------------------------------------------
## input
minRev_time_notatP_rfun <- function(n, setldel) {
  # n = number of minor revisions
  # setldel = settlement delay
  sort(rtri(n, min = setldel/6, max = setldel, mode = setldel/6))
}

## implementation and output
minor_test <- claim_minRev_time(
  claims = test_claims,
  minRev_list = minor, # we will update the previous minor list
  rfun_notatP = minRev_time_notatP_rfun
)

# view the minor revision history of the 10th claim in the 1st occurrence period
# observe that we have now updated the time of minor revisions
minor_test[[1]][[10]]

## -----------------------------------------------------------------------------
## input
# package default function for the size of minor revisions
dflt.minRev_size <- function(
  # n = number of minor revisions
  # minRev_time = epochs of the minor revisions (from claim notification)
  # majRev_time_2nd = epoch of 2nd major revision (from claim notification)
  # setldel = settlement delay
  n, minRev_time, majRev_time_2nd, setldel) {

  k <- length(minRev_time)
  minRev_factor <- vector(length = k)

  if (k >= 1) {
    for (i in 1:k) {
      curr <- minRev_time[i]
      if (curr <= setldel/3) {
        meanlog <- 0.15
      } else if (curr <= (2/3) * setldel) {
        meanlog <- 0
      } else {
        meanlog <- -0.1
      }
      sdlog <- ifelse(curr > majRev_time_2nd, 0.05, 0.1)
      minRev_factor[i] <- stats::rlnorm(n = 1, meanlog, sdlog)
    }
  }

  minRev_factor
}

## -----------------------------------------------------------------------------
# parameter function for minor revision at payments
minRev_size_param_atP <- function(major, minor, setldel) {
  list(minRev_time = minor$minRev_time_atP,
       majRev_time_2nd = ifelse(
         # so it always holds minRev_time < majRev_time_2nd
         is.na(major$majRev_time[2]), setldel + 1, major$majRev_time[2]),
       setldel = setldel)
}

# parameter function for minor revisions NOT at payments
minRev_size_param_notatP <- function(major, minor, setldel) {
  list(minRev_time = minor$minRev_time_notatP,
       majRev_time_2nd = ifelse(
         # so it always holds minRev_time < majRev_time_2nd
         is.na(major$majRev_time[2]), setldel + 1, major$majRev_time[2]),
       setldel = setldel)
}

## -----------------------------------------------------------------------------
## implementation and output
minor <- claim_minRev_size(
  claims = test_claims,
  majRev_list = major,
  minRev_list = minor,
  rfun = dflt.minRev_size,
  paramfun_atP = minRev_size_param_atP,
  paramfun_notatP = minRev_size_param_notatP
)

# view the minor revision history of the 10th claim in the 1st occurrence period
# observe that we have now updated the size of minor revisions
minor[[1]][[10]]

## -----------------------------------------------------------------------------
## input
paramfun_atP <- function(claim_size, ...) {
  c(min = pmin(1, pmax(log(claim_size / 15000), 0.5)),
    max = pmin(1, pmax(log(claim_size / 15000), 0.5)) + 1)
}
paramfun_notatP <- paramfun_atP

## implementation and output
claim_size_vect <- unlist(test_claims$claim_size_list)
minor_test <- claim_minRev_size(
  test_claims, major, minor,
  rfun = stats::runif, paramfun_atP, paramfun_notatP,
  claim_size = claim_size_vect)
minor_test[[1]][[10]]

## -----------------------------------------------------------------------------
# exclude inflation (by default)
result <- claim_history(test_claims, major, minor)
# include inflation
result_inflated <- claim_history(
  test_claims, major, minor, 
  base_inflation_vector = rep((1 + 0.02)^(1/4) - 1, times = 80))

## -----------------------------------------------------------------------------
data <- generate_incurred_dataset(test_claims_object, result)
str(data)
head(data, n = 9)

data_inflated <- generate_incurred_dataset(test_claims_object, result_inflated)
str(data_inflated)
head(data_inflated, n = 9)

## -----------------------------------------------------------------------------
str(test_incurred_dataset_noInf)
str(test_incurred_dataset_inflated)

## -----------------------------------------------------------------------------
square_inc <- output_incurred(result)
square_cum <- output_incurred(result, incremental = F)
square_inflated_inc <- output_incurred(result_inflated)
square_inflated_cum <- output_incurred(result_inflated, incremental = F)

yearly_inc <- output_incurred(result, aggregate_level = 4)
yearly_cum <- output_incurred(result, aggregate_level = 4, incremental = F)
yearly_cum

# apply standard actuarial reserving techniques using the `ChainLadder` package
# selected <- attr(ChainLadder::ata(yearly_cum), "vwtd")

## -----------------------------------------------------------------------------
# output the past cumulative triangle
cumtri <- output_incurred(result, aggregate_level = 4, 
                          incremental = FALSE, future = FALSE)
# calculate the age to age factors
selected <- attr(ChainLadder::ata(cumtri), "vwtd")
# complete the triangle
CL_prediction <- cumtri
J <- nrow(cumtri)
for (i in 2:J) {
  for (j in (J - i + 2):J) {
    CL_prediction[i, j] <- CL_prediction[i, j - 1] * selected[j - 1]
  }
}

CL_prediction

