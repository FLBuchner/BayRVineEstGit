library(cmdstanr)

###############################################################################
# functions needed for sample_tau function
###############################################################################
normalizeRVineMatrix <- function(RVM) {

  oldOrder <- diag(RVM$Matrix)
  Matrix <- reorderRVineMatrix(RVM$Matrix)

  return(RVineMatrix(Matrix,
                     RVM$family,
                     RVM$par,
                     RVM$par2,
                     names = rev(RVM$names[oldOrder]),
                     check.pars = FALSE))
}

reorderRVineMatrix <- function(Matrix, oldOrder = NULL) {

  if (length(oldOrder) == 0) {
    oldOrder <- diag(Matrix)
  }
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)

  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
  }
  return(Matrix)
}

##################################################################################################
sample_tau = function(RVM, udata, iter_w = 150, iter_s = 2000, init = 0.8, s_bound = 0.95, refresh_s = 100, adapt_delta_s = 0.8){

  STAN <- cmdstan_model("STAN.stan")
  n <- nrow(udata)
  d <- dim(RVM$Matrix)[1]

  o <- diag(RVM$Matrix)
  if (any(o != length(o):1)) {
    oldRVM   <- RVM
    RVM      <- normalizeRVineMatrix(RVM)
    udata <- udata[, o[length(o):1]]
  }

  #############
  ## From the package: the different parameters/inputs we need)
  #############

  family <- as.vector(RVM$family)
  family[is.na(family)] <- 0
  par2 <- as.vector(RVM$par2)
  par2[is.na(par2)] <- 0
  condirect   <- as.vector(as.numeric(RVM$CondDistr$direct))
  conindirect <- as.vector(as.numeric(RVM$CondDistr$indirect))
  maxmat      <- as.vector(RVM$MaxMat)
  matri       <- as.vector(RVM$Matrix)
  matri[is.na(matri)]   <- 0
  maxmat[is.na(maxmat)] <- 0
  condirect[is.na(condirect)]     <- 0
  conindirect[is.na(conindirect)] <- 0

  famvec <- famVector(RVM)
  L_bound   <- c()
  U_bound   <- c()

  for (i in 1:(d*(d-1)/2)) {
    if (famvec[i] %in% c(3, 13, 4, 14, 6, 16)) {
      L_bound[i]   <- 0
      U_bound[i]   <- s_bound
    } else if (famvec[i] %in% c(23, 33, 24, 34, 26, 36)) {
      L_bound[i]   <- -s_bound
      U_bound[i]   <- 0
    } else {
      L_bound[i]   <- -s_bound
      U_bound[i]   <- s_bound
    }
  }

  #############
  ## Sample from STAN
  #############

  x_grid <- c(-10^(5:2), seq(-36, 36, l = 100), 10^(2:5))
  # frankTauVals <- 1 - 4/frankParGrid + 4/frankParGrid * copula::debye1(frankParGrid)
  y_grid <- c(-0.99996000, -0.99960007, -0.99600658, -0.96065797, -0.89396585, -0.89188641, -0.88972402, -0.88747361,
              -0.88512973, -0.88268644, -0.88013730, -0.87747533, -0.87469288, -0.87178163, -0.86873246, -0.86553538,
              -0.86217942, -0.85865252, -0.85494133, -0.85103114, -0.84690560, -0.84254657, -0.83793380, -0.83304468,
              -0.82785385, -0.82233280, -0.81644934, -0.81016705, -0.80344454, -0.79623459, -0.78848313, -0.78012801,
              -0.77109744, -0.76130821, -0.75066334, -0.73904940, -0.72633308, -0.71235706, -0.69693500, -0.67984550,
              -0.66082501, -0.63955977, -0.61567712, -0.58873731, -0.55822774, -0.52356346, -0.48410024, -0.43916961,
              -0.38814802, -0.33057147, -0.26629772, -0.19569472, -0.11979813, -0.04035073,  0.04035073,  0.11979813,
               0.19569472,  0.26629772,  0.33057147,  0.38814802,  0.43916961,  0.48410024,  0.52356346,  0.55822774,
               0.58873731,  0.61567712,  0.63955977,  0.66082501,  0.67984550,  0.69693500,  0.71235706,  0.72633308,
               0.73904940,  0.75066334,  0.76130821,  0.77109744,  0.78012801,  0.78848313,  0.79623459,  0.80344454,
               0.81016705,  0.81644934,  0.82233280,  0.82785385,  0.83304468,  0.83793380,  0.84254657,  0.84690560,
               0.85103114,  0.85494133,  0.85865252,  0.86217942,  0.86553538,  0.86873246,  0.87178163,  0.87469288,
               0.87747533,  0.88013730,  0.88268644,  0.88512973,  0.88747361,  0.88972402,  0.89188641,  0.89396585,
               0.96065797,  0.99600658,  0.99960007,  0.99996000)

  data_stan <- list(n=n, d=d, L_bound=L_bound, U_bound=U_bound, para2=par2, udata=udata,
                 family=as.integer(family), maxmat=as.integer(maxmat), matri=as.integer(matri),
                 condirect=as.integer(condirect), conindirect=as.integer(conindirect),
                 x_grid=x_grid, y_grid=y_grid)

  fit <- STAN$sample(iter_warmup = iter_w, iter_sampling = iter_s, refresh = refresh_s,
                     data=data_stan, seed = 317, chains = 4, parallel_chains = 4,
                     save_warmup = TRUE, adapt_delta = adapt_delta_s, init = init)

  return(fit)
}
