library(kde1d)
library(bayesplot)
library(ggpubr)
library(VineCopula)
library(latex2exp)

estimate_mode=function(x) {
  tryCatch( { d=density(x)
  return(d$x[which.max(d$y)])
  }, error = function(e) {print(e)
    return(mean(x))})
}


tauCIMLE <- function(data, fit, quantile) {
  H <- 0
  n <- nrow(data)

  for (i in 1:n) {
    H <- H + RVineHessian(data[i,], fit$RVM)$hessian
  }

  H <- H/n

  g_val1  <- par2tauderC(fit$RVM$par[5,4])
  g_val2  <- par2tauderG(fit$RVM$par[5,3])
  g_val3  <- par2tauderJ(fit$RVM$par[5,2])
  g_val4  <- par2tauderC(fit$RVM$par[5,1])
  g_val5  <- par2tauderG(fit$RVM$par[4,3])
  g_val6  <- par2tauderJ(fit$RVM$par[4,2])
  g_val7  <- par2tauderC(fit$RVM$par[4,1])
  g_val8  <- par2tauderC(fit$RVM$par[3,2])
  g_val9  <- par2tauderG(fit$RVM$par[3,1])
  g_val10 <- par2tauderG(fit$RVM$par[2,1])

  S21  <- g_val1^2 * abs(H[1,1])
  S22  <- g_val2^2 * abs(H[2,2])
  S23  <- g_val3^2 * abs(H[3,3])
  S24  <- g_val4^2 * abs(H[4,4])
  S25  <- g_val5^2 * abs(H[5,5])
  S26  <- g_val6^2 * abs(H[6,6])
  S27  <- g_val7^2 * abs(H[7,7])
  S28  <- g_val8^2 * abs(H[8,8])
  S29  <- g_val9^2 * abs(H[9,9])
  S210 <- g_val10^2 * abs(H[10,10])

  q <- qnorm(quantile)

  tau_mle <- c(fit$RVM$tau[5,4], fit$RVM$tau[5,3], fit$RVM$tau[5,2],
               fit$RVM$tau[5,1], fit$RVM$tau[4,3], fit$RVM$tau[4,2],
               fit$RVM$tau[4,1], fit$RVM$tau[3,2], fit$RVM$tau[3,1],
               fit$RVM$tau[2,1])
  L <- c(fit$RVM$tau[5,4] - q*sqrt(S21)/sqrt(n), fit$RVM$tau[5,3] - q*sqrt(S22)/sqrt(n),
         fit$RVM$tau[5,2] - q*sqrt(S23)/sqrt(n), fit$RVM$tau[5,1] - q*sqrt(S24)/sqrt(n),
         fit$RVM$tau[4,3] - q*sqrt(S25)/sqrt(n), fit$RVM$tau[4,2] - q*sqrt(S26)/sqrt(n),
         fit$RVM$tau[4,1] - q*sqrt(S27)/sqrt(n), fit$RVM$tau[3,2] - q*sqrt(S28)/sqrt(n),
         fit$RVM$tau[3,1] - q*sqrt(S29)/sqrt(n), fit$RVM$tau[2,1] - q*sqrt(S210)/sqrt(n))
  U <- c(fit$RVM$tau[5,4] + q*sqrt(S21)/sqrt(n), fit$RVM$tau[5,3] + q*sqrt(S22)/sqrt(n),
         fit$RVM$tau[5,2] + q*sqrt(S23)/sqrt(n), fit$RVM$tau[5,1] + q*sqrt(S24)/sqrt(n),
         fit$RVM$tau[4,3] + q*sqrt(S25)/sqrt(n), fit$RVM$tau[4,2] + q*sqrt(S26)/sqrt(n),
         fit$RVM$tau[4,1] + q*sqrt(S27)/sqrt(n), fit$RVM$tau[3,2] + q*sqrt(S28)/sqrt(n),
         fit$RVM$tau[3,1] + q*sqrt(S29)/sqrt(n), fit$RVM$tau[2,1] + q*sqrt(S210)/sqrt(n))

  return(c(L, U))
}

par2tauderG <- function(par){
  return(1/(par^2))
}

par2tauderC <- function(par){
  return(-1/(2+par)^2)
}

par2tauderJ <- function(par){
  return(((par-2)*(1 - trigamma(1/par)/par^2 - trigamma((2+par)/(2*par)) + par)
          - (-2 - 2 * digamma(1) + 2 * log(2) + digamma(1/par) + digamma((2+par)/(2*par)) + par)) / (-2+par)^2)
}

name_tau <- function(RVM, edge = FALSE) {
  M <- RVM$Matrix
  d <- dim(RVM$Matrix)[1]
  names <- c()
  tau <- ""
  if (edge == FALSE) {
    tau <- "tau_"
  }

  if (d >= 10) {
    for (i in d:2) {
      for (j in (i-1):1) {
        name <- paste0(tau, M[i,j], ",", M[j,j])
        if (i < d) {
          k <- d - i
          name <- paste0(name, ";")
          for (l in 1:k) {
            if (l == 1) {
              name <- paste0(name, M[i+l,j])
            } else {
              name <- paste0(name, ",", M[i+l,j])
            }
          }
        }
        names <- c(names, name)
      }
    }
  } else {
    for (i in d:2) {
      for (j in (i-1):1) {
        name <- paste0(tau, M[i,j], M[j,j])
        if (i < d) {
          k <- d - i
          name <- paste0(name, ";")
          for (l in 1:k) {
            name <- paste0(name, M[i+l,j])
          }
        }
        names <- c(names, name)
      }
    }
  }
  return(names)
}

name_tau_split <- function(RVM) {
  M <- RVM$Matrix
  d <- dim(RVM$Matrix)[1]
  conditioned  <- c()
  conditioning <- c()

  for (i in d:2) {
    for (j in (i-1):1) {
      conditioned <- c(conditioned, paste0(M[i,j], ",", M[j,j]))
      if (i < d) {
        aux <- c()
        k <- d - i
        for (l in 1:k) {
          if (l == 1) {
            aux <- paste0(aux, M[i+l,j])
          } else {
            aux <- paste0(aux, ",", M[i+l,j])
          }
        }
        conditioning <- c(conditioning, aux)
      } else {
        conditioning <- c(conditioning, "")
      }
    }
  }
  df <- cbind(conditioned, conditioning)
  return(df)
}

tauVector <- function(RVM) {
  d <- dim(RVM$Matrix)[1]
  M <- RVM$tau
  tauV <- c()

  for (i in d:2) {
    for (j in (i-1):1) {
      tauV <- c(tauV, M[i,j])
    }
  }
  return(tauV)
}

famVector <- function(RVM) {
  d <- dim(RVM$Matrix)[1]
  M <- RVM$family
  V <- c()

  for (i in d:2) {
    for (j in (i-1):1) {
      V <- c(V, M[i,j])
    }
  }
  return(V)
}

famVectorName <- function(RVM) {
  d <- dim(RVM$Matrix)[1]
  M <- RVM$family
  famV <- c()

  for (i in d:2) {
    for (j in (i-1):1) {
      famV <- c(famV, M[i,j])
    }
  }

  for (i in 1:length(famV)) {
    if (famV[i] == 0) {
      famV[i] <- "I"
    } else if (famV[i] == 1) {
      famV[i] <- "N"
    } else if (famV[i] == 2) {
      famV[i] <- "t"
    } else if (famV[i] == 3) {
      famV[i] <- "C"
    } else if (famV[i] == 4) {
      famV[i] <- "G"
    } else if (famV[i] == 5) {
      famV[i] <- "F"
    } else if (famV[i] == 6) {
      famV[i] <- "J"
    } else if (famV[i] == 13) {
      famV[i] <- "SC"
    } else if (famV[i] == 14) {
      famV[i] <- "SG"
    } else if (famV[i] == 16) {
      famV[i] <- "SJ"
    } else if (famV[i] == 23) {
      famV[i] <- "C90"
    } else if (famV[i] == 24) {
      famV[i] <- "G90"
    } else if (famV[i] == 26) {
      famV[i] <- "J90"
    } else if (famV[i] == 33) {
      famV[i] <- "C270"
    } else if (famV[i] == 34) {
      famV[i] <- "G270"
    } else if (famV[i] == 36) {
      famV[i] <- "J270"
    }
  }
  return(famV)
}

dtau <- function(d) {
  tau <- c()
  for(i in 1:(d*(d-1)/2)) {
    tau <- c(tau, paste0("tau[", i, "]"))
  }
  return(tau)
}

PIT_kde <- function(data) {
  udata <- data

  for (i in colnames(data)) {
    fit <- kde1d(data[,i])
    udata[,i] <- pkde1d(data[,i], fit)
  }

  return(udata)
}

parMatrix <- function(par, d) {
  mat <- matrix(0, d, d)
  k <- 0

  for (i in d:2) {
    mat[i, 1:(i-1)] <- par[(k+i-1):(k+1)]
    k <- k + (i-1)
  }

  return(mat)
}

medianTauMatrix <- function(fit, d, tau) {
  mat <- matrix(0, d, d)
  k <- 0

  for (i in d:2) {
    mat[i, 1:(i-1)] <- fit$summary(tau)[(k+i-1):(k+1),3][[1]]
    k <- k + (i-1)
  }

  return(mat)
}

tauIntervalsPlot <- function(tau, tau_title, true_tau, fit_list, S=20) {
  #color_scheme_set(c("#e37222","#a2ad00","#98c6ea","#0065bd","#C680BB","#EF9067"))

  tau_array <- array(dim = c(500, 4, S), dimnames = list("iteration" = 1:500, "chain" = 1:4, "variable" = 1:S))

  for (i in 1:S) {
    tau_array[,,i] <- fit_list[[i]]$draws(tau)
  }

  p <- mcmc_intervals(tau_array) +
    vline_at(true_tau, color = "red") +
    ggtitle(paste("credible intervals for ", tau_title)) +
    ylab("Number of simulated data set")

  return(p)
}

tauDensityPlot <- function(fit_list, true_tau, sim_tau, tau, N=3, nrow = 4){
  #color_scheme_set(c("#e37222","#C680BB","#98c6ea","#0065bd","#a2ad00","#EF9067"))
  plot_list <- list()

  for (i in 1:N) {
    plot_list[[i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_dens(fit_list[[i]]$draws("tau[1]")) +
          geom_vline(aes(xintercept = true_tau[1], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[1], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[1]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                              values = c("original data" = 1,
                                         "median" = 1,
                                         "vinereg" = 1),
                              guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          ylab("tau_12") + xlab("")
      } else {
        mcmc_dens(fit_list[[i]]$draws("tau[1]")) +
          geom_vline(aes(xintercept = true_tau[1], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[1], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[1]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                                values = c("original data" = 1,
                                           "median" = 1,
                                           "vinereg" = 1),
                                guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          xlab("")
      }
    })
  }
  for (i in 1:N) {
    plot_list[[N+i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_dens(fit_list[[i]]$draws("tau[2]")) +
          geom_vline(aes(xintercept = true_tau[2], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[2], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[2]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                                values = c("original data" = 1,
                                           "median" = 1,
                                           "vinereg" = 1),
                                guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          ylab("tau_23") + xlab("")
      } else {
        mcmc_dens(fit_list[[i]]$draws("tau[2]")) +
          geom_vline(aes(xintercept = true_tau[2], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[2], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[2]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                                values = c("original data" = 1,
                                           "median" = 1,
                                           "vinereg" = 1),
                                guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          xlab("")
      }
    })
  }
  for (i in 1:N) {
    plot_list[[2*N+i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_dens(fit_list[[i]]$draws("tau[3]")) +
          geom_vline(aes(xintercept = true_tau[3], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[3], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[3]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                                values = c("original data" = 1,
                                           "median" = 1,
                                           "vinereg" = 1),
                                guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          ylab("tau_13_2") + xlab(paste("Replication ", i, sep = ""))
      } else {
        mcmc_dens(fit_list[[i]]$draws("tau[3]")) +
          geom_vline(aes(xintercept = true_tau[3], linetype = "original data"), color = "red") +
          geom_vline(aes(xintercept = sim_tau[3], linetype = "vinereg"), color = "green") +
          geom_vline(aes(xintercept = median(fit_list[[i]]$draws("tau[3]")), linetype = "median"), color = "#C680BB") +
          scale_linetype_manual(name = "tau",
                                values = c("original data" = 1,
                                           "median" = 1,
                                           "vinereg" = 1),
                                guide = guide_legend(override.aes = list(colour = c("#C680BB", "red", "green")))) +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          xlab(paste("Replication ", i, sep = ""))
      }
    })
  }

  return(ggarrange(plotlist = plot_list, ncol = N, nrow = nrow, common.legend = TRUE, legend = "right"))
}

tauTracePlot <- function(fit_list, true_tau, sim_tau, N=4, nrow = 3){
  #color_scheme_set(c("#98c6ea","#C680BB","#e37222","#a2ad00","#0065bd","#EF9067"))
  plot_list <- list()

  for (i in 1:N) {
    plot_list[[i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_trace(fit_list[[i]]$draws("tau[1]")) +
          hline_at(true_tau[1], color = "red") +
          hline_at(sim_tau[1], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          theme(axis.title.x=element_blank()) +
          ylab("tau_12")
      } else {
        mcmc_trace(fit_list[[i]]$draws("tau[1]")) +
          hline_at(true_tau[1], color = "red") +
          hline_at(sim_tau[1], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank())
      }
    })
  }
  for (i in 1:N) {
    plot_list[[N+i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_trace(fit_list[[i]]$draws("tau[2]")) +
          hline_at(true_tau[2], color = "red") +
          hline_at(sim_tau[2], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          theme(axis.title.x=element_blank()) +
          ylab("tau_23")
      } else {
        mcmc_trace(fit_list[[i]]$draws("tau[2]")) +
          hline_at(true_tau[2], color = "red") +
          hline_at(sim_tau[2], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank())
      }
    })
  }
  for (i in 1:N) {
    plot_list[[2*N+i]] <- local({
      i <- i
      if (i == 1) {
        mcmc_trace(fit_list[[i]]$draws("tau[3]")) +
          hline_at(true_tau[3], color = "red") +
          hline_at(sim_tau[3], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          ylab("tau_13_2") +
          xlab(paste("Replication ", i, sep = ""))
      } else {
        mcmc_trace(fit_list[[i]]$draws("tau[3]")) +
          hline_at(true_tau[3], color = "red") +
          hline_at(sim_tau[3], color = "green") +
          scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
          theme(axis.title.y=element_blank()) +
          xlab(paste("Replication ", i, sep = ""))
      }
    })
  }

  return(ggarrange(plotlist = plot_list, ncol = N, nrow = nrow, common.legend = TRUE, legend = "right"))
}
