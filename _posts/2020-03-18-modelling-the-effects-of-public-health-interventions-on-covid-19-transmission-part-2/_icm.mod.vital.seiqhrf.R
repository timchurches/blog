departures.seiqhrf.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  groups <- dat$param$groups
  group <- dat$attr$group


  # Susceptible departures ------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$ds.rate, dat$param$ds.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$ds.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$ds.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$ds.flow.g2[at] <- nDeparturesG2
    }
  }

# Exposed Departures ---------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "e")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$de.rate, dat$param$de.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$de.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$de.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$de.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$de.flow.g2[at] <- nDeparturesG2
    }
  }


  # Infected Departures ---------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$di.rate, dat$param$di.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$di.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$di.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$di.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$di.flow.g2[at] <- nDeparturesG2
    }
  }

  # Quarantined Departures ---------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "q")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$dq.rate, dat$param$dq.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$dq.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$dq.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$dq.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$dq.flow.g2[at] <- nDeparturesG2
    }
  }
  
  # Hospitalised Departures ---------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "h")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$dh.rate, dat$param$dh.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$dh.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$dh.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$dh.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$dh.flow.g2[at] <- nDeparturesG2
    }
  }
  

  # Recovered Departures --------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "r")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$dr.rate, dat$param$dr.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        nDepartures <- sum(group[idsDpt] == 1)
        nDeparturesG2 <- sum(group[idsDpt] == 2)
        dat$attr$active[idsDpt] <- 0
      }
    } else {
      nDepartures <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDepartures)] <- 0
      if (groups == 2) {
        nDeparturesG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeparturesG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$dr.flow <- c(0, nDepartures)
    if (groups == 2) {
      dat$epi$dr.flow.g2 <- c(0, nDeparturesG2)
    }
  } else {
    dat$epi$dr.flow[at] <- nDepartures
    if (groups == 2) {
      dat$epi$dr.flow.g2[at] <- nDeparturesG2
    }
  }

  return(dat)
}

arrivals.seiqhrf.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  a.rate.g2 <- dat$param$a.rate.g2
  a.rand <- dat$control$a.rand
  groups <- dat$param$groups
  nOld <- dat$epi$num[at - 1]

  # checking params, should be in control.icm or params.icm eventually
  type <- dat$control$type
  nsteps <- dat$control$nsteps
  
  if (!(length(a.rate) == 1 || length(a.rate == nsteps))) {
    stop("Length of a.rate must be 1 or the value of nsteps")
  }
  if (!is.null(a.rate.g2) && 
      !(length(a.rate.g2) == 1 || length(a.rate.g2 == nsteps))) {
    stop("Length of a.rate.g2 must be 1 or the value of nsteps")
  }
  
  a.prop.e <- dat$param$a.prop.e
  if (!(length(a.prop.e) == 1 || length(a.prop.e == nsteps))) {
    stop("Length of a.prop.e must be 1 or the value of nsteps")
  }
  a.prop.i <- dat$param$a.prop.i
  if (!(length(a.prop.i) == 1 || length(a.prop.i == nsteps))) {
    stop("Length of a.prop.i must be 1 or the value of nsteps")
  }
  a.prop.q <- dat$param$a.prop.q
  if (!(length(a.prop.q) == 1 || length(a.prop.q == nsteps))) {
    stop("Length of a.prop.q must be 1 or the value of nsteps")
  }

  a.prop.e.g2 <- dat$param$a.prop.e.g2
  if (!is.null(a.prop.e.g2) &&
      !(length(a.prop.e.g2) == 1 || length(a.prop.e.g2 == nsteps))) {
    stop("Length of a.prop.e.g2 must be 1 or the value of nsteps")
  }
  a.prop.i.g2 <- dat$param$a.prop.i.g2
  if (!is.null(a.prop.i.g2) &&
      !(length(a.prop.i.g2) == 1 || length(a.prop.i.g2 == nsteps))) {
    stop("Length of a.prop.i.g2 must be 1 or the value of nsteps")
  }
  a.prop.q.g2 <- dat$param$a.prop.q.g2
  if (!is.null(a.prop.q.g2) &&
      !(length(a.prop.q.g2) == 1 || length(a.prop.q.g2 == nsteps))) {
    stop("Length of a.prop.q.g2 must be 1 or the value of nsteps")
  }
  
  # Process -----------------------------------------------------------------
  nArrivals <- nArrivals.g2 <- 0

  if (groups == 1) {
    if (a.rand == TRUE) {
      nArrivals <- sum(rbinom(nOld, 1, a.rate))
    }
    if (a.rand == FALSE) {
      nArrivals <- round(nOld * a.rate)
    }
  }
  if (groups == 2) {
    nOldG2 <- dat$epi$num.g2[at - 1]
    if (a.rand == TRUE) {
      if (is.na(a.rate.g2)) {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivals.g2 <- sum(rbinom(nOld, 1, a.rate))
      } else {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivals.g2 <- sum(rbinom(nOldG2, 1, a.rate.g2))
      }
    }
    if (a.rand == FALSE) {
      if (is.na(a.rate.g2)) {
        nArrivals <- round(nOld * a.rate)
        nArrivals.g2 <- round(nOld * a.rate)
      } else {
        nArrivals <- round(nOld * a.rate)
        nArrivals.g2 <- round(nOldG2 * a.rate.g2)
      }
    }
  }

  
  ## Set attributes
  totArrivals <- 0
  totArrivals.g2 <- 0
  
  # partition arrivals into compartments
  if (length(a.prop.e) > 1) {
    nArrivals.e <- round(nArrivals*(a.prop.e[at]))
    totArrivals <- totArrivals + nArrivals.e
    if (!is.null(a.prop.e.g2)) {
      nArrivals.e.g2 <- round(nArrivals.g2*(a.prop.e.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.e.g2
    } else {
      nArrivals.e.g2 <- round(nArrivals.g2*(a.prop.e.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.e.g2
    }
  } else {
    nArrivals.e <- round(nArrivals*a.prop.e)
    totArrivals <- totArrivals + nArrivals.e
    if (!is.null(a.prop.e.g2)) {
      nArrivals.e.g2 <- round(nArrivals.g2*(a.prop.e.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.e.g2
    } else {
      nArrivals.e.g2 <- round(nArrivals.g2*(a.prop.e.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.e.g2
    }
  }

  if (length(a.prop.i) > 1) {
    nArrivals.i <- round(nArrivals*(a.prop.i[at]))
    totArrivals <- totArrivals + nArrivals.i
    if (!is.null(a.prop.i.g2)) {
      nArrivals.i.g2 <- round(nArrivals.g2*(a.prop.i.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.i.g2
    } else {
      nArrivals.i.g2 <- round(nArrivals.g2*(a.prop.i.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.i.g2
    }
  } else {
    nArrivals.i <- round(nArrivals*a.prop.i)
    totArrivals <- totArrivals + nArrivals.i
    if (!is.null(a.prop.i.g2)) {
      nArrivals.i.g2 <- round(nArrivals.g2*(a.prop.i.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.i.g2
    } else {
      nArrivals.i.g2 <- round(nArrivals.g2*(a.prop.i.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.i.g2
    }
  }

  if (length(a.prop.q) > 1) {
    nArrivals.q <- round(nArrivals*(a.prop.q[at]))
    totArrivals <- totArrivals + nArrivals.q
    if (!is.null(a.prop.q.g2)) {
      nArrivals.q.g2 <- round(nArrivals.g2*(a.prop.q.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.q.g2
    } else {
      nArrivals.q.g2 <- round(nArrivals.g2*(a.prop.q.g2[at]))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.q.g2
    }
  } else {
    nArrivals.q <- round(nArrivals*a.prop.q)
    totArrivals <- totArrivals + nArrivals.q
    if (!is.null(a.prop.q.g2)) {
      nArrivals.q.g2 <- round(nArrivals.g2*(a.prop.q.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.q.g2
    } else {
      nArrivals.q.g2 <- round(nArrivals.g2*(a.prop.q.g2))
      totArrivals.g2 <- totArrivals.g2 + nArrivals.q.g2
    }
  }

  # debug
  print("totArrivals:")
  print(totArrivals)
  print("totArrivals.g2:")
  print(totArrivals.g2)
  print("----")
  
  # group 1
  dat$attr$active <- c(dat$attr$active, rep(1, totArrivals))
  dat$attr$group <- c(dat$attr$group, rep(1, totArrivals))
  dat$attr$status <- c(dat$attr$status,
                       rep("e", nArrivals.e),
                       rep("i", nArrivals.i),
                       rep("q", nArrivals.q),
                       rep("s", totArrivals - nArrivals.e - nArrivals.i - nArrivals.q))
  dat$attr$expTime <- c(dat$attr$expTime, rep(NA, totArrivals))
  dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArrivals))
  dat$attr$quarTime <- c(dat$attr$quarTime, rep(NA, totArrivals))
  dat$attr$hospTime <- c(dat$attr$ihospTime, rep(NA, totArrivals))
  dat$attr$recovTime <- c(dat$attr$recovTime, rep(NA, totArrivals))
  dat$attr$fatTime <- c(dat$attr$fatTime, rep(NA, totArrivals))

  # group 2
  if (length(totArrivals.g2) > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, totArrivals.g2))
    dat$attr$group <- c(dat$attr$group, rep(2, totArrivals.g2))
    dat$attr$status <- c(dat$attr$status,
                         rep("e", nArrivals.e.g2),
                         rep("i", nArrivals.i.g2),
                         rep("q", nArrivals.q.g2),
                         rep("s", totArrivals.g2 - nArrivals.e.g2 - 
                                   nArrivals.i.g2 - nArrivals.q.g2))
    dat$attr$expTime <- c(dat$attr$expTime, rep(NA, totArrivals.g2))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArrivals.g2))
    dat$attr$quarTime <- c(dat$attr$quarTime, rep(NA, totArrivals.g2))
    dat$attr$hospTime <- c(dat$attr$ihospTime, rep(NA, totArrivals.g2))
    dat$attr$recovTime <- c(dat$attr$recovTime, rep(NA, totArrivals.g2))
    dat$attr$fatTime <- c(dat$attr$fatTime, rep(NA, totArrivals.g2))
  }
  
  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$a.flow <- c(0, totArrivals)
    dat$epi$a.e.flow <- c(0, nArrivals.e)
    dat$epi$a.i.flow <- c(0, nArrivals.i)
    dat$epi$a.q.flow <- c(0, nArrivals.q)
  } else {
    dat$epi$a.flow[at] <- totArrivals
    dat$epi$a.e.flow[at] <- c(0, nArrivals.e)
    dat$epi$a.i.flow[at] <- c(0, nArrivals.i)
    dat$epi$a.q.flow[at] <- c(0, nArrivals.q)
  }
  if (length(totArrivals.g2) > 0 && dat$param$groups == 2) {
    if (at == 2) {
      dat$epi$a.flow.g2 <- c(0, totArrivals.g2)
      dat$epi$a.e.flow.g2 <- c(0, nArrivals.e.g2)
      dat$epi$a.i.flow.g2 <- c(0, nArrivals.i.g2)
      dat$epi$a.q.flow.g2 <- c(0, nArrivals.q.g2)
    } else {
      dat$epi$a.flow.g2[at] <- totArrivals.g2
      dat$epi$a.e.flow.g2[at] <- c(0, nArrivals.e.g2)
      dat$epi$a.i.flow.g2[at] <- c(0, nArrivals.i.g2)
      dat$epi$a.q.flow.g2[at] <- c(0, nArrivals.q.g2)
    }
  }

  return(dat)
}
