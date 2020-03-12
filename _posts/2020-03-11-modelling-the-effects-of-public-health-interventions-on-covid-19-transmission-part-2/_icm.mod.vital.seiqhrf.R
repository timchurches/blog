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
