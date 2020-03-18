infection.seiqhr.icm <- function(dat, at) {
  type <- dat$control$type

# Transmission from infected  
  ## Expected acts
  if (dat$param$groups == 1) {
    acts <- round(dat$param$act.rate.i * dat$epi$num[at - 1] / 2)
  }
  if (dat$param$groups == 2) {
    if (dat$param$balance == "g1") {
      acts <- round(dat$param$act.rate.i *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
    }
    if (dat$param$balance == "g2") {
      acts <- round(dat$param$act.rate.i.g2 *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
    }
  }

  ## Edgelist
  if (dat$param$groups == 1) {
    p1 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
  } else {
    p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1),
                  acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2),
                  acts, replace = TRUE)
  }

  del <- NULL
  if (length(p1) > 0 & length(p2) > 0) {
    del <- data.frame(p1, p2)
    if (dat$param$groups == 1) {
      while (any(del$p1 == del$p2)) {
        del$p2 <- ifelse(del$p1 == del$p2,
                         ssample(which(dat$attr$active == 1), 1), del$p2)
      }
    }

    ## Discordant edgelist (del)
    del$p1.stat <- dat$attr$status[del$p1]
    del$p2.stat <- dat$attr$status[del$p2]
    serodis <- (del$p1.stat == "s" & del$p2.stat == "i") |
               (del$p1.stat == "i" & del$p2.stat == "s")
    del <- del[serodis == TRUE, ]

    ## Transmission on edgelist
    if (nrow(del) > 0) {
      if (dat$param$groups == 1) {
        del$tprob <- dat$param$inf.prob.i
      } else {
        del$tprob <- ifelse(del$p1.stat == "s", dat$param$inf.prob.i,
                                                dat$param$inf.prob.i.g2)
      }
      if (!is.null(dat$param$inter.eff.i) && at >= dat$param$inter.start.i &&
          at <= dat$param$inter.stop.i) {
        del$tprob <- del$tprob * (1 - dat$param$inter.eff.i)
      }
      del$trans <- rbinom(nrow(del), 1, del$tprob)
      del <- del[del$trans == TRUE, ]
      if (nrow(del) > 0) {
        if (dat$param$groups == 1) {
          newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
          nExp.i <- length(newIds)
        }
        if (dat$param$groups == 2) {
          newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
          newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
          nExp.i <- length(newIdsg1)
          nExpg2.i <- length(newIdsg2)
          newIds <- c(newIdsg1, newIdsg2)
        }
        dat$attr$status[newIds] <- "e"
        dat$attr$expTime[newIds] <- at
      } else {
        nExp.i <- nExpg2.i <- 0
      }
    } else {
      nExp.i <- nExpg2.i <- 0
    }
  } else {
    nExp.i <- nExpg2.i <- 0
  }

  if (type == "SEIQHR") {  
  # Transmission from quarantined  
    ## Expected acts
    if (dat$param$groups == 1) {
      acts <- round(dat$param$act.rate.q * dat$epi$num[at - 1] / 2)
    }
    if (dat$param$groups == 2) {
      if (dat$param$balance == "g1") {
        acts <- round(dat$param$act.rate.q *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      }
      if (dat$param$balance == "g2") {
        acts <- round(dat$param$act.rate.q.g2 *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      }
    }
  
    ## Edgelist
    if (dat$param$groups == 1) {
      p1 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
    } else {
      p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1),
                    acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2),
                    acts, replace = TRUE)
    }
  
    del <- NULL
    if (length(p1) > 0 & length(p2) > 0) {
      del <- data.frame(p1, p2)
      if (dat$param$groups == 1) {
        while (any(del$p1 == del$p2)) {
          del$p2 <- ifelse(del$p1 == del$p2,
                           ssample(which(dat$attr$active == 1), 1), del$p2)
        }
      }
  
      ## Discordant edgelist (del)
      del$p1.stat <- dat$attr$status[del$p1]
      del$p2.stat <- dat$attr$status[del$p2]
      # serodiscordance
      serodis <- (del$p1.stat == "s" & del$p2.stat == "q") |
                 (del$p1.stat == "q" & del$p2.stat == "s")
      del <- del[serodis == TRUE, ]
  
      ## Transmission on edgelist
      if (nrow(del) > 0) {
        if (dat$param$groups == 1) {
          del$tprob <- dat$param$inf.prob.q
        } else {
          del$tprob <- ifelse(del$p1.stat == "s", dat$param$inf.prob.q,
                                                  dat$param$inf.prob.q.g2)
        }
        if (!is.null(dat$param$inter.eff.q) && at >= dat$param$inter.start.q &&
            at <= dat$param$inter.stop.q) {
          del$tprob <- del$tprob * (1 - dat$param$inter.eff.q)
        }
        del$trans <- rbinom(nrow(del), 1, del$tprob)
        del <- del[del$trans == TRUE, ]
        if (nrow(del) > 0) {
          if (dat$param$groups == 1) {
            newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
            nExp.q <- length(newIds)
          }
          if (dat$param$groups == 2) {
            newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
            newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
            nExp.q <- length(newIdsg1)
            nExpg2.q <- length(newIdsg2)
            newIds <- c(newIdsg1, newIdsg2)
          }
          dat$attr$status[newIds] <- "e"
          dat$attr$expTime[newIds] <- at
        } else {
          nExp.q <- nExpg2.q <- 0
        }
      } else {
        nExp.q <- nExpg2.q <- 0
      }
    } else {
      nExp.q <- nExpg2.q <- 0
    }
  }

  ## Output
  if (type == "SEIQHR") {  
    if (at == 2) {
      dat$epi$se.flow <- c(0, nExp.i + nExp.q)
    } else {
      dat$epi$se.flow[at] <- nExp.i + nExp.q
    }
    if (dat$param$groups == 2) {
      if (at == 2) {
        dat$epi$se.flow.g2 <- c(0, nExpg2.i + nExpg2.q )
      } else {
        dat$epi$se.flow.g2[at] <- nExpg2.i + nExpg2.q
      }
    }
  } else {
    if (at == 2) {
      dat$epi$se.flow <- c(0, nExp.i)
    } else {
      dat$epi$se.flow[at] <- nExp.i
    }
    if (dat$param$groups == 2) {
      if (at == 2) {
        dat$epi$se.flow.g2 <- c(0, nExpg2.i)
      } else {
        dat$epi$se.flow.g2[at] <- nExpg2.i
      }
    }
  }
  return(dat)

}


progress.seiqhr.icm <- function(dat, at) {

  #print(at)
  #print(dat$control$type)
  #print("-------")
  
  # Conditions --------------------------------------------------------------
  if (!(dat$control$type %in% c("SIR", "SIS", "SEIR", "SEIQHR"))) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status

  groups <- dat$param$groups
  group <- dat$attr$group

  type <- dat$control$type
  recovState <- ifelse(type %in% c("SIR", "SEIR", "SEIQHR"), "r", "s")
  progState <- "i"
  quarState <- "q"
  hospState <- "h"
  
  # --- progress ----
  prog.rand <- dat$control$prog.rand
  prog.rate <- dat$param$prog.rate
  prog.rate.g2 <- dat$param$prog.rate.g2
  prog.dist.mu <- dat$param$prog.dist.mu
  prog.dist.sigma <- dat$param$prog.dist.sigma
  prog.dist.mu.g2 <- dat$param$prog.dist.mu.g2
  prog.dist.sigma.g2 <- dat$param$prog.dist.sigma.g2
  
  nProg <- nProgG2 <- 0
  idsElig <- which(active == 1 & status == "e")
  nElig <- length(idsElig)

  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(prog.rate, prog.rate.g2)
    ratesElig <- rates[gElig]

    if (prog.rand == TRUE) {
      vecProg <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecProg) > 0) {
        idsProg <- idsElig[vecProg]
        nProg <- sum(group[idsProg] == 1)
        nProgG2 <- sum(group[idsProg] == 2)
        status[idsProg] <- progState
      }
    } else {
      vecTimeSinceExp <- at - dat$attr$expTime[idsElig]
      gammaRatesElig <- EpiEstim::discr_si(vecTimeSinceExp, prog.dist.mu, prog.dist.sigma) 
      nProg <- sum(gammaRatesElig[gElig == 1] > 0.0)
      if (nProg > 0) {
        if (FALSE & at <= 30) {
          print(idsElig[gElig == 1])
          print(vecTimeSinceExp[gElig == 1])
          print(nProg)
          print(round(sum(gammaRatesElig[gElig == 1])))
          print(sum(gElig == 1))
          print(ssample(idsElig[gElig == 1], 
                     nProg, prob = gammaRatesElig[gElig == 1]))
          print("------")
        }  
        status[ssample(idsElig[gElig == 1], 
                     nProg, prob = gammaRatesElig[gElig == 1])] <- progState
      }
      if (groups == 2) {
        nProgG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        if (nProgG2 > 0) {
            status[ssample(idsElig[gElig == 2], 
                       nProgG2, prob = gammaRatesElig[gElig == 2])] <- progState
        }
      }
    }
  }
  dat$attr$status <- status
  dat$attr$infTime[idsProg] <- at

  if (type == "SEIQHR") {  
    # ----- quarantine ------- 
    quar.rand <- dat$control$quar.rand
    quar.rate <- dat$param$quar.rate
    quar.rate.g2 <- dat$param$quar.rate.g2
  
    nQuar <- nQuarG2 <- 0
    idsElig <- which(active == 1 & status == "i")
    nElig <- length(idsElig)
  
    if (nElig > 0) {
  
      gElig <- group[idsElig]
      rates <- c(quar.rate, quar.rate.g2)
      ratesElig <- rates[gElig]
  
      if (quar.rand == TRUE) {
        vecQuar <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecQuar) > 0) {
          idsQuar <- idsElig[vecQuar]
          nQuar <- sum(group[idsQuar] == 1)
          nQuarG2 <- sum(group[idsQuar] == 2)
          status[idsQuar] <- quarState
        }
      } else {
        nQuar <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
        status[ssample(idsElig[gElig == 1], nQuar)] <- quarState
        if (groups == 2) {
          nQuarG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          status[ssample(idsElig[gElig == 2], nQuar)] <- quarState
        }
      }
    }
    dat$attr$status <- status
  
    # ----- hospitalise ------- 
    hosp.rand <- dat$control$hosp.rand
    hosp.rate <- dat$param$hosp.rate
    hosp.rate.g2 <- dat$param$hosp.rate.g2
  
    nHosp <- nHospG2 <- 0
    idsElig <- which(active == 1 & (status == "i" | status == "q"))
    nElig <- length(idsElig)
  
    if (nElig > 0) {
  
      gElig <- group[idsElig]
      rates <- c(hosp.rate, hosp.rate.g2)
      ratesElig <- rates[gElig]
  
      if (hosp.rand == TRUE) {
        vecHosp <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecHosp) > 0) {
          idsHosp <- idsElig[vecHosp]
          nHosp <- sum(group[idsHosp] == 1)
          nHospG2 <- sum(group[idsHosp] == 2)
          status[idsHosp] <- hospState
        }
      } else {
        nHosp <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
        status[ssample(idsElig[gElig == 1], nHosp)] <- hospState
        if (groups == 2) {
          nHospG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          status[ssample(idsElig[gElig == 2], nHosp)] <- hospState
        }
      }
    }
    dat$attr$status <- status
  }
  
  # ----- recover ------- 
  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.g2 <- dat$param$rec.rate.g2

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & (status == "i" | status == "q" | status == "h"))
  nElig <- length(idsElig)

  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(rec.rate, rec.rate.g2)
    ratesElig <- rates[gElig]

    if (rec.rand == TRUE) {
      vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecRecov) > 0) {
        idsRecov <- idsElig[vecRecov]
        nRecov <- sum(group[idsRecov] == 1)
        nRecovG2 <- sum(group[idsRecov] == 2)
        status[idsRecov] <- recovState
      }
    } else {
      nRecov <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      status[ssample(idsElig[gElig == 1], nRecov)] <- recovState
      if (groups == 2) {
        nRecovG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        status[ssample(idsElig[gElig == 2], nRecov)] <- recovState
      }
    }
  }
  dat$attr$status <- status


  # Output ------------------------------------------------------------------
  outName_a <- ifelse(type %in% c("SIR", "SEIR"), "ir.flow", "is.flow")
  outName_a[2] <- paste0(outName_a, ".g2")
  if (type %in% c("SEIR", "SEIQHR")) {
    outName_b <- "ei.flow"
    outName_b[2] <- paste0(outName_b, ".g2")
  }
  if (type == "SEIQHR") {
    outName_c <- "iq.flow"
    outName_c[2] <- paste0(outName_c, ".g2")
    outName_d <- "iq2h.flow"
    outName_d[2] <- paste0(outName_d, ".g2")
  }
  
  ## Summary statistics
  if (at == 2) {
    dat$epi[[outName_a[1]]] <- c(0, nRecov)
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]] <- c(0, nProg) 
    }
    if (type == "SEIQHR") {
      dat$epi[[outName_c[1]]] <- c(0, nQuar) 
      dat$epi[[outName_d[1]]] <- c(0, nHosp) 
    }
  } else {
    dat$epi[[outName_a[1]]][at] <- nRecov
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]][at] <- nProg 
    }
    if (type == "SEIQHR") {
      dat$epi[[outName_c[1]]][at] <- nQuar 
      dat$epi[[outName_d[1]]][at] <- nHosp 
    }
  }
  if (groups == 2) {
    if (at == 2) {
      dat$epi[[outName_a[2]]] <- c(0, nRecovG2)
      if (type %in% c("SEIR", "SEIQHR")) {
        dat$epi[[outName_b[2]]] <- c(0, nProgG2) 
      }
      if (type == "SEIQHR") {
        dat$epi[[outName_c[2]]] <- c(0, nQuarG2) 
        dat$epi[[outName_d[2]]] <- c(0, nHospG2) 
      }
    } else {
      dat$epi[[outName_a[2]]][at] <- nRecovG2
      if (type %in% c("SEIR", "SEIQHR")) {
        dat$epi[[outName_b[2]]][at] <- nProgG2 
      }
    if (type == "SEIQHR") {
      dat$epi[[outName_c[2]]][at] <- nQuarG2 
      dat$epi[[outName_d[2]]][at] <- nHospG2 
    }
    }
  }

  return(dat)
}

###############

