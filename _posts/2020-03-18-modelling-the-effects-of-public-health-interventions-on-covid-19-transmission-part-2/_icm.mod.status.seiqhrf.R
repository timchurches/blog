infection.seiqhrf.icm <- function(dat, at) {
  type <- dat$control$type
  
  # the following checks need to be moved to control.icm in due course
  nsteps <- dat$control$nsteps

  act.rate.i <- dat$param$act.rate.i
  if (!(length(act.rate.i) == 1 || length(act.rate.i == nsteps))) {
    stop("Length of act.rate.i must be 1 or the value of nsteps")
  }
  act.rate.i.g2 <- dat$param$act.rate.i.g2
  if (!is.null(act.rate.i.g2) && 
      !(length(act.rate.i.g2) == 1 || length(act.rate.i.g2 == nsteps))) {
    stop("Length of act.rate.i.g2 must be 1 or the value of nsteps")
  }
  inf.prob.i <- dat$param$inf.prob.i
  if (!(length(inf.prob.i) == 1 || length(inf.prob.i == nsteps))) {
    stop("Length of inf.prob.i must be 1 or the value of nsteps")
  }
  inf.prob.i.g2 <- dat$param$inf.prob.i.g2
  if (!is.null(inf.prob.i.g2) &&
      !(length(inf.prob.i.g2) == 1 || length(inf.prob.i.g2 == nsteps))) {
    stop("Length of inf.prob.i.g2 must be 1 or the value of nsteps")
  }
  if (type %in% c("SEIQHR", "SEIQHRF")) {  
    quar.rate <- dat$param$quar.rate
    if (!(length(quar.rate) == 1 || length(quar.rate == nsteps))) {
      stop("Length of quar.rate must be 1 or the value of nsteps")
    }
    quar.rate.g2 <- dat$param$quar.rate.g2
    if (!is.null(quar.rate.g2) &&
        !(length(quar.rate.g2) == 1 || length(quar.rate.g2 == nsteps))) {
      stop("Length of quar.rate.g2 must be 1 or the value of nsteps")
    }
    disch.rate <- dat$param$disch.rate
    if (!(length(disch.rate) == 1 || length(disch.rate == nsteps))) {
      stop("Length of disch.rate must be 1 or the value of nsteps")
    }
    disch.rate.g2 <- dat$param$disch.rate.g2
    if (!is.null(disch.rate.g2) &&
        !(length(disch.rate.g2) == 1 || length(disch.rate.g2 == nsteps))) {
      stop("Length of disch.rate.g2 must be 1 or the value of nsteps")
    }
  }
  
  if (type %in% c("SEIQHR", "SEIQHRF")) {  
    act.rate.e <- dat$param$act.rate.e
    if (!(length(act.rate.e) == 1 || length(act.rate.e == nsteps))) {
      stop("Length of act.rate.e must be 1 or the value of nsteps")
    }
    act.rate.e.g2 <- dat$param$act.rate.e.g2
    if (!is.null(act.rate.e.g2) &&
        !(length(act.rate.e.g2) == 1 || length(act.rate.e.g2 == nsteps))) {
      stop("Length of act.rate.e.g2 must be 1 or the value of nsteps")
    }
    inf.prob.e <- dat$param$inf.prob.e
    if (!(length(inf.prob.e) == 1 || length(inf.prob.e == nsteps))) {
      stop("Length of inf.prob.e must be 1 or the value of nsteps")
    }
    inf.prob.e.g2 <- dat$param$inf.prob.e.g2
    if (!is.null(inf.prob.e.g2) &&
        !(length(inf.prob.e.g2) == 1 || length(inf.prob.e.g2 == nsteps))) {
      stop("Length of inf.prob.e.g2 must be 1 or the value of nsteps")
    }
    
    act.rate.q <- dat$param$act.rate.q
    if (!(length(act.rate.q) == 1 || length(act.rate.q == nsteps))) {
      stop("Length of act.rate.q must be 1 or the value of nsteps")
    }
    act.rate.q.g2 <- dat$param$act.rate.q.g2
    if (!is.null(act.rate.q.g2) &&
        !(length(act.rate.q.g2) == 1 || length(act.rate.q.g2 == nsteps))) {
      stop("Length of act.rate.q.g2 must be 1 or the value of nsteps")
    }
    inf.prob.q <- dat$param$inf.prob.q
    if (!(length(inf.prob.q) == 1 || length(inf.prob.q == nsteps))) {
      stop("Length of inf.prob.q must be 1 or the value of nsteps")
    }
    inf.prob.q.g2 <- dat$param$inf.prob.q.g2
    if (!is.null(inf.prob.q.g2) &&
        !(length(inf.prob.q.g2) == 1 || length(inf.prob.q.g2 == nsteps))) {
      stop("Length of inf.prob.q.g2 must be 1 or the value of nsteps")
    }
  }

# Transmission from infected  
  ## Expected acts
  if (dat$param$groups == 1) {
    if (length(act.rate.i) > 1) {
      acts <- round(act.rate.i[at - 1] * dat$epi$num[at - 1] / 2)
    } else {
      acts <- round(act.rate.i * dat$epi$num[at - 1] / 2)
    }
  }
  if (dat$param$groups == 2) {
    if (dat$param$balance == "g1") {
      if (length(act.rate.i) > 1) {
        acts <- round(act.rate.i[at - 1] *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      } else {
        acts <- round(act.rate.i *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      }
    }
    if (dat$param$balance == "g2") {
      if (length(act.rate.i.g2) > 1) {
        acts <- round(act.rate.i.g2[at - 1] *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      } else {
        acts <- round(act.rate.i.g2 *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
      }
    }
  }

  ## Edgelist
  if (dat$param$groups == 1) {
    p1 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
  } else {
    p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1 & dat$attr$status != "f"),
                  acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2 & dat$attr$status != "f"),
                  acts, replace = TRUE)
  }

  del <- NULL
  if (length(p1) > 0 & length(p2) > 0) {
    del <- data.frame(p1, p2)
    if (dat$param$groups == 1) {
      while (any(del$p1 == del$p2)) {
        del$p2 <- ifelse(del$p1 == del$p2,
                         ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), 1), del$p2)
      }
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
      if (length(inf.prob.i) > 1) {
        del$tprob <- inf.prob.i[at]
      } else {
        del$tprob <- inf.prob.i
      }
    } else {
      if (length(inf.prob.i) > 1) {
        del$tprob <- ifelse(del$p1.stat == "s", inf.prob.i[at],
                                                inf.prob.i.g2[at])
      } else {
        del$tprob <- ifelse(del$p1.stat == "s", inf.prob.i,
                                                inf.prob.i.g2)
      }
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

  if (type %in% c("SEIQHRF")) {  
    
  # Transmission from exposed  
    ## Expected acts
    if (dat$param$groups == 1) {
      if (length(act.rate.e) > 1) {
        acts <- round(act.rate.e[at - 1] * dat$epi$num[at - 1] / 2)
      } else {
        acts <- round(act.rate.e * dat$epi$num[at - 1] / 2)
      }
    }
    if (dat$param$groups == 2) {
      if (dat$param$balance == "g1") {
        if (length(act.rate.e.g2) > 1) {
          acts <- round(act.rate.e.g2[at - 1] *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        } else {
          acts <- round(act.rate.e.g2 *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        }
      }
      if (dat$param$balance == "g2") {
        if (length(act.rate.e.g2) > 1) {
          acts <- round(act.rate.e.g2[at - 1] *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        } else {
          acts <- round(act.rate.e.g2 *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        }
      }
    }
  
    ## Edgelist
    if (dat$param$groups == 1) {
      p1 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
    } else {
      p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1 & dat$attr$status != "f"),
                    acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2 & dat$attr$status != "f"),
                    acts, replace = TRUE)
    }
  
    del <- NULL
    if (length(p1) > 0 & length(p2) > 0) {
      del <- data.frame(p1, p2)
      if (dat$param$groups == 1) {
        while (any(del$p1 == del$p2)) {
          del$p2 <- ifelse(del$p1 == del$p2,
                           ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), 1), del$p2)
        }
      }
  
      ## Discordant edgelist (del)
      del$p1.stat <- dat$attr$status[del$p1]
      del$p2.stat <- dat$attr$status[del$p2]
      # serodiscordance
      serodis <- (del$p1.stat == "s" & del$p2.stat == "e") |
                 (del$p1.stat == "e" & del$p2.stat == "s")
      del <- del[serodis == TRUE, ]
  
      ## Transmission on edgelist
      if (nrow(del) > 0) {
        if (dat$param$groups == 1) {
          if (length(inf.prob.e) > 1) {
            del$tprob <- inf.prob.e[at]
          } else {
            del$tprob <- inf.prob.e
          }
        } else {
          if (length(inf.prob.e) > 1) {
            del$tprob <- ifelse(del$p1.stat == "s", inf.prob.e[at],
                                                    inf.prob.e.g2[at])
          } else {
            del$tprob <- ifelse(del$p1.stat == "s", inf.prob.e,
                                                    inf.prob.e.g2)
          }
        }
        if (!is.null(dat$param$inter.eff.e) && at >= dat$param$inter.start.e &&
            at <= dat$param$inter.stop.e) {
          del$tprob <- del$tprob * (1 - dat$param$inter.eff.e)
        }
        del$trans <- rbinom(nrow(del), 1, del$tprob)
        del <- del[del$trans == TRUE, ]
        if (nrow(del) > 0) {
          if (dat$param$groups == 1) {
            newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
            nExp.e <- length(newIds)
          }
          if (dat$param$groups == 2) {
            newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
            newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
            nExp.e <- length(newIdsg1)
            nExpg2.e <- length(newIdsg2)
            newIds <- c(newIdsg1, newIdsg2)
          }
          dat$attr$status[newIds] <- "e"
          dat$attr$expTime[newIds] <- at
        } else {
          nExp.e <- nExpg2.e <- 0
        }
      } else {
        nExp.e <- nExpg2.e <- 0
      }
    } else {
      nExp.e <- nExpg2.e <- 0
    }
    
    
  # Transmission from quarantined  
    ## Expected acts
    if (dat$param$groups == 1) {
      if (length(act.rate.q) > 1) {
        acts <- round(act.rate.q[at - 1] * dat$epi$num[at - 1] / 2)
      } else {
        acts <- round(act.rate.q * dat$epi$num[at - 1] / 2)
      }
    }
    if (dat$param$groups == 2) {
      if (dat$param$balance == "g1") {
        if (length(act.rate.q.g2) > 1) {
          acts <- round(act.rate.q.g2[at - 1] *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        } else {
          acts <- round(act.rate.q.g2 *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        }
      }
      if (dat$param$balance == "g2") {
        if (length(act.rate.q.g2) > 1) {
          acts <- round(act.rate.q.g2[at - 1] *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        } else {
          acts <- round(act.rate.q.g2 *
                        (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
        }
      }
    }
  
    ## Edgelist
    if (dat$param$groups == 1) {
      p1 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
    } else {
      p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1 & dat$attr$status != "f"),
                    acts, replace = TRUE)
      p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2 & dat$attr$status != "f"),
                    acts, replace = TRUE)
    }
  
    del <- NULL
    if (length(p1) > 0 & length(p2) > 0) {
      del <- data.frame(p1, p2)
      if (dat$param$groups == 1) {
        while (any(del$p1 == del$p2)) {
          del$p2 <- ifelse(del$p1 == del$p2,
                           ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), 1), del$p2)
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
          if (length(inf.prob.q) > 1) {
            del$tprob <- inf.prob.q[at]
          } else {
            del$tprob <- inf.prob.q
          }
        } else {
          if (length(inf.prob.q) > 1) {
            del$tprob <- ifelse(del$p1.stat == "s", inf.prob.q[at],
                                                    inf.prob.q.g2[at])
          } else {
            del$tprob <- ifelse(del$p1.stat == "s", inf.prob.q,
                                                    inf.prob.q.g2)
          }
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
  if (type %in% c("SEIQHR", "SEIQHRF")) {  
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

# utility functions
cum_discr_si <- function(vecTimeSinceExp, scale, shape) {
  vlen <- length(vecTimeSinceExp)
  if (vlen > 0) {
    probVec <- numeric(vlen)
    for (p in 1:vlen) {
      probVec[p] <- pweibull(vecTimeSinceExp[p], shape=shape, scale=scale)
    }
  } else {
      probVec <- 0    
  }
  return(probVec)
}

progress.seiqhrf.icm <- function(dat, at) {

  #print(at)
  #print(dat$control$type)
  #print("-------")
  
  # Conditions --------------------------------------------------------------
  if (!(dat$control$type %in% c("SIR", "SIS", "SEIR", "SEIQHR", "SEIQHRF"))) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status

  groups <- dat$param$groups
  group <- dat$attr$group

  type <- dat$control$type
  recovState <- ifelse(type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF"), "r", "s")
  progState <- "i"
  quarState <- "q"
  hospState <- "h"
  fatState <- "f"
  
  # --- progress from exposed to infectious ----
  prog.rand <- dat$control$prog.rand
  prog.rate <- dat$param$prog.rate
  prog.rate.g2 <- dat$param$prog.rate.g2
  prog.dist.scale <- dat$param$prog.dist.scale
  prog.dist.shape <- dat$param$prog.dist.shape
  prog.dist.scale.g2 <- dat$param$prog.dist.scale.g2
  prog.dist.shape.g2 <- dat$param$prog.dist.shape.g2
  
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
        dat$attr$infTime[idsProg] <- at
      }
    } else {
      vecTimeSinceExp <- at - dat$attr$expTime[idsElig]
      gammaRatesElig <- pweibull(vecTimeSinceExp, prog.dist.shape, scale=prog.dist.scale) 
      nProg <- round(sum(gammaRatesElig[gElig == 1], na.rm=TRUE))
      if (nProg > 0) {
        ids2bProg <- ssample(idsElig[gElig == 1], 
                      nProg, prob = gammaRatesElig[gElig == 1])
        status[ids2bProg] <- progState
        dat$attr$infTime[ids2bProg] <- at
        # debug
        if (FALSE & at <= 30) {
          print(paste("at:", at))
          print("idsElig:")
          print(idsElig[gElig == 1])
          print("vecTimeSinceExp:")
          print(vecTimeSinceExp[gElig == 1])
          print("gammaRatesElig:")
          print(gammaRatesElig)
          print(paste("nProg:",nProg))
          print(paste("sum of elig rates:", round(sum(gammaRatesElig[gElig == 1]))))
          print(paste("sum(gElig == 1):", sum(gElig == 1)))
          print("ids progressed:")
          print(ids2bProg)
          print("probs of ids to be progressed:")
          print(gammaRatesElig[which(idsElig %in% ids2bProg)]) 
          print("days since exposed of ids to be progressed:")
          print(vecTimeSinceExp[which(idsElig %in% ids2bProg)]) 
          print("------")
        }  
      }
      if (groups == 2) {
        nProgG2 <- round(sum(gammaRatesElig[gElig == 2], na.rm=TRUE))
        if (nProgG2 > 0) {
          ids2bProgG2 <- ssample(idsElig[gElig == 2], 
                        nProgG2, prob = gammaRatesElig[gElig == 2])
          status[ids2bProgG2] <- progState
          dat$attr$infTime[ids2bProgG2] <- at
        }
      }
    }
  }
  dat$attr$status <- status

  if (type %in% c("SEIQHR", "SEIQHRF")) {  
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

      if (length(quar.rate) > 1) {
          qrate <- quar.rate[at]
      } else {
          qrate <- quar.rate
      }
      if (length(quar.rate.g2) > 1) {
          qrate.g2 <- quar.rate.g2[at]
      } else {
          qrate.g2 <- quar.rate.g2
      }
      rates <- c(qrate, qrate.g2)
      ratesElig <- rates[gElig]
      if (quar.rand == TRUE) {
        vecQuar <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecQuar) > 0) {
          idsQuar <- idsElig[vecQuar]
          nQuar <- sum(group[idsQuar] == 1)
          nQuarG2 <- sum(group[idsQuar] == 2)
          status[idsQuar] <- quarState
          dat$attr$quarTime[idsQuar] <- at
        }
      } else {
        nQuar <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
        idsQuar <- ssample(idsElig[gElig == 1], nQuar)
        status[idsQuar] <- quarState
        dat$attr$quarTime[idsQuar] <- at
        if (groups == 2) {
          nQuarG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          idsQuarG2 <- ssample(idsElig[gElig == 2], nQuarG2)
          status[idsQuarG2] <- quarState
          dat$attr$quarTime[idsQuarG2] <- at
        }
      }
    }
    dat$attr$status <- status
  
    # ----- need to be hospitalised ------- 
    hosp.rand <- dat$control$hosp.rand
    hosp.rate <- dat$param$hosp.rate
    hosp.rate.g2 <- dat$param$hosp.rate.g2
  
    nHosp <- nHospG2 <- 0
    idsElig <- which(active == 1 & (status == "i" | status == "q"))
    nElig <- length(idsElig)
    idsHosp <- numeric(0)
    
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
        idsHosp <- ssample(idsElig[gElig == 1], nHosp)
        status[idsHosp] <- hospState
        if (groups == 2) {
          nHospG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          idsHospG2 <- ssample(idsElig[gElig == 2], nHospG2)
          status[idsHospG2] <- hospState
          idsHosp <- c(idsHosp, idsHospG2)
        }
      }
    }
    dat$attr$status <- status
    dat$attr$hospTime[idsHosp] <- at

    # ----- discharge from need to be hospitalised ------- 
    disch.rand <- dat$control$disch.rand
    disch.rate <- dat$param$disch.rate
    disch.rate.g2 <- dat$param$disch.rate.g2
  
    nDisch <- nDischG2 <- 0
    idsElig <- which(active == 1 & status == "h")
    nElig <- length(idsElig)
    idsDisch <- numeric(0)
  
    if (nElig > 0) {
  
      gElig <- group[idsElig]
      rates <- c(disch.rate, disch.rate.g2)

      if (length(disch.rate) > 1) {
          dcrate <- disch.rate[at]
      } else {
          dcrate <- disch.rate
      }
      if (length(disch.rate.g2) > 1) {
          dcrate.g2 <- disch.rate.g2[at]
      } else {
          dcrate.g2 <- disch.rate.g2
      }
      
      rates <- c(dcrate, dcrate.g2)
      ratesElig <- rates[gElig]
  
      if (disch.rand == TRUE) {
        vecDisch <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecDisch) > 0) {
          idsDisch <- idsElig[vecDisch]
          nDisch <- sum(group[idsDisch] == 1)
          nDischG2 <- sum(group[idsDisch] == 2)
          status[idsDisch] <- recovState
        }
      } else {
        nDisch <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
        idsDisch <- ssample(idsElig[gElig == 1], nDisch)
        status[idsDisch] <- recovState
        if (groups == 2) {
          nDischG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          idsDischG2 <- ssample(idsElig[gElig == 2], nDischG2)
          status[idsDischG2] <- recovState
          idsDisch <- c(idsDisch, idsDischG2)
        }
      }
    }
    dat$attr$status <- status
    dat$attr$dischTime[idsDisch] <- at
  }
  
  # ----- recover ------- 
  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.g2 <- dat$param$rec.rate.g2
  rec.dist.scale <- dat$param$rec.dist.scale
  rec.dist.shape <- dat$param$rec.dist.shape
  rec.dist.scale.g2 <- dat$param$rec.dist.scale.g2
  rec.dist.shape.g2 <- dat$param$rec.dist.shape.g2

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & (status == "i" | status == "q" | status == "h"))
  nElig <- length(idsElig)
  idsRecov <- numeric(0)
  
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
      vecTimeSinceExp <- at - dat$attr$expTime[idsElig]
      vecTimeSinceExp[is.na(vecTimeSinceExp)] <- 0
      gammaRatesElig <- pweibull(vecTimeSinceExp, rec.dist.shape, scale=rec.dist.scale) 
      nRecov <- round(sum(gammaRatesElig[gElig == 1], na.rm=TRUE))
      if (nRecov > 0) {
        idsRecov <- ssample(idsElig[gElig == 1], 
                      nRecov, prob = gammaRatesElig[gElig == 1])
        status[idsRecov] <- recovState
        # debug
        if (FALSE & at <= 30) {
          print(paste("at:", at))
          print("idsElig:")
          print(idsElig[gElig == 1])
          print("vecTimeSinceExp:")
          print(vecTimeSinceExp[gElig == 1])
          print("gammaRatesElig:")
          print(gammaRatesElig)
          print(paste("nRecov:",nRecov))
          print(paste("sum of elig rates:", round(sum(gammaRatesElig[gElig == 1]))))
          print(paste("sum(gElig == 1):", sum(gElig == 1)))
          print("ids recovered:")
          print(idsRecov)
          print("probs of ids to be progressed:")
          print(gammaRatesElig[which(idsElig %in% idsRecov)]) 
          print("days since exposed of ids to be Recovered:")
          print(vecTimeSinceExp[which(idsElig %in% idsRecov)]) 
          print("------")
        }  

      }
      if (groups == 2) {
        nRecovG2 <- round(sum(gammaRatesElig[gElig == 2], na.rm=TRUE))
        if (nRecovG2 > 0) {
          idsRecovG2 <- ssample(idsElig[gElig == 2], 
                        nRecovG2, prob = gammaRatesElig[gElig == 2])
          status[idsRecovG2] <- recovState
          idsRecov <- c(idsRecov, idsRecovG2)
        }
      }
    }
  }
  dat$attr$status <- status
  dat$attr$recovTime[idsRecov] <- at
  
  fatEnable <- TRUE
  if (fatEnable & type %in% c("SEIQHRF")) {  
    # ----- case fatality ------- 
    fat.rand <- dat$control$fat.rand
    fat.rate.base <- dat$param$fat.rate.base
    fat.rate.base.g2 <- dat$param$fat.rate.base.g2
    fat.rate.base.g2 <- ifelse(is.null(fat.rate.base.g2), 
                               0, fat.rate.base.g2)
    fat.rate.overcap <- dat$param$fat.rate.overcap
    fat.rate.overcap.g2 <- dat$param$fat.rate.overcap.g2
    fat.rate.overcap.g2 <- ifelse(is.null(fat.rate.overcap.g2), 
                                  0, fat.rate.overcap.g2)
    hosp.cap <- dat$param$hosp.cap
    fat.tcoeff <- dat$param$fat.tcoeff
    
    nFat <- nFatG2 <- 0
    idsElig <- which(active == 1 & status =="h")
    nElig <- length(idsElig)
  
    if (nElig > 0) {
      gElig <- group[idsElig]
      timeInHospElig <- at - dat$attr$hospTime[idsElig]
      rates <- c(fat.rate.base, fat.rate.base.g2)
      h.num.yesterday <- 0
      if (!is.null(dat$epi$h.num[at - 1])) {
        h.num.yesterday <- dat$epi$h.num[at - 1]
        if (h.num.yesterday > hosp.cap) {
          blended.rate <- ((hosp.cap * fat.rate.base) + 
                ((h.num.yesterday - hosp.cap) * fat.rate.overcap)) / 
                  h.num.yesterday
          blended.rate.g2 <- ((hosp.cap * fat.rate.base.g2) + 
                ((h.num.yesterday - hosp.cap) * fat.rate.overcap.g2)) / 
                  h.num.yesterday
          rates <- c(blended.rate, blended.rate.g2)
        }  
      } 
      ratesElig <- rates[gElig]
      ratesElig <- ratesElig + timeInHospElig*fat.tcoeff*ratesElig

      if (fat.rand == TRUE) {
        vecFat <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecFat) > 0) {
          idsFat <- idsElig[vecFat]
          nFat <- sum(group[idsFat] == 1)
          nFatG2 <- sum(group[idsFat] == 2)
          status[idsFat] <- fatState
          dat$attr$fatTime[idsFat] <- at
        }
      } else {
        nFat <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
        idsFat <- ssample(idsElig[gElig == 1], nFat)
        status[idsFat] <- fatState
        dat$attr$fatTime[idsFat] <- at
        if (groups == 2) {
          nFatG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
          idsFatG2 <- ssample(idsElig[gElig == 2], nFatG2)
          status[idsFatG2] <- fatState
          dat$attr$fatTime[idsFatG2] <- at
        }
      }
    }
    dat$attr$status <- status
  }

  # Output ------------------------------------------------------------------
  outName_a <- ifelse(type %in% c("SIR", "SEIR"), "ir.flow", "is.flow")
  outName_a[2] <- paste0(outName_a, ".g2")
  if (type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
    outName_b <- "ei.flow"
    outName_b[2] <- paste0(outName_b, ".g2")
  }
  if (type %in% c("SEIQHR", "SEIQHRF")) {
    outName_c <- "iq.flow"
    outName_c[2] <- paste0(outName_c, ".g2")
    outName_d <- "iq2h.flow"
    outName_d[2] <- paste0(outName_d, ".g2")
  }
  if (type %in% c("SEIQHRF")) {
    outName_e <- "hf.flow"
    outName_e[2] <- paste0(outName_e, ".g2")
  }
  ## Summary statistics
  if (at == 2) {
    dat$epi[[outName_a[1]]] <- c(0, nRecov)
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]] <- c(0, nProg) 
    }
    if (type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi[[outName_c[1]]] <- c(0, nQuar) 
      dat$epi[[outName_d[1]]] <- c(0, nHosp) 
    }
    if (fatEnable & type %in% c("SEIQHRF")) {
      dat$epi[[outName_e[1]]] <- c(0, nFat) 
    }
  } else {
    dat$epi[[outName_a[1]]][at] <- nRecov
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]][at] <- nProg 
    }
    if (type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi[[outName_c[1]]][at] <- nQuar 
      dat$epi[[outName_d[1]]][at] <- nHosp 
    }
    if (fatEnable & type %in% c("SEIQHRF")) {
      dat$epi[[outName_e[1]]][at] <- nFat 
    }
  }
  if (groups == 2) {
    if (at == 2) {
      dat$epi[[outName_a[2]]] <- c(0, nRecovG2)
      if (type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi[[outName_b[2]]] <- c(0, nProgG2) 
      }
      if (type %in% c("SEIQHR", "SEIQHRF")) {
        dat$epi[[outName_c[2]]] <- c(0, nQuarG2) 
        dat$epi[[outName_d[2]]] <- c(0, nHospG2) 
      }
      if (type %in% c("SEIQHRF")) {
        dat$epi[[outName_e[2]]] <- c(0, nFatG2) 
      }
    } else {
      dat$epi[[outName_a[2]]][at] <- nRecovG2
      if (type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi[[outName_b[2]]][at] <- nProgG2 
      }
      if (type %in% c("SEIQHR", "SEIQHRF")) {
        dat$epi[[outName_c[2]]][at] <- nQuarG2 
        dat$epi[[outName_d[2]]][at] <- nHospG2 
      }
      if (type %in% c("SEIQHRF")) {
        dat$epi[[outName_e[2]]][at] <- nFatG2 
      }
    }
  }

  return(dat)
}

###############

