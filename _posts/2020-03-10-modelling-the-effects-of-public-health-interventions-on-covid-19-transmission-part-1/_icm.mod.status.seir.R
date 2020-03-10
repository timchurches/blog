infection.icm <- function(dat, at) {

  ## Expected acts
  if (dat$param$groups == 1) {
    acts <- round(dat$param$act.rate * dat$epi$num[at - 1] / 2)
  }
  if (dat$param$groups == 2) {
    if (dat$param$balance == "g1") {
      acts <- round(dat$param$act.rate *
                      (dat$epi$num[at - 1] + dat$epi$num.g2[at - 1]) / 2)
    }
    if (dat$param$balance == "g2") {
      acts <- round(dat$param$act.rate.g2 *
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


    ## Discordant edgelist
    del$p1.stat <- dat$attr$status[del$p1]
    del$p2.stat <- dat$attr$status[del$p2]
    serodis <- (del$p1.stat == "s" & del$p2.stat == "i") |
               (del$p1.stat == "i" & del$p2.stat == "s")
    del <- del[serodis == TRUE, ]


    ## Transmission on edgelist
    if (nrow(del) > 0) {
      if (dat$param$groups == 1) {
        del$tprob <- dat$param$inf.prob
      } else {
        del$tprob <- ifelse(del$p1.stat == "s", dat$param$inf.prob,
                                                dat$param$inf.prob.g2)
      }
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$tprob <- del$tprob * (1 - dat$param$inter.eff)
      }
      del$trans <- rbinom(nrow(del), 1, del$tprob)
      del <- del[del$trans == TRUE, ]
      if (nrow(del) > 0) {
        if (dat$param$groups == 1) {
          newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
          nInf <- length(newIds)
        }
        if (dat$param$groups == 2) {
          newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
          newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
          nInf <- length(newIdsg1)
          nInfg2 <- length(newIdsg2)
          newIds <- c(newIdsg1, newIdsg2)
        }
        dat$attr$status[newIds] <- "e"
        dat$attr$infTime[newIds] <- at
      } else {
        nInf <- nInfg2 <- 0
      }
    } else {
      nInf <- nInfg2 <- 0
    }
  } else {
    nInf <- nInfg2 <- 0
  }


  ## Output
  if (at == 2) {
    dat$epi$se.flow <- c(0, nInf)
  } else {
    dat$epi$se.flow[at] <- nInf
  }
  if (dat$param$groups == 2) {
    if (at == 2) {
      dat$epi$se.flow.g2 <- c(0, nInfg2)
    } else {
      dat$epi$se.flow.g2[at] <- nInfg2
    }
  }

  return(dat)

}


progress.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (!(dat$control$type %in% c("SIR", "SIS", "SEIR"))) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status

  groups <- dat$param$groups
  group <- dat$attr$group

  type <- dat$control$type
  recovState <- ifelse(type %in% c("SIR", "SEIR"), "r", "s")
  progState <- "i"
  
  # --- progress ----
  prog.rand <- dat$control$prog.rand
  prog.rate <- dat$param$prog.rate
  prog.rate.g2 <- dat$param$prog.rate.g2
 
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
      nProg <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      status[ssample(idsElig[gElig == 1], nProg)] <- progState
      if (groups == 2) {
        nProgG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        status[ssample(idsElig[gElig == 2], nProg)] <- progState
      }
    }
  }
  dat$attr$status <- status
  
  # ----- recover ------- 
  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.g2 <- dat$param$rec.rate.g2

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & status == "i")
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
  if (type == "SEIR") {
    outName_b <- "ei.flow"
    outName_b[2] <- paste0(outName_b, ".g2")
  }
  
  ## Summary statistics
  if (at == 2) {
    dat$epi[[outName_a[1]]] <- c(0, nRecov)
    if (type == "SEIR") {
      dat$epi[[outName_b[1]]] <- c(0, nProg) 
    }
  } else {
    dat$epi[[outName_a[1]]][at] <- nRecov
    if (type == "SEIR") {
      dat$epi[[outName_b[1]]][at] <- nProg 
    }
  }
  if (groups == 2) {
    if (at == 2) {
      dat$epi[[outName_a[2]]] <- c(0, nRecovG2)
      if (type == "SEIR") {
        dat$epi[[outName_b[2]]] <- c(0, nProgG2) 
      }
    } else {
      dat$epi[[outName_a[2]]][at] <- nRecovG2
      if (type == "SEIR") {
        dat$epi[[outName_b[2]]][at] <- nProgG2 
      }
    }
  }

  return(dat)
}

###############

init_status.icm <- function(dat) {

  # Variables ---------------------------------------------------------------
  type <- dat$control$type
  group <- dat$attr$group
  nGroups <- dat$param$groups

  nG1 <- sum(group == 1)
  nG2 <- sum(group == 2)

  e.num <- dat$init$e.num
  i.num <- dat$init$i.num
  r.num <- dat$init$r.num
  e.num.g2 <- dat$init$e.num.g2
  i.num.g2 <- dat$init$i.num.g2
  r.num.g2 <- dat$init$r.num.g2


  # Status ------------------------------------------------------------------
  status <- rep("s", nG1 + nG2)
  if (type == "SEIR") {
    status[sample(which(group == 1 & status == "s"), size = e.num)] <- "e"
    if (nGroups == 2) {
      status[sample(which(group == 2 & status == "s"), size = e.num.g2)] <- "e"
    }
  }
  status[sample(which(group == 1), size = i.num)] <- "i"
  if (nGroups == 2) {
    status[sample(which(group == 2), size = i.num.g2)] <- "i"
  }
  if (type %in% c("SIR", "SEIR")) {
    status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
    if (nGroups == 2) {
      status[sample(which(group == 2 & status == "s"), size = r.num.g2)] <- "r"
    }
  }
  dat$attr$status <- status

  # Infection Time ----------------------------------------------------------
  idsInf <- which(status == "i")
  infTime <- rep(NA, length(status))

  # If vital=TRUE, infTime is a uniform draw over the duration of infection
  if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
    infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate) + 2
  } else {
    if (dat$control$type == "SI" || dat$param$rec.rate == 0) {
      # infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                 length(idsInf), replace = TRUE)
    } else {
      if (nGroups == 1) {
        infTime[idsInf] <- ssample(1:(-round(1 / dat$param$rec.rate) + 2),
                                   length(idsInf), replace = TRUE)
      }
      if (nGroups == 2) {
        infG1 <- which(status == "i" & group == 1)
        infTime[infG1] <- ssample(1:(-round(1 / dat$param$rec.rate) + 2),
                                  length(infG1), replace = TRUE)
        infG2 <- which(status == "i" & group == 2)
        infTime[infG2] <- ssample(1:(-round(1 / dat$param$rec.rate.g2) + 2),
                                  length(infG2), replace = TRUE)
      }
    }
  }
  dat$attr$infTime <- infTime

  return(dat)
}
