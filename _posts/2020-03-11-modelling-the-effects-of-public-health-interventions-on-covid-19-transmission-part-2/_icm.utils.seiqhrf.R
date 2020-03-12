get_prev.seiqhrf.icm <- function(dat, at) {

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(dat$attr$active == 1 &
                         dat$attr$status == "s" &
                         dat$attr$group == 1)
    dat$epi$i.num <- sum(dat$attr$active == 1 &
                         dat$attr$status == "i" &
                         dat$attr$group == 1)
    if (dat$control$type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi$e.num <- sum(dat$attr$active == 1 &
                             dat$attr$status == "e" &
                             dat$attr$group == 1)
    }        
    if (dat$control$type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF")) {
      dat$epi$r.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "r" &
                           dat$attr$group == 1)
    }
    if (dat$control$type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi$q.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "q" &
                           dat$attr$group == 1)
      dat$epi$h.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "h" &
                           dat$attr$group == 1)
    }
    if (dat$control$type =="SEIQHRF") {
      dat$epi$f.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "f" &
                           dat$attr$group == 1)
    }
    if (dat$control$type == "SIR") {
      dat$epi$num <- dat$epi$s.num +
                     dat$epi$i.num +
                     dat$epi$r.num
    } else if (dat$control$type == "SEIR") {
      dat$epi$num <- dat$epi$s.num +
                     dat$epi$e.num +
                     dat$epi$i.num +
                     dat$epi$r.num
    } else if (dat$control$type == "SEIQHR") {
      dat$epi$num <- dat$epi$s.num +
                     dat$epi$e.num +
                     dat$epi$i.num +
                     dat$epi$q.num +
                     dat$epi$h.num +
                     dat$epi$r.num
    } else if (dat$control$type == "SEIQHRF") {
      dat$epi$num <- dat$epi$s.num +
                     dat$epi$e.num +
                     dat$epi$i.num +
                     dat$epi$q.num +
                     dat$epi$h.num +
                     dat$epi$r.num +
                     dat$epi$f.num
    } else {
      dat$epi$num <- dat$epi$s.num + dat$epi$i.num
    }         
    if (dat$param$groups == 2) {
      dat$epi$s.num.g2 <- sum(dat$attr$active == 1 &
                              dat$attr$status == "s" &
                              dat$attr$group == 2)
      if (dat$control$type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
          dat$epi$e.num.g2 <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "e" &
                                  dat$attr$group == 2)
      }
      if (dat$control$type %in% c("SEIQHR", "SEIQHRF")) {
          dat$epi$q.num.g2 <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "q" &
                                  dat$attr$group == 2)
          dat$epi$h.num.g2 <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "h" &
                                  dat$attr$group == 2)
      }
      if (dat$control$type %in% c("SEIQHRF")) {
          dat$epi$f.num.g2 <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "f" &
                                  dat$attr$group == 2)
      }
      dat$epi$i.num.g2 <- sum(dat$attr$active == 1 &
                              dat$attr$status == "i" &
                              dat$attr$group == 2)
      dat$epi$num.g2 <- dat$epi$s.num.g2 + dat$epi$i.num.g2
      if (dat$control$type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi$r.num.g2 <- sum(dat$attr$active == 1 &
                                dat$attr$status == "r" &
                                dat$attr$group == 2)
      }
      if (dat$control$type == "SIR") {
        dat$epi$num.g2 <- dat$epi$s.num.g2 +
                          dat$epi$i.num.g2 +
                          dat$epi$r.num.g2
      } else if (dat$control$type == "SEIR") {
        dat$epi$num.g2 <- dat$epi$s.num.g2 +
                          dat$epi$e.num.g2 +
                          dat$epi$i.num.g2 +
                          dat$epi$r.num.g2
      } else if (dat$control$type == "SEIQHR") {
        dat$epi$num.g2 <- dat$epi$s.num.g2 +
                          dat$epi$e.num.g2 +
                          dat$epi$i.num.g2 +
                          dat$epi$q.num.g2 +
                          dat$epi$h.num.g2 +
                          dat$epi$r.num.g2
      } else if (dat$control$type == "SEIQHRF") {
        dat$epi$num.g2 <- dat$epi$s.num.g2 +
                          dat$epi$e.num.g2 +
                          dat$epi$i.num.g2 +
                          dat$epi$q.num.g2 +
                          dat$epi$h.num.g2 +
                          dat$epi$r.num.g2 +
                          dat$epi$f.num.g2
      } else {
        dat$epi$num.g2 <- dat$epi$s.num.g2 + dat$epi$i.num.g2
      }         
    }
  } else {
    dat$epi$s.num[at] <- sum(dat$attr$active == 1 &
                             dat$attr$status == "s" &
                             dat$attr$group == 1)
    if (dat$control$type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi$e.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$status == "e" &
                                 dat$attr$group == 1)
    }
    dat$epi$i.num[at] <- sum(dat$attr$active == 1 &
                             dat$attr$status == "i" &
                             dat$attr$group == 1)
    if (dat$control$type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF")) {
      dat$epi$r.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "r" &
                               dat$attr$group == 1)
    }
    if (dat$control$type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi$q.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "q" &
                               dat$attr$group == 1)
      dat$epi$h.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "h" &
                               dat$attr$group == 1)
    }
    if (dat$control$type %in% c("SEIQHRF")) {
      dat$epi$f.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "f" &
                               dat$attr$group == 1)
    }
    if (dat$control$type == "SIR") {
        dat$epi$num[at] <- dat$epi$s.num[at] +
                             dat$epi$i.num[at] +
                             dat$epi$r.num[at]
    } else if (dat$control$type == "SEIR") {
        dat$epi$num[at] <- dat$epi$s.num[at] +
                             dat$epi$e.num[at] +
                             dat$epi$i.num[at] +
                             dat$epi$r.num[at]
    } else if (dat$control$type == "SEIQHR") {
        dat$epi$num[at] <- dat$epi$s.num[at] +
                             dat$epi$e.num[at] +
                             dat$epi$i.num[at] +
                             dat$epi$q.num[at] +
                             dat$epi$h.num[at] +
                             dat$epi$r.num[at]
    } else if (dat$control$type == "SEIQHRF") {
        dat$epi$num[at] <- dat$epi$s.num[at] +
                             dat$epi$e.num[at] +
                             dat$epi$i.num[at] +
                             dat$epi$q.num[at] +
                             dat$epi$h.num[at] +
                             dat$epi$r.num[at] +
                             dat$epi$f.num[at]
    } else {
        dat$epi$num[at] <- dat$epi$s.num[at] + dat$epi$i.num[at]
    }
        
    if (dat$param$groups == 2) {
      dat$epi$s.num.g2[at] <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "s" &
                                  dat$attr$group == 2)
      if (dat$control$type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi$i.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "e" &
                                    dat$attr$group == 2)
      }
      dat$epi$i.num.g2[at] <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "i" &
                                  dat$attr$group == 2)
      if (dat$control$type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF")) {
        dat$epi$r.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "r" &
                                    dat$attr$group == 2)
      }
      if (dat$control$type %in% c("SEIQHR", "SEIQHRF")) {
        dat$epi$q.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "q" &
                                    dat$attr$group == 2)
        dat$epi$h.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "h" &
                                    dat$attr$group == 2)
      }
      if (dat$control$type %in% c("SEIQHRF")) {
        dat$epi$f.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "f" &
                                    dat$attr$group == 2)
      }
      if (dat$control$type == "SIR") {
        dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
                              dat$epi$i.num.g2[at] +
                              dat$epi$r.num.g2[at]
      } else if (dat$control$type == "SEIR") {
        dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
                              dat$epi$e.num.g2[at] +
                              dat$epi$i.num.g2[at] +
                              dat$epi$r.num.g2[at]
      } else if (dat$control$type == "SEIQHR") {
        dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
                              dat$epi$e.num.g2[at] +
                              dat$epi$i.num.g2[at] +
                              dat$epi$q.num.g2[at] +
                              dat$epi$h.num.g2[at] +
                              dat$epi$r.num.g2[at]
      } else if (dat$control$type == "SEIQHRF") {
        dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
                              dat$epi$e.num.g2[at] +
                              dat$epi$i.num.g2[at] +
                              dat$epi$q.num.g2[at] +
                              dat$epi$h.num.g2[at] +
                              dat$epi$r.num.g2[at] +
                              dat$epi$f.num.g2[at]
      } else {
        dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
                              dat$epi$i.num.g2[at]
      }
    }
  }

  return(dat)
}
