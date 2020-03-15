icm.seiqhrf <- function(param, init, control) {

  crosscheck.icm(param, init, control)
  verbose.icm(control, type = "startup")

  # Simulation loop start
  for (s in 1:control$nsims) {

    ## Initialization module
    if (!is.null(control[["initialize.FUN"]])) {
      dat <- do.call(control[["initialize.FUN"]], list(param, init, control))
    }


    # Timestep loop
    for (at in 2:control$nsteps) {

      ## User Modules
      um <- control$user.mods
      if (length(um) > 0) {
        for (i in 1:length(um)) {
          dat <- do.call(control[[um[i]]], list(dat, at))
        }
      }

      ## Infection
      if (!is.null(control[["infection.FUN"]])) {
        dat <- do.call(control[["infection.FUN"]], list(dat, at))
      }


      ## Recovery
      if (!is.null(control[["recovery.FUN"]])) {
        dat <- do.call(control[["recovery.FUN"]], list(dat, at))
      }


      ## Departure Module
      if (!is.null(control[["departures.FUN"]])) {
        dat <- do.call(control[["departures.FUN"]], list(dat, at))
      }


      ## Arrival Module
      if (!is.null(control[["arrivals.FUN"]])) {
        dat <- do.call(control[["arrivals.FUN"]], list(dat, at))
      }


      ## Outputs
      if (!is.null(control[["get_prev.FUN"]])) {
        dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
      }


      ## Track progress
      verbose.icm(dat, type = "progress", s, at)
    }

    # Set output
    if (s == 1) {
      out <- saveout.seiqhrf.icm(dat, s)
    } else {
      out <- saveout.seiqhrf.icm(dat, s, out)
    }



  } # Simulation loop end


  class(out) <- "icm"
  return(out)
}
