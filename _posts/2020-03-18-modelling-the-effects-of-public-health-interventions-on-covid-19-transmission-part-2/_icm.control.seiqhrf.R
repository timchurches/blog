control.icm <- function(type, nsteps, nsims = 1, 
                        rec.rand = TRUE, quar.rand = TRUE, hosp.rand = TRUE, disch.rand = TRUE,
                        fat.rand = TRUE, a.rand = TRUE, d.rand = TRUE, initialize.FUN = initialize.icm,
                        infection.FUN = infection.icm, recovery.FUN = recovery.icm,
                        departures.FUN = departures.icm, arrivals.FUN = arrivals.icm,
                        get_prev.FUN = get_prev.icm, verbose = FALSE,
                        verbose.int = 0, skip.check = FALSE, ncores=1, ...) {

  # Get arguments
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    if (as.logical(mget(arg) != "")) {
      p[arg] <- list(get(arg))
    }
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  if ("births.FUN" %in% names(dot.args)) {
    p$arrivals.FUN <- dot.args$births.FUN
    p$births.FUN <- dot.args$births.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the birth function births.FUN to arrivals.FUN. See documentation for details.")
  }
  if ("deaths.FUN" %in% names(dot.args)) {
    p$departures.FUN <- dot.args$deaths.FUN
    p$deaths.FUN <- dot.args$deaths.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the death function deaths.FUN to departures.FUN. See documentation for details.")
  }


  ## Module classification
  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)


  ## Defaults and checks
  if (is.null(p$type) | !(p$type %in% c("SI", "SIS", "SIR", "SEIR", "SEIQHR", "SEIQHRF"))) {
    stop("Specify type as \"SI\", \"SIS\", \"SIR\", \"SEIR\", \"SEIQHR\", or \"SEIQHRF\" ", call. = FALSE)
  }
  if (is.null(p$nsteps)) {
    stop("Specify nsteps", call. = FALSE)
  }


  ## Output
  class(p) <- c("control.icm", "list")
  return(p)
}
