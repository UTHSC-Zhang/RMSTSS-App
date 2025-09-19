
# R/generate_sets.R
# Generate datasets to disk (CSV/TXT/RDS/RData) + manifest. No L/tau recorded.

`%||%` <- function(x, y) if (is.null(x)) y else x

.collect_metadata <- function(dat, recipe, files = list(),
                              seed, scenario_id, rep, sc_params = list()) {
  has_arm <- "arm" %in% names(dat)
  list(
    dataset_id = sprintf("sc%03d_r%02d", scenario_id, rep),
    scenario_id = scenario_id,
    rep = rep,
    seed_used = seed %||% NA_integer_,
    n = nrow(dat),
    n_treat = if (has_arm) sum(dat$arm == 1L) else NA_integer_,  # may be NA if no arm
    n_control = if (has_arm) sum(dat$arm == 0L) else NA_integer_,
    event_rate = mean(dat$status == 1L),
    achieved_censoring = attr(dat, "achieved_censoring", exact = TRUE) %||% mean(dat$status == 0L),
    model = recipe$event_time$model,
    baseline = recipe$event_time$baseline,
    effects  = recipe$event_time$effects,
    treatment = recipe$treatment %||% NULL,
    censoring = recipe$censoring %||% NULL,
    covariates = if (!is.null(recipe$covariates$defs))
      lapply(recipe$covariates$defs, function(d) d[c("name","type","dist","params")]) else NULL,
    allocation = recipe$treatment$allocation %||% NA_character_,
    params = sc_params,
    files = files,
    created_at = as.character(Sys.time())
  )
}

#' Generate simulated datasets across scenario combinations (TXT/CSV/RDS/RData)
#'
#' Builds a grid of scenarios from a base recipe and a named list of variations,
#' simulates one or more replicates per scenario, and writes the datasets to files
#' in a target folder as \code{.txt} (tab), \code{.csv}, \code{.rds}, and/or \code{.RData}.
#' A manifest is written as \code{manifest.rds}. No YAML and **no L/tau**.
#'
#' @param base_recipe A recipe list (use \code{validate_recipe()} if needed).
#' @param vary Named list; keys are dotted paths inside the recipe (e.g., \code{"n"},
#'   \code{"censoring.target"}, \code{"event_time.effects.treatment"}).
#' @param out_dir Directory to write datasets and \code{manifest.rds} (created if missing).
#' @param formats Character vector subset of \code{c("txt","csv","rds","rdata")}.
#' @param n_reps Integer; number of replicates per scenario.
#' @param seed_base Optional integer; per-rep seed computed as
#'   \code{seed_base + scenario_id*1000 + rep}.
#' @param filename_template Base filename (no extension) with placeholders:
#'   \code{"{scenario_id}"}, \code{"{rep}"} and any dotted path used in \code{vary}.
#' @return Invisibly returns the manifest \code{data.frame} and writes \code{manifest.rds}.
#' @examples
#' covs <- list(list(name="x", type="continuous", dist="normal", params=list(mean=0, sd=1)))
#' rec <- recipe_quick_aft(60, "aft_lognormal",
#'          baseline=list(mu=2.7, sigma=0.6), treat_effect=-0.2,
#'          covariates=covs, target_censoring=0.2)
#' dir.create("checks", showWarnings = FALSE)
#' man <- generate_recipe_sets(rec, vary=list(n=c(60,80)), out_dir="checks",
#'          formats=c("csv","rds"), n_reps=1, seed_base=123)
#' @export
generate_recipe_sets <- function(base_recipe,
                                 vary = list(),
                                 out_dir,
                                 formats = c("txt","csv","rds","rdata"),
                                 n_reps = 1L,
                                 seed_base = NULL,
                                 filename_template = "sc{scenario_id}_r{rep}") {
  if (missing(out_dir) || !nzchar(out_dir)) stop("out_dir is required")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  formats <- unique(tolower(formats))
  allowed <- c("txt","csv","rds","rdata")
  if (!length(formats)) formats <- allowed
  if (!all(formats %in% allowed)) stop("formats must be subset of: ", paste(allowed, collapse=", "))
  if (!is.list(vary)) stop("vary must be a named list (possibly empty)")

  base_recipe <- validate_recipe(base_recipe)
  grids <- recipe_grid(base_recipe, vary)
  if (!length(grids)) grids <- list(base_recipe)

  `%||%` <- function(x, y) if (is.null(x)) y else x
  sanitize <- function(x) { x <- gsub("[^A-Za-z0-9._-]+", "_", x); x <- gsub("_+", "_", x); x <- sub("^_", "", x); x <- sub("_$", "", x); x }
  get_by_path <- function(obj, path) { parts <- strsplit(path, ".", fixed = TRUE)[[1]]; cur <- obj; for (p in parts) cur <- cur[[p]]; cur }
  decl_levels <- function(recipe) {
    out <- list(); defs <- recipe$covariates$defs
    for (d in defs) if (identical(d$type, "categorical") && d$dist == "categorical" && !is.null(d$params$labels))
      out[[d$name]] <- d$params$labels
    out
  }
  build_name <- function(template, tokens) {
    out <- template
    for (nm in names(tokens)) {
      ph <- paste0("{", nm, "}")
      val <- tokens[[nm]]
      if (is.numeric(val) && length(val) == 1) val <- format(val, digits = 6, trim = TRUE, scientific = FALSE)
      if (length(val) > 1) val <- paste(val, collapse = "-")
      out <- gsub(ph, as.character(val), out, fixed = TRUE)
    }
    sanitize(out)
  }

  manifest_rows <- list()
  sc_id <- 0L
  for (sc in seq_along(grids)) {
    r <- grids[[sc]]; sc_id <- sc_id + 1L
    sc_params <- list(); if (length(vary)) for (k in names(vary)) sc_params[[k]] <- get_by_path(r, k)
    levels_decl <- decl_levels(r)
    for (rep in seq_len(n_reps)) {
      seed_use <- if (!is.null(seed_base)) as.integer(seed_base + sc_id * 1000 + rep) else NULL
      dat <- simulate_from_recipe(r, seed = seed_use)

      if (length(levels_decl)) for (nm in names(levels_decl)) if (!is.null(dat[[nm]]))
        dat[[nm]] <- factor(dat[[nm]], levels = levels_decl[[nm]])

      base_name <- build_name(filename_template, c(list(scenario_id = sc_id, rep = rep), sc_params))

      paths <- list(txt = NA_character_, csv = NA_character_, rds = NA_character_, rdata = NA_character_)
      if ("txt"   %in% formats) { p <- file.path(out_dir, paste0(base_name, ".txt"));   utils::write.table(dat, p, sep="\t", row.names=FALSE, quote=FALSE); paths$txt   <- p }
      if ("csv"   %in% formats) { p <- file.path(out_dir, paste0(base_name, ".csv"));   utils::write.csv(dat, p, row.names=FALSE);                              paths$csv   <- p }
      if ("rds"   %in% formats) { p <- file.path(out_dir, paste0(base_name, ".rds"));   saveRDS(dat, p);                                                       paths$rds   <- p }
      if ("rdata" %in% formats) { p <- file.path(out_dir, paste0(base_name, ".RData")); dat_obj <- dat; meta <- list(seed=seed_use, achieved_censoring=attr(dat,"achieved_censoring"), scenario_id=sc_id, rep=rep, params=sc_params, n=nrow(dat)); save(dat_obj, meta, file=p); paths$rdata <- p }

      meta_full <- .collect_metadata(
         dat        = dat,
         recipe     = r,
         files      = paths,
         seed       = seed_use,
         scenario_id = sc_id,
         rep        = rep,
         sc_params  = sc_params
      )

      row <- data.frame(
         scenario_id = sc_id, rep = rep, seed = seed_use %||% NA_integer_,
         achieved_censoring = meta_full$achieved_censoring, n = nrow(dat),
         file_txt = paths$txt, file_csv = paths$csv, file_rds = paths$rds, file_rdata = paths$rdata,
         stringsAsFactors = FALSE
      )
      if (length(sc_params)) {
         add <- as.list(sc_params); names(add) <- paste0("p__", sanitize(names(add)))
         for (nm in names(add)) row[[nm]] <- as.character(add[[nm]])
      }
      row$meta <- I(list(meta_full))
      manifest_rows[[length(manifest_rows) + 1L]] <- row
    }
  }
  manifest <- do.call(rbind, manifest_rows)
  saveRDS(manifest, file.path(out_dir, "manifest.rds"))
  invisible(manifest)
}

#' Expand a recipe over a grid of values (list-only)
#' @param base A recipe list.
#' @param vary Named list of vectors; keys are dotted paths.
#' @return A list of recipe lists.
#' @examples
#' r <- recipe_quick_aft(60, "aft_lognormal", baseline=list(mu=2.7, sigma=0.6),
#'        treat_effect=-0.2, covariates=list(list(name="x", type="continuous", dist="normal",
#'        params=list(mean=0, sd=1))))
#' recipe_grid(r, list(n=c(60,80), "event_time.effects.treatment"=c(-0.2,-0.4)))
#' @export
recipe_grid <- function(base, vary) {
  if (!length(vary)) return(list(base))
  keys <- names(vary); grid <- expand.grid(vary, stringsAsFactors = FALSE)
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    r <- base
    for (k in keys) {
      parts <- strsplit(k, ".", fixed = TRUE)[[1]]
      r <- .set_by_path(r, parts, grid[[i, k]])
    }
    out[[i]] <- r
  }
  out
}

.set_by_path <- function(obj, parts, value) {
  if (length(parts) == 1) { obj[[parts[1]]] <- value; return(obj) }
  p <- parts[1]
  obj[[p]] <- .set_by_path(obj[[p]], parts[-1], value)
  obj
}
