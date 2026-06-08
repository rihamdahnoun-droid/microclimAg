#' microclimAg: Analyse du Microclimat Agricole et des Ilots de Chaleur Ruraux
#'
#' @description
#' Le package \code{microclimAg} fournit des outils pour analyser les variations
#' microclimatiques dans les paysages agricoles. Il integre des donnees de
#' teledetection thermique (LST), des indices de vegetation (NDVI), des donnees
#' meteorologiques et des informations sur l'occupation du sol.
#'
#' @section Fonctions principales:
#' \describe{
#'   \item{\code{\link{import_weather_data}}}{Import des donnees meteorologiques}
#'   \item{\code{\link{import_station_locations}}}{Import des stations meteo}
#'   \item{\code{\link{calculate_ndvi}}}{Calcul du NDVI}
#'   \item{\code{\link{calculate_lst_anomaly}}}{Calcul des anomalies LST}
#'   \item{\code{\link{cluster_thermal_zones}}}{Classification des zones thermiques}
#'   \item{\code{\link{plot_microclimate_map}}}{Cartographie microclimatique}
#' }
#'
#' @docType package
#' @name microclimAg-package
#' @aliases microclimAg
"_PACKAGE"

#' Simuler des donnees meteorologiques
#'
#' @description
#' Genere un jeu de donnees meteorologiques simule pour les tests et demonstrations.
#' Utilise quand aucune donnee reelle n'est disponible.
#'
#' @param n_days Nombre de jours a simuler (defaut: 365)
#' @param n_stations Nombre de stations meteo (defaut: 5)
#' @param start_date Date de debut au format "YYYY-MM-DD" (defaut: "2023-01-01")
#' @param seed Graine aleatoire pour reproductibilite (defaut: 42)
#'
#' @return Un data.frame avec les colonnes: date, station_id, temperature,
#'   humidity, wind_speed, radiation
#'
#' @examples
#' weather_sim <- simulate_weather_data(n_days = 30, n_stations = 3)
#' head(weather_sim)
#'
#' @export
simulate_weather_data <- function(n_days = 365, n_stations = 5,
                                  start_date = "2023-01-01", seed = 42) {
  set.seed(seed)

  dates       <- seq(as.Date(start_date), by = "day", length.out = n_days)
  station_ids <- paste0("ST", sprintf("%02d", seq_len(n_stations)))

  grid        <- expand.grid(date = dates, station_id = station_ids,
                             stringsAsFactors = FALSE)
  grid$date   <- as.Date(grid$date)
  n_rows      <- nrow(grid)

  day_of_year     <- as.numeric(format(grid$date, "%j"))
  temp_seasonal   <- 20 + 12 * sin((day_of_year - 80) * 2 * pi / 365)

  grid$temperature <- round(temp_seasonal + stats::rnorm(n_rows, 0, 2.5), 1)
  grid$humidity    <- round(pmax(10, pmin(95,
                                          60 - 0.8 * temp_seasonal + stats::rnorm(n_rows, 0, 10))), 1)
  grid$wind_speed  <- round(pmax(0, stats::rnorm(n_rows, 3, 1.5)), 2)
  grid$radiation   <- round(pmax(0,
                                 250 + 150 * sin((day_of_year - 80) * 2 * pi / 365) +
                                   stats::rnorm(n_rows, 0, 30)), 1)

  station_offsets        <- stats::rnorm(n_stations, 0, 1.5)
  names(station_offsets) <- station_ids
  grid$temperature       <- round(grid$temperature +
                                    station_offsets[grid$station_id], 1)

  grid[order(grid$date, grid$station_id), ]
}

#' Simuler des donnees raster (LST et NDVI)
#'
#' @description
#' Genere des rasters simules de LST et NDVI pour les demonstrations.
#'
#' @param nrow Nombre de lignes du raster (defaut: 50)
#' @param ncol Nombre de colonnes du raster (defaut: 50)
#' @param n_dates Nombre de dates (couches temporelles) (defaut: 3)
#' @param type Type de raster: "lst" ou "ndvi" (defaut: "lst")
#' @param seed Graine aleatoire (defaut: 42)
#'
#' @return Un objet SpatRaster (package terra)
#'
#' @examples
#' lst_sim  <- simulate_raster_data(type = "lst",  n_dates = 2)
#' ndvi_sim <- simulate_raster_data(type = "ndvi", n_dates = 1)
#'
#' @export
simulate_raster_data <- function(nrow = 50, ncol = 50, n_dates = 3,
                                 type = "lst", seed = 42) {
  set.seed(seed)

  type <- match.arg(type, c("lst", "ndvi"))

  layers <- list()
  for (i in seq_len(n_dates)) {

    x_grad <- matrix(rep(seq_len(ncol), each = nrow),  nrow = nrow)
    y_grad <- matrix(rep(seq_len(nrow), times = ncol), nrow = nrow)

    if (type == "lst") {
      vals <- 30 + 10 * (x_grad / ncol) - 5 * (y_grad / nrow) +
        matrix(stats::rnorm(nrow * ncol, 0, 2), nrow = nrow)
      # Zone chaude simulee
      cx <- round(ncol * 0.6); cy <- round(nrow * 0.4)
      for (r in max(1L, cy - 8L):min(nrow, cy + 8L)) {
        for (c in max(1L, cx - 8L):min(ncol, cx + 8L)) {
          d <- sqrt((r - cy)^2 + (c - cx)^2)
          vals[r, c] <- vals[r, c] + pmax(0, 8 - d)
        }
      }
    } else {
      vals <- 0.3 - 0.3 * (x_grad / ncol) + 0.2 * (y_grad / nrow) +
        matrix(stats::rnorm(nrow * ncol, 0, 0.05), nrow = nrow)
      vals <- pmax(-0.2, pmin(0.9, vals))
    }

    # Creer le raster via nrow/ncol/nlyrs puis assigner les valeurs
    r_obj <- terra::rast(nrows = nrow, ncols = ncol, nlyrs = 1,
                         xmin = -5.5, xmax = -5.0,
                         ymin = 33.5, ymax = 34.0,
                         crs  = "EPSG:4326")
    terra::values(r_obj) <- as.vector(t(vals))
    names(r_obj) <- paste0(type, "_date", i)
    layers[[i]] <- r_obj
  }

  if (n_dates == 1L) return(layers[[1L]])
  terra::rast(layers)
}
