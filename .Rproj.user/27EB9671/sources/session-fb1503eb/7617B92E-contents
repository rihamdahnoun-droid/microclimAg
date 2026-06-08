#' Télécharger des données LST (simulé)
#'
#' @description
#' Simule le téléchargement de données LST (Land Surface Temperature) depuis
#' MODIS ou Landsat. En production, cette fonction devrait utiliser l'API
#' NASA Earthdata ou Google Earth Engine. Ici, des données réalistes sont simulées.
#'
#' @param start_date Date de début "YYYY-MM-DD"
#' @param end_date Date de fin "YYYY-MM-DD"
#' @param bbox Vecteur numérique: c(xmin, xmax, ymin, ymax) en degrés WGS84
#' @param source Source des données: "MODIS" ou "Landsat" (défaut: "MODIS")
#' @param resolution Résolution spatiale en degrés (défaut: 0.01)
#' @param seed Graine aléatoire pour reproductibilité (défaut: 42)
#'
#' @return Un SpatRaster multi-couches (une couche par période)
#'
#' @details
#' Les données MODIS LST (MOD11A1/MYD11A1) ont une résolution native de 1km.
#' Les données Landsat (Band 10) ont une résolution de 30m. Dans cette version
#' de démonstration, des rasters simulés sont retournés.
#'
#' @examples
#' lst <- download_lst_data(
#'   start_date = "2023-06-01",
#'   end_date   = "2023-08-31",
#'   bbox       = c(-5.5, -5.0, 33.5, 34.0),
#'   source     = "MODIS"
#' )
#' terra::plot(lst, main = "LST simulee (°C)")
#'
#' @export
download_lst_data <- function(start_date = "2023-06-01",
                              end_date   = "2023-08-31",
                              bbox       = c(-5.5, -5.0, 33.5, 34.0),
                              source     = "MODIS",
                              resolution = 0.01,
                              seed       = 42) {

  source <- match.arg(source, c("MODIS", "Landsat"))

  message(sprintf("Simulation du telechargement %s LST...", source))
  message(sprintf("Periode: %s a %s", start_date, end_date))
  message(sprintf("Zone: xmin=%.2f xmax=%.2f ymin=%.2f ymax=%.2f", bbox[1], bbox[2], bbox[3], bbox[4]))

  # Calculer nombre de dates mensuelles dans la période
  dates_seq <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  n_dates <- length(dates_seq)
  message(sprintf("Nombre de scenes: %d", n_dates))

  # Dimensions du raster basées sur bbox et résolution
  nrow_r <- max(10, round((bbox[4] - bbox[3]) / resolution))
  ncol_r <- max(10, round((bbox[2] - bbox[1]) / resolution))

  # Simuler rasters LST
  set.seed(seed)
  layers <- list()

  for (i in seq_len(n_dates)) {
    mois <- lubridate::month(dates_seq[i])
    # Saisonnalité thermique : été chaud, hiver frais (Maroc)
    temp_base <- 25 + 12 * sin((mois - 3) * pi / 6)

    x_grad <- matrix(rep(seq(0, 1, length.out = ncol_r), each = nrow_r), nrow = nrow_r)
    y_grad <- matrix(rep(seq(0, 1, length.out = nrow_r), times = ncol_r), nrow = nrow_r)

    vals <- temp_base + 8 * x_grad - 4 * y_grad +
      matrix(stats::rnorm(nrow_r * ncol_r, 0, 1.5), nrow = nrow_r)

    # Ajouter un îlot de chaleur simulé
    cx <- round(ncol_r * 0.65)
    cy <- round(nrow_r * 0.45)
    radius <- round(min(nrow_r, ncol_r) * 0.12)
    for (r in max(1, cy - radius):min(nrow_r, cy + radius)) {
      for (c in max(1, cx - radius):min(ncol_r, cx + radius)) {
        dist <- sqrt((r - cy)^2 + (c - cx)^2)
        if (dist < radius) vals[r, c] <- vals[r, c] + (radius - dist) * 0.8
      }
    }

    r_layer <- terra::rast(vals)
    terra::ext(r_layer) <- terra::ext(bbox[1], bbox[2], bbox[3], bbox[4])
    terra::crs(r_layer) <- "EPSG:4326"
    names(r_layer) <- format(dates_seq[i], "%Y-%m")
    layers[[i]] <- r_layer
  }

  lst_raster <- terra::rast(layers)
  message(sprintf("LST telecharge: %d couches, resolution ~%.4f deg", terra::nlyr(lst_raster), resolution))
  lst_raster
}


#' Nettoyer les données LST
#'
#' @description
#' Nettoie un raster LST en supprimant les pixels nuageux, corrigeant les
#' valeurs aberrantes et interpolant les valeurs manquantes.
#'
#' @param lst_raster Un objet SpatRaster contenant les données LST
#' @param min_temp Température minimale valide en °C (défaut: -10)
#' @param max_temp Température maximale valide en °C (défaut: 65)
#' @param cloud_threshold Seuil de détection nuages (défaut: NULL, désactivé)
#' @param interpolate_na Interpoler les valeurs manquantes (défaut: TRUE)
#' @param smooth Lisser le raster (fenêtre de lissage en cellules, défaut: NULL)
#'
#' @return Un SpatRaster nettoyé
#'
#' @examples
#' lst <- simulate_raster_data(type = "lst", n_dates = 3)
#' lst_clean <- clean_lst_data(lst, min_temp = 5, max_temp = 55)
#' terra::plot(lst_clean[[1]], main = "LST nettoyee")
#'
#' @export
clean_lst_data <- function(lst_raster,
                           min_temp = -10,
                           max_temp = 65,
                           cloud_threshold = NULL,
                           interpolate_na = TRUE,
                           smooth = NULL) {

  if (!inherits(lst_raster, "SpatRaster")) {
    stop("lst_raster doit etre un objet SpatRaster (package terra)")
  }

  message("Nettoyage des donnees LST...")
  n_total <- terra::ncell(lst_raster) * terra::nlyr(lst_raster)

  # 1. Masquer valeurs hors plage physique
  lst_clean <- terra::app(lst_raster, function(x) {
    x[x < min_temp | x > max_temp] <- NA
    x
  })

  n_masked <- sum(is.na(terra::values(lst_clean))) - sum(is.na(terra::values(lst_raster)))
  message(sprintf("  Pixels aberrants masques: %d (%.1f%%)",
                  max(0, n_masked), max(0, n_masked) / n_total * 100))

  # 2. Interpolation spatiale des NA par fenêtre glissante
  if (interpolate_na) {
    n_na_before <- sum(is.na(terra::values(lst_clean)))
    if (n_na_before > 0) {
      lst_clean <- terra::focal(lst_clean, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
      n_na_after <- sum(is.na(terra::values(lst_clean)))
      message(sprintf("  NA interpoles: %d -> %d restants", n_na_before, n_na_after))
    }
  }

  # 3. Lissage spatial optionnel
  if (!is.null(smooth)) {
    lst_clean <- terra::focal(lst_clean, w = smooth, fun = "mean", na.rm = TRUE)
    message(sprintf("  Lissage applique (fenetre: %dx%d)", smooth, smooth))
  }

  message("Nettoyage termine.")
  lst_clean
}


#' Calculer les anomalies thermiques LST
#'
#' @description
#' Calcule les anomalies thermiques par rapport à la moyenne saisonnière ou
#' via un z-score. Permet d'identifier les zones anormalement chaudes ou froides.
#'
#' @param lst_raster Un SpatRaster multi-couches de LST
#' @param method Méthode de calcul: "difference" (écart à la moyenne) ou
#'   "zscore" (z-score normalisé) (défaut: "difference")
#' @param reference Couche de référence: "mean" (moyenne de toutes les couches)
#'   ou un index entier de couche (défaut: "mean")
#'
#' @return Un SpatRaster d'anomalies thermiques
#'
#' @examples
#' lst <- simulate_raster_data(type = "lst", n_dates = 6)
#' anom <- calculate_lst_anomaly(lst, method = "zscore")
#' terra::plot(anom, main = "Anomalies LST (Z-score)")
#'
#' @export
calculate_lst_anomaly <- function(lst_raster,
                                  method    = "difference",
                                  reference = "mean") {

  method <- match.arg(method, c("difference", "zscore"))

  if (terra::nlyr(lst_raster) < 2) {
    stop("Il faut au moins 2 couches pour calculer des anomalies.")
  }

  message(sprintf("Calcul anomalies LST (methode: %s)...", method))

  if (method == "difference") {
    if (reference == "mean") {
      ref_raster <- terra::app(lst_raster, fun = "mean", na.rm = TRUE)
    } else {
      ref_raster <- lst_raster[[as.integer(reference)]]
    }
    anomalies <- lst_raster - ref_raster

  } else if (method == "zscore") {
    lst_mean <- terra::app(lst_raster, fun = "mean", na.rm = TRUE)
    lst_sd   <- terra::app(lst_raster, fun = "sd",   na.rm = TRUE)
    # Éviter division par zéro
    lst_sd[lst_sd < 0.001] <- 0.001
    anomalies <- (lst_raster - lst_mean) / lst_sd
  }

  names(anomalies) <- paste0("anomalie_", names(lst_raster))

  # Statistiques
  vals <- terra::values(anomalies[[1]])
  message(sprintf("  Anomalie min: %.2f, max: %.2f, ecart-type: %.2f",
                  min(vals, na.rm = TRUE),
                  max(vals, na.rm = TRUE),
                  stats::sd(vals, na.rm = TRUE)))

  anomalies
}
