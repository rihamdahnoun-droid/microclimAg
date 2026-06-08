#' Interpoler la température spatialement
#'
#' @description
#' Interpole les températures des stations météo pour créer une surface
#' continue. Méthodes disponibles: IDW (Inverse Distance Weighting) ou
#' krigeage simplifié.
#'
#' @param stations_sf Objet sf avec les stations et températures
#' @param temp_col Nom de la colonne température dans stations_sf (défaut: "temperature")
#' @param method Méthode: "idw" ou "kriging" (défaut: "idw")
#' @param resolution Résolution de la grille de sortie en degrés (défaut: 0.01)
#' @param idw_power Puissance IDW (défaut: 2)
#'
#' @return Un SpatRaster de température interpolée
#'
#' @examples
#' # Créer des données de stations simulées
#' set.seed(42)
#' stations_df <- data.frame(
#'   station_id  = paste0("ST", 1:10),
#'   longitude   = runif(10, -5.5, -5.0),
#'   latitude    = runif(10, 33.5, 34.0),
#'   temperature = runif(10, 28, 40)
#' )
#' stations_sf <- sf::st_as_sf(stations_df, coords = c("longitude", "latitude"), crs = 4326)
#' temp_raster <- interpolate_temperature(stations_sf, method = "idw")
#' terra::plot(temp_raster, main = "Temperature interpolee IDW (°C)")
#'
#' @export
interpolate_temperature <- function(stations_sf,
                                    temp_col  = "temperature",
                                    method    = "idw",
                                    resolution = 0.01,
                                    idw_power = 2) {

  method <- match.arg(method, c("idw", "kriging"))

  if (!temp_col %in% names(stations_sf)) {
    stop(sprintf("Colonne '%s' introuvable dans stations_sf", temp_col))
  }

  coords <- sf::st_coordinates(stations_sf)
  temps  <- stations_sf[[temp_col]]
  valid  <- !is.na(temps)
  coords <- coords[valid, ]
  temps  <- temps[valid]

  # Créer grille cible
  bbox <- sf::st_bbox(stations_sf)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], by = resolution)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], by = resolution)
  grid_xy <- expand.grid(x = x_seq, y = y_seq)

  message(sprintf("Interpolation %s: %d stations -> grille %dx%d",
                  toupper(method), sum(valid), length(x_seq), length(y_seq)))

  if (method == "idw") {
    # IDW manuel
    z_interp <- apply(grid_xy, 1, function(pt) {
      d <- sqrt((coords[,1] - pt[1])^2 + (coords[,2] - pt[2])^2)
      if (any(d == 0)) return(temps[which(d == 0)[1]])
      w <- 1 / d^idw_power
      sum(w * temps) / sum(w)
    })

  } else if (method == "kriging") {
    # Krigeage simplifié (variogramme exponentiel)
    range_kr <- mean(c(diff(range(coords[,1])), diff(range(coords[,2])))) * 0.5
    sill_kr  <- stats::var(temps)

    z_interp <- apply(grid_xy, 1, function(pt) {
      d <- sqrt((coords[,1] - pt[1])^2 + (coords[,2] - pt[2])^2)
      gamma_vals <- sill_kr * (1 - exp(-d / range_kr))
      w <- pmax(0, sill_kr - gamma_vals)
      if (sum(w) == 0) return(mean(temps))
      sum(w * temps) / sum(w)
    })
  }

  # Créer raster
  mat <- matrix(z_interp, nrow = length(y_seq), ncol = length(x_seq), byrow = TRUE)
  r_out <- terra::rast(mat)
  terra::ext(r_out) <- terra::ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
  terra::crs(r_out) <- "EPSG:4326"
  names(r_out) <- paste0("temp_", method)

  message(sprintf("  Plage: %.1f - %.1f °C", min(z_interp), max(z_interp)))
  r_out
}


#' Analyser la corrélation spatiale (Moran's I)
#'
#' @description
#' Calcule l'indice de Moran's I pour détecter l'autocorrélation spatiale
#' dans les données thermiques. Un I positif indique des clusters thermiques,
#' un I négatif indique une dispersion.
#'
#' @param features_df Data.frame avec colonnes x, y, et une variable cible
#' @param variable Nom de la variable à analyser (défaut: "lst")
#' @param k_neighbours Nombre de voisins pour la matrice de poids (défaut: 5)
#' @param n_sample Nombre de points max à utiliser (défaut: 300)
#'
#' @return Une liste avec: moran_i, p_value, interpretation, et le dataframe
#'
#' @examples
#' lst  <- simulate_raster_data(type = "lst", n_dates = 1)
#' features <- extract_microclimate_features(lst, n_sample = 200)
#' result <- analyze_spatial_correlation(features, variable = "lst")
#' print(result$interpretation)
#'
#' @export
analyze_spatial_correlation <- function(features_df,
                                        variable      = "lst",
                                        k_neighbours  = 5,
                                        n_sample      = 300) {

  if (!variable %in% names(features_df)) {
    stop(sprintf("Variable '%s' introuvable dans features_df", variable))
  }

  df <- features_df[!is.na(features_df[[variable]]), ]

  # Limiter le nombre de points pour la performance
  if (nrow(df) > n_sample) {
    df <- df[sample(nrow(df), n_sample), ]
  }

  message(sprintf("Moran's I sur '%s' (%d points, %d voisins)...", variable, nrow(df), k_neighbours))

  coords <- as.matrix(df[, c("x", "y")])
  z <- df[[variable]]

  # Matrice de poids k-voisins les plus proches (calcul direct)
  n <- nrow(coords)
  W <- matrix(0, n, n)

  for (i in seq_len(n)) {
    d <- sqrt(rowSums(sweep(coords, 2, coords[i,], "-")^2))
    d[i] <- Inf
    knn_idx <- order(d)[seq_len(min(k_neighbours, n-1))]
    W[i, knn_idx] <- 1
  }

  # Standardiser W par rangée
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1
  W_std <- W / row_sums

  # Calcul Moran's I
  z_centered <- z - mean(z)
  S0 <- sum(W_std)
  I <- (n / S0) * (sum(W_std * outer(z_centered, z_centered)) / sum(z_centered^2))

  # Valeur attendue et variance sous H0
  E_I <- -1 / (n - 1)
  Var_I <- (n^2 * (n - 3) * n) / ((n - 1)^2 * (n + 1) * (n + 3) * S0^2)
  Var_I <- max(Var_I, 1e-10)

  z_score <- (I - E_I) / sqrt(Var_I)
  p_value <- 2 * stats::pnorm(abs(z_score), lower.tail = FALSE)

  # Interprétation
  interp <- if (I > 0.3 && p_value < 0.05) {
    "Forte autocorrelation positive : clusters thermiques significatifs (ilots de chaleur)"
  } else if (I > 0.1 && p_value < 0.05) {
    "Autocorrelation positive moderee : tendance au regroupement thermique"
  } else if (I < -0.1 && p_value < 0.05) {
    "Autocorrelation negative : distribution dispersee des temperatures"
  } else {
    "Pas d'autocorrelation significative : distribution spatiale aleatoire"
  }

  message(sprintf("  Moran's I = %.4f, p = %.4f", I, p_value))
  message(sprintf("  %s", interp))

  list(
    moran_i       = round(I, 4),
    expected_i    = round(E_I, 4),
    z_score       = round(z_score, 4),
    p_value       = round(p_value, 4),
    significant   = p_value < 0.05,
    interpretation = interp,
    data          = df
  )
}


#' Classifier les zones thermiques par clustering
#'
#' @description
#' Applique un clustering K-means sur les données LST pour identifier des
#' zones thermiques (fraîches, modérées, chaudes). Retourne un raster classifié.
#'
#' @param lst_raster SpatRaster de LST
#' @param n_clusters Nombre de clusters (défaut: 3)
#' @param layer_index Index de la couche LST (défaut: 1)
#' @param seed Graine aléatoire (défaut: 42)
#' @param max_iter Iterations max K-means (défaut: 100)
#'
#' @return Une liste: $raster (SpatRaster classifié), $stats (statistiques par classe),
#'   $labels (étiquettes des classes)
#'
#' @examples
#' lst    <- simulate_raster_data(type = "lst", n_dates = 1)
#' result <- cluster_thermal_zones(lst, n_clusters = 3)
#' terra::plot(result$raster, main = "Zones thermiques K-means",
#'             col = c("steelblue", "gold", "firebrick"))
#' print(result$stats)
#'
#' @export
cluster_thermal_zones <- function(lst_raster,
                                  n_clusters  = 3,
                                  layer_index = 1,
                                  seed        = 42,
                                  max_iter    = 100) {

  message(sprintf("Clustering thermique K-means (k=%d)...", n_clusters))

  lst_layer <- lst_raster[[min(layer_index, terra::nlyr(lst_raster))]]
  vals <- terra::values(lst_layer)
  valid_idx <- !is.na(vals)

  if (sum(valid_idx) < n_clusters * 10) {
    stop("Pas assez de pixels valides pour le clustering.")
  }

  set.seed(seed)
  km <- stats::kmeans(vals[valid_idx], centers = n_clusters,
                      nstart = 10, iter.max = max_iter)

  # Réordonner les clusters par température croissante
  cluster_means <- tapply(vals[valid_idx], km$cluster, mean)
  rank_order    <- rank(cluster_means)
  km$cluster    <- rank_order[km$cluster]

  # Reconstruire le raster
  cluster_vals <- rep(NA_real_, terra::ncell(lst_layer))
  cluster_vals[valid_idx] <- km$cluster
  cluster_raster <- terra::rast(matrix(cluster_vals,
                                       nrow = terra::nrow(lst_layer),
                                       ncol = terra::ncol(lst_layer)))
  terra::ext(cluster_raster) <- terra::ext(lst_layer)
  terra::crs(cluster_raster) <- terra::crs(lst_layer)
  names(cluster_raster) <- "thermal_cluster"

  # Étiquettes
  if (n_clusters == 3) {
    labels <- c("Zone fraiche", "Zone moderee", "Zone chaude")
  } else if (n_clusters == 4) {
    labels <- c("Tres fraiche", "Fraiche", "Chaude", "Tres chaude")
  } else {
    labels <- paste("Cluster", 1:n_clusters)
  }

  # Statistiques par cluster
  sorted_centers <- sort(tapply(vals[valid_idx], km$cluster, mean))
  stats_df <- data.frame(
    cluster = 1:n_clusters,
    label   = labels[1:n_clusters],
    temp_mean = round(as.numeric(sorted_centers), 2),
    temp_sd   = round(tapply(vals[valid_idx], km$cluster, stats::sd)[order(cluster_means)], 2),
    n_pixels  = as.integer(table(km$cluster)[order(cluster_means)])
  )

  message("  Statistiques par cluster:")
  for (i in 1:nrow(stats_df)) {
    message(sprintf("  [%s] moy=%.1f°C, sd=%.2f, n=%d",
                    stats_df$label[i], stats_df$temp_mean[i],
                    stats_df$temp_sd[i], stats_df$n_pixels[i]))
  }

  list(raster = cluster_raster, stats = stats_df, labels = labels)
}


#' Modéliser les relations microclimatiques
#'
#' @description
#' Modélise la relation entre la température de surface (LST) et des
#' variables explicatives (NDVI, altitude, occupation du sol) par régression
#' linéaire.
#'
#' @param features_df Data.frame avec les variables microclimatiques
#' @param target Nom de la variable cible (défaut: "lst")
#' @param predictors Vecteur de noms de prédicteurs (défaut: NULL = auto)
#' @param method Méthode: "lm" (régression linéaire) (défaut: "lm")
#'
#' @return Une liste avec: model, summary, coefficients, r_squared, importance
#'
#' @examples
#' lst     <- simulate_raster_data(type = "lst",  n_dates = 1)
#' ndvi    <- simulate_raster_data(type = "ndvi", n_dates = 1)
#' lc      <- import_landcover()
#' features <- extract_microclimate_features(lst, ndvi, lc)
#' model    <- model_temperature_relationships(features)
#' print(model$coefficients)
#' cat("R² =", model$r_squared, "\n")
#'
#' @export
model_temperature_relationships <- function(features_df,
                                            target     = "lst",
                                            predictors = NULL,
                                            method     = "lm") {

  if (!target %in% names(features_df)) {
    stop(sprintf("Variable cible '%s' introuvable.", target))
  }

  # Sélectionner prédicteurs numériques disponibles
  if (is.null(predictors)) {
    candidate <- setdiff(names(features_df), c(target, "x", "y", "land_cover_name"))
    predictors <- candidate[sapply(features_df[candidate], is.numeric)]
    predictors <- predictors[sapply(features_df[predictors],
                                    function(x) sum(!is.na(x)) > 10)]
  }

  if (length(predictors) == 0) {
    stop("Aucun predicteur numerique disponible.")
  }

  message(sprintf("Regression '%s' ~ %s", target, paste(predictors, collapse=" + ")))

  # Formuler et ajuster le modèle
  formula_str <- paste(target, "~", paste(predictors, collapse = " + "))
  df_clean <- features_df[, c(target, predictors)]
  df_clean <- df_clean[complete.cases(df_clean), ]

  if (nrow(df_clean) < length(predictors) + 5) {
    stop("Pas assez de donnees completes pour la regression.")
  }

  model <- stats::lm(stats::as.formula(formula_str), data = df_clean)
  summ  <- summary(model)

  # Importance des variables (valeur absolue des coefficients standardisés)
  scaled_df <- as.data.frame(scale(df_clean))
  model_scaled <- stats::lm(stats::as.formula(formula_str), data = scaled_df)
  importance <- abs(stats::coef(model_scaled)[-1])
  importance <- sort(importance, decreasing = TRUE)

  message(sprintf("  R² = %.4f, R² ajuste = %.4f", summ$r.squared, summ$adj.r.squared))
  message("  Variable la plus importante: ", names(importance)[1])

  list(
    model        = model,
    summary      = summ,
    coefficients = as.data.frame(summ$coefficients),
    r_squared    = round(summ$r.squared, 4),
    adj_r_squared = round(summ$adj.r.squared, 4),
    importance   = round(importance, 4),
    n_obs        = nrow(df_clean),
    predictors   = predictors
  )
}
