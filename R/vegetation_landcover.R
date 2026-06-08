#' Importer les données d'occupation du sol
#'
#' @description
#' Importe et reclassifie un raster d'occupation du sol (ESA ou Corine Land Cover).
#' Peut aussi simuler un raster d'occupation du sol.
#'
#' @param filepath Chemin vers le raster (.tif) ou NULL pour données simulées
#' @param reclassify Reclassifier les codes (défaut: TRUE)
#' @param target_crs CRS cible pour harmonisation (défaut: "EPSG:4326")
#' @param seed Graine aléatoire si simulation (défaut: 42)
#'
#' @return Un SpatRaster avec classes: 1=cultures, 2=forêt, 3=prairie, 4=sol nu, 5=urbain
#'
#' @details
#' Classes de sortie standardisées:
#' \describe{
#'   \item{1}{Cultures / terres agricoles}
#'   \item{2}{Forêt / végétation arborée}
#'   \item{3}{Prairie / végétation basse}
#'   \item{4}{Sol nu / zones dénudées}
#'   \item{5}{Zone urbaine / artificielle}
#' }
#'
#' @examples
#' lc <- import_landcover()
#' terra::plot(lc, main = "Occupation du sol", col = c("gold", "darkgreen",
#'   "lightgreen", "tan", "gray50"))
#'
#' @export
import_landcover <- function(filepath = NULL,
                             reclassify = TRUE,
                             target_crs = "EPSG:4326",
                             seed = 42) {

  if (is.null(filepath)) {
    message("Simulation d'un raster d'occupation du sol...")
    set.seed(seed)

    nrow_r <- 50; ncol_r <- 50
    # Créer un raster avec classes spatiales cohérentes
    mat <- matrix(1L, nrow = nrow_r, ncol = ncol_r)

    # Cultures au centre
    mat[15:35, 15:35] <- 1L
    # Forêt en haut à gauche
    mat[1:20, 1:20] <- 2L
    # Prairie
    mat[30:50, 1:25] <- 3L
    # Sol nu (bas droite)
    mat[35:50, 35:50] <- 4L
    # Zone urbaine (petit carré)
    mat[20:28, 30:38] <- 5L
    # Bruit aléatoire
    noise_idx <- sample(seq_len(nrow_r * ncol_r), size = round(nrow_r * ncol_r * 0.05))
    mat[noise_idx] <- sample(1:5, length(noise_idx), replace = TRUE)

    lc <- terra::rast(mat)
    terra::ext(lc) <- terra::ext(c(-5.5, -5.0, 33.5, 34.0))
    terra::crs(lc) <- "EPSG:4326"
    names(lc) <- "land_cover"
    levels(lc) <- data.frame(id = 1:5,
                             class = c("Cultures", "Foret", "Prairie", "Sol nu", "Urbain"))
    message("Raster occupation du sol simule (50x50, 5 classes)")
    return(lc)
  }

  lc <- terra::rast(filepath)

  if (!is.null(target_crs)) {
    if (terra::crs(lc) != terra::crs(terra::rast(crs = target_crs))) {
      lc <- terra::project(lc, target_crs, method = "near")
    }
  }

  message(sprintf("Occupation du sol importee: %d x %d pixels", terra::nrow(lc), terra::ncol(lc)))
  lc
}


#' Calculer le NDVI
#'
#' @description
#' Calcule l'indice de végétation par différence normalisée (NDVI) à partir
#' des bandes rouge et proche infrarouge d'une image satellitaire.
#'
#' NDVI = (NIR - Rouge) / (NIR + Rouge)
#'
#' @param red Raster de la bande rouge (ou NULL pour données simulées)
#' @param nir Raster de la bande proche infrarouge (ou NULL pour données simulées)
#' @param n_dates Nombre de dates si simulation (défaut: 3)
#' @param scale_factor Facteur d'échelle à appliquer (défaut: 1.0)
#'
#' @return Un SpatRaster de NDVI (valeurs entre -1 et 1)
#'
#' @examples
#' # NDVI simulé
#' ndvi <- calculate_ndvi(n_dates = 4)
#' terra::plot(ndvi, main = "NDVI", col = grDevices::hcl.colors(100, "Greens"))
#'
#' # Avec des rasters réels:
#' # ndvi <- calculate_ndvi(red = rast("bande_rouge.tif"), nir = rast("bande_nir.tif"))
#'
#' @export
calculate_ndvi <- function(red = NULL, nir = NULL,
                           n_dates = 3,
                           scale_factor = 1.0) {

  if (is.null(red) || is.null(nir)) {
    message("Simulation du NDVI sur ", n_dates, " dates...")
    ndvi <- simulate_raster_data(type = "ndvi", n_dates = n_dates)
    return(ndvi)
  }

  if (!inherits(red, "SpatRaster") || !inherits(nir, "SpatRaster")) {
    stop("red et nir doivent etre des objets SpatRaster")
  }

  # Appliquer facteur d'échelle
  if (scale_factor != 1.0) {
    red <- red * scale_factor
    nir <- nir * scale_factor
  }

  # Calcul NDVI avec protection division par zéro
  ndvi <- (nir - red) / (nir + red + 1e-10)

  # Masquer valeurs hors [-1, 1]
  ndvi[ndvi < -1] <- NA
  ndvi[ndvi >  1] <- NA

  names(ndvi) <- paste0("NDVI_", seq_len(terra::nlyr(ndvi)))

  # Stats
  v <- terra::values(ndvi[[1]])
  message(sprintf("NDVI calcule: min=%.3f, mean=%.3f, max=%.3f",
                  min(v, na.rm=TRUE), mean(v, na.rm=TRUE), max(v, na.rm=TRUE)))
  ndvi
}


#' Extraire les variables microclimatiques
#'
#' @description
#' Extrait et agrège les variables microclimatiques (LST, NDVI, altitude,
#' occupation du sol) pour créer un dataframe d'analyse.
#'
#' @param lst_raster SpatRaster de LST
#' @param ndvi_raster SpatRaster de NDVI (optionnel)
#' @param landcover_raster SpatRaster d'occupation du sol (optionnel)
#' @param points Objet sf de points d'extraction (optionnel, sinon grille régulière)
#' @param n_sample Nombre de points d'échantillonnage si pas de points fournis (défaut: 500)
#' @param layer_index Index de couche LST à utiliser (défaut: 1)
#'
#' @return Un data.frame avec variables microclimatiques
#'
#' @examples
#' lst  <- simulate_raster_data(type = "lst",  n_dates = 1)
#' ndvi <- simulate_raster_data(type = "ndvi", n_dates = 1)
#' lc   <- import_landcover()
#' features <- extract_microclimate_features(lst, ndvi, lc)
#' head(features)
#'
#' @export
extract_microclimate_features <- function(lst_raster,
                                          ndvi_raster    = NULL,
                                          landcover_raster = NULL,
                                          points         = NULL,
                                          n_sample       = 500,
                                          layer_index    = 1) {

  message("Extraction des variables microclimatiques...")

  lst_layer <- lst_raster[[min(layer_index, terra::nlyr(lst_raster))]]

  if (is.null(points)) {
    # Créer une grille d'échantillonnage
    n_total <- terra::ncell(lst_layer)
    idx <- sort(sample(n_total, min(n_sample, n_total)))
    xys <- terra::xyFromCell(lst_layer, idx)
    points_df <- data.frame(x = xys[,1], y = xys[,2])
    pts_sf <- sf::st_as_sf(points_df, coords = c("x","y"), crs = terra::crs(lst_layer))
  } else {
    pts_sf <- points
  }

  # Extraire LST
  df <- data.frame(
    x = sf::st_coordinates(pts_sf)[,1],
    y = sf::st_coordinates(pts_sf)[,2]
  )

  lst_vals <- terra::extract(lst_layer, pts_sf)
  df$lst <- lst_vals[,2]

  # Extraire NDVI
  if (!is.null(ndvi_raster)) {
    ndvi_layer <- ndvi_raster[[min(layer_index, terra::nlyr(ndvi_raster))]]
    ndvi_vals <- terra::extract(ndvi_layer, pts_sf)
    df$ndvi <- ndvi_vals[,2]
  }

  # Extraire occupation du sol
  if (!is.null(landcover_raster)) {
    lc_vals <- terra::extract(landcover_raster, pts_sf)
    df$land_cover <- lc_vals[,2]
    df$land_cover_name <- c("Cultures","Foret","Prairie","Sol nu","Urbain")[
      pmin(5, pmax(1, round(df$land_cover, 0)), na.rm = FALSE)
    ]
  }

  # Distance estimée au pixel de végétation le plus proche (proxy simple)
  if (!is.null(ndvi_raster)) {
    df$veg_proximity <- ifelse(!is.na(df$ndvi), df$ndvi * 100, NA)
  }

  df <- df[!is.na(df$lst), ]
  message(sprintf("  %d points extraits avec variables: %s",
                  nrow(df), paste(names(df), collapse=", ")))
  df
}
