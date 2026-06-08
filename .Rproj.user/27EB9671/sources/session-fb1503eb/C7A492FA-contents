#' Importer des données météorologiques
#'
#' @description
#' Importe des données météorologiques depuis un fichier CSV ou Excel.
#' Gère les dates, la fusion multi-stations et les données manquantes.
#'
#' @param filepath Chemin vers le fichier (CSV ou Excel .xlsx/.xls)
#' @param date_col Nom de la colonne date (défaut: "date")
#' @param temp_col Nom de la colonne température (défaut: "temperature")
#' @param hum_col Nom de la colonne humidité (défaut: "humidity")
#' @param wind_col Nom de la colonne vent (défaut: "wind_speed")
#' @param rad_col Nom de la colonne rayonnement (défaut: "radiation")
#' @param station_col Nom de la colonne station (défaut: "station_id")
#' @param fill_missing Méthode pour combler les NA: "mean", "linear", "none"
#'   (défaut: "mean")
#'
#' @return Un data.frame avec les données météo standardisées
#'
#' @examples
#' # Utiliser des données simulées (sans fichier réel)
#' weather <- simulate_weather_data(n_days = 30, n_stations = 3)
#'
#' # Ou depuis un fichier CSV:
#' # weather <- import_weather_data("mes_donnees.csv")
#'
#' @export
import_weather_data <- function(filepath = NULL,
                                date_col = "date",
                                temp_col = "temperature",
                                hum_col  = "humidity",
                                wind_col = "wind_speed",
                                rad_col  = "radiation",
                                station_col = "station_id",
                                fill_missing = "mean") {

  # Si pas de fichier, utiliser données simulées
  if (is.null(filepath)) {
    message("Aucun fichier fourni. Utilisation de donnees simulees.")
    df <- simulate_weather_data()
    return(df)
  }

  # Détecter le format
  ext <- tolower(tools::file_ext(filepath))

  if (ext == "csv") {
    df <- utils::read.csv(filepath, stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(filepath)
    df <- as.data.frame(df)
  } else {
    stop("Format non supporte. Utilisez CSV ou Excel (.xlsx, .xls)")
  }

  # Renommer colonnes si trouvées
  col_map <- c(date_col, temp_col, hum_col, wind_col, rad_col, station_col)
  std_names <- c("date", "temperature", "humidity", "wind_speed", "radiation", "station_id")

  for (i in seq_along(col_map)) {
    if (col_map[i] %in% names(df)) {
      names(df)[names(df) == col_map[i]] <- std_names[i]
    }
  }

  # Convertir dates
  if ("date" %in% names(df)) {
    df$date <- tryCatch(
      as.Date(df$date),
      error = function(e) {
        warning("Conversion date echouee, tentative format lubridate...")
        lubridate::ymd(df$date)
      }
    )
  }

  # Gestion des valeurs manquantes
  num_cols <- intersect(c("temperature", "humidity", "wind_speed", "radiation"), names(df))

  if (fill_missing == "mean") {
    for (col in num_cols) {
      if (any(is.na(df[[col]]))) {
        col_mean <- mean(df[[col]], na.rm = TRUE)
        df[[col]][is.na(df[[col]])] <- col_mean
        message(sprintf("NA remplaces par la moyenne dans '%s'", col))
      }
    }
  } else if (fill_missing == "linear") {
    for (col in num_cols) {
      if (any(is.na(df[[col]]))) {
        df[[col]] <- stats::approx(seq_len(nrow(df)), df[[col]],
                                   xout = seq_len(nrow(df)))$y
      }
    }
  }

  # Statistiques de base
  message(sprintf("Donnees importees: %d lignes, %d stations, periode: %s a %s",
                  nrow(df),
                  ifelse("station_id" %in% names(df), length(unique(df$station_id)), 1),
                  ifelse("date" %in% names(df), as.character(min(df$date, na.rm=TRUE)), "?"),
                  ifelse("date" %in% names(df), as.character(max(df$date, na.rm=TRUE)), "?")))

  df
}


#' Importer les localisations des stations météorologiques
#'
#' @description
#' Importe les localisations des stations météo depuis un shapefile ou un CSV
#' avec coordonnées GPS. Convertit en objet spatial sf.
#'
#' @param filepath Chemin vers le fichier (shapefile .shp ou CSV avec lon/lat)
#' @param lon_col Nom de la colonne longitude (pour CSV, défaut: "longitude")
#' @param lat_col Nom de la colonne latitude (pour CSV, défaut: "latitude")
#' @param id_col Nom de la colonne identifiant (défaut: "station_id")
#' @param crs Code CRS (défaut: 4326 = WGS84)
#' @param plot Afficher une carte des stations (défaut: TRUE)
#'
#' @return Un objet sf avec les localisations des stations
#'
#' @examples
#' # Créer des stations simulées
#' stations_df <- data.frame(
#'   station_id = paste0("ST", 1:5),
#'   longitude  = c(-5.2, -5.1, -5.3, -5.0, -5.4),
#'   latitude   = c(33.7, 33.8, 33.6, 33.9, 33.75),
#'   nom        = paste("Station", 1:5)
#' )
#' stations_sf <- import_station_locations(stations_df)
#'
#' @export
import_station_locations <- function(filepath = NULL,
                                     lon_col = "longitude",
                                     lat_col = "latitude",
                                     id_col  = "station_id",
                                     crs     = 4326,
                                     plot    = TRUE) {

  # Si data.frame directement passé (pour usage en mémoire)
  if (is.data.frame(filepath)) {
    df <- filepath
  } else if (is.null(filepath)) {
    message("Aucun fichier fourni. Creation de stations simulees.")
    df <- data.frame(
      station_id = paste0("ST", 1:5),
      longitude  = c(-5.20, -5.10, -5.30, -5.05, -5.40),
      latitude   = c(33.70, 33.80, 33.65, 33.90, 33.75),
      nom        = paste("Station", 1:5)
    )
  } else {
    ext <- tolower(tools::file_ext(filepath))

    if (ext == "shp") {
      sf_obj <- sf::st_read(filepath, quiet = TRUE)
      if (sf::st_crs(sf_obj)$epsg != crs) {
        sf_obj <- sf::st_transform(sf_obj, crs)
      }
      if (plot) {
        graphics::plot(sf_obj["geometry"],
                       main = "Stations meteorologiques",
                       pch = 16, col = "steelblue")
      }
      return(sf_obj)
    } else if (ext == "csv") {
      df <- utils::read.csv(filepath, stringsAsFactors = FALSE)
    } else {
      stop("Format non supporte. Utilisez .shp ou .csv")
    }
  }

  # Renommer colonnes
  if (lon_col %in% names(df)) names(df)[names(df) == lon_col] <- "longitude"
  if (lat_col %in% names(df)) names(df)[names(df) == lat_col] <- "latitude"

  # Vérifications
  if (!all(c("longitude", "latitude") %in% names(df))) {
    stop("Colonnes longitude/latitude introuvables. Verifiez les noms de colonnes.")
  }

  if (any(abs(df$longitude) > 180) || any(abs(df$latitude) > 90)) {
    warning("Coordonnees hors limites WGS84. Verifiez le CRS.")
  }

  # Convertir en objet sf
  sf_obj <- sf::st_as_sf(df,
                         coords = c("longitude", "latitude"),
                         crs = crs)

  message(sprintf("%d stations importees (CRS: EPSG:%d)", nrow(sf_obj), crs))

  # Afficher carte simple
  if (plot) {
    graphics::plot(sf::st_geometry(sf_obj),
                   main = "Stations meteorologiques",
                   pch = 16, col = "steelblue", cex = 1.5)
    if (id_col %in% names(sf_obj)) {
      graphics::text(sf::st_coordinates(sf_obj),
                     labels = sf_obj[[id_col]],
                     pos = 3, cex = 0.7)
    }
  }

  sf_obj
}
