#' Comparer LST et NDVI
#'
#' @description
#' Analyse la relation entre la température de surface (LST) et le NDVI.
#' Produit un scatterplot avec régression et calcule le coefficient de corrélation.
#'
#' @param lst_raster SpatRaster de LST
#' @param ndvi_raster SpatRaster de NDVI
#' @param layer_index Index de couche à utiliser (défaut: 1)
#' @param n_sample Nombre de pixels à échantillonner (défaut: 1000)
#' @param plot Afficher le graphique (défaut: TRUE)
#'
#' @return Une liste avec: correlation, p_value, plot (ggplot), data
#'
#' @examples
#' lst  <- simulate_raster_data(type = "lst",  n_dates = 3)
#' ndvi <- simulate_raster_data(type = "ndvi", n_dates = 3)
#' result <- compare_lst_vs_ndvi(lst, ndvi, layer_index = 1)
#' cat("Correlation LST-NDVI:", result$correlation, "\n")
#'
#' @export
compare_lst_vs_ndvi <- function(lst_raster,
                                ndvi_raster,
                                layer_index = 1,
                                n_sample    = 1000,
                                plot        = TRUE) {

  lst_layer  <- lst_raster[[min(layer_index, terra::nlyr(lst_raster))]]
  ndvi_layer <- ndvi_raster[[min(layer_index, terra::nlyr(ndvi_raster))]]

  # Harmoniser les résolutions si nécessaire
  if (!identical(terra::ext(lst_layer), terra::ext(ndvi_layer)) ||
      terra::ncell(lst_layer) != terra::ncell(ndvi_layer)) {
    message("Resolutions differentes, rééchantillonnage du NDVI sur LST...")
    ndvi_layer <- terra::resample(ndvi_layer, lst_layer, method = "bilinear")
  }

  lst_vals  <- as.vector(terra::values(lst_layer))
  ndvi_vals <- as.vector(terra::values(ndvi_layer))

  valid <- !is.na(lst_vals) & !is.na(ndvi_vals)
  lst_vals  <- lst_vals[valid]
  ndvi_vals <- ndvi_vals[valid]

  if (length(lst_vals) > n_sample) {
    idx <- sample(length(lst_vals), n_sample)
    lst_vals  <- lst_vals[idx]
    ndvi_vals <- ndvi_vals[idx]
  }

  # Corrélation
  cor_test <- stats::cor.test(ndvi_vals, lst_vals)
  r <- round(cor_test$estimate, 4)
  p <- round(cor_test$p.value, 4)

  message(sprintf("Correlation NDVI-LST: r = %.4f, p = %.4f", r, p))
  interp <- if (r < -0.5) "Forte relation negative: vegetation reduit la chaleur"
  else if (r < -0.2) "Relation negative moderee: vegetation tend a refroidir"
  else if (r > 0.5) "Forte relation positive: zones chaudes plus vegetalisees"
  else "Faible correlation LST-NDVI"
  message(sprintf("  %s", interp))

  df_plot <- data.frame(lst = lst_vals, ndvi = ndvi_vals)

  # Graphique ggplot2
  p_obj <- ggplot2::ggplot(df_plot, ggplot2::aes(x = ndvi, y = lst)) +
    ggplot2::geom_point(alpha = 0.3, color = "#1D9E75", size = 1.2) +
    ggplot2::geom_smooth(method = "lm", color = "#993C1D", linewidth = 1.2, se = TRUE) +
    ggplot2::labs(
      title    = "Relation NDVI - Temperature de surface (LST)",
      subtitle = sprintf("r = %.3f (p = %.4f) | n = %d pixels", r, p, length(lst_vals)),
      x        = "NDVI",
      y        = "LST (°C)",
      caption  = "microclimAg"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (plot) print(p_obj)

  list(
    correlation    = r,
    p_value        = p,
    interpretation = interp,
    plot           = p_obj,
    data           = df_plot
  )
}


#' Cartographier le microclimat
#'
#' @description
#' Crée des cartes thématiques du microclimat: LST, anomalies thermiques,
#' îlots de chaleur, NDVI, zones de clustering.
#'
#' @param raster_data SpatRaster à cartographier (LST, NDVI, anomalies, clusters)
#' @param type Type de carte: "lst", "ndvi", "anomaly", "cluster" (défaut: "lst")
#' @param layer_index Index de couche (défaut: 1)
#' @param title Titre de la carte (défaut: NULL = automatique)
#' @param palette Palette de couleurs (défaut: NULL = automatique selon type)
#' @param output_file Chemin pour sauvegarder (PNG/PDF, défaut: NULL)
#' @param facet Afficher toutes les couches (défaut: FALSE)
#'
#' @return Un objet ggplot2 invisible
#'
#' @examples
#' lst <- simulate_raster_data(type = "lst", n_dates = 3)
#' p   <- plot_microclimate_map(lst, type = "lst")
#'
#' ndvi <- simulate_raster_data(type = "ndvi", n_dates = 1)
#' p2   <- plot_microclimate_map(ndvi, type = "ndvi")
#'
#' @export
plot_microclimate_map <- function(raster_data,
                                  type        = "lst",
                                  layer_index = 1,
                                  title       = NULL,
                                  palette     = NULL,
                                  output_file = NULL,
                                  facet       = FALSE) {

  type <- match.arg(type, c("lst", "ndvi", "anomaly", "cluster"))

  # Palettes selon type
  palettes <- list(
    lst     = c("#313695","#4575b4","#74add1","#abd9e9","#fee090","#fdae61","#f46d43","#d73027","#a50026"),
    ndvi    = c("#d73027","#f46d43","#fee08b","#d9ef8b","#66bd63","#1a9850","#006837"),
    anomaly = c("#313695","#4575b4","#74add1","#ffffff","#fdae61","#f46d43","#d73027"),
    cluster = c("#4575b4","#fee090","#d73027","#74add1","#a50026","#313695")
  )
  if (is.null(palette)) palette <- palettes[[type]]

  # Titres automatiques
  auto_titles <- c(
    lst     = "Temperature de surface - LST (°C)",
    ndvi    = "Indice de vegetation NDVI",
    anomaly = "Anomalie thermique LST",
    cluster = "Zones thermiques (K-means)"
  )
  if (is.null(title)) title <- auto_titles[type]

  # Préparer les données pour ggplot
  if (facet && terra::nlyr(raster_data) > 1) {
    # Toutes les couches avec facet
    df_list <- list()
    for (i in seq_len(min(terra::nlyr(raster_data), 6))) {
      layer <- raster_data[[i]]
      df_i <- as.data.frame(layer, xy = TRUE)
      names(df_i)[3] <- "value"
      df_i$layer <- names(raster_data)[i]
      df_list[[i]] <- df_i
    }
    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colors = palette, na.value = "white") +
      ggplot2::facet_wrap(~layer) +
      ggplot2::labs(title = title, x = "Longitude", y = "Latitude", fill = type) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     axis.text  = ggplot2::element_text(size = 7))

  } else {
    layer <- raster_data[[min(layer_index, terra::nlyr(raster_data))]]
    df <- as.data.frame(layer, xy = TRUE)
    names(df)[3] <- "value"

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colors = palette, na.value = "white") +
      ggplot2::labs(title = title, x = "Longitude", y = "Latitude",
                    fill = ifelse(type == "lst", "LST (°C)",
                                  ifelse(type == "ndvi", "NDVI",
                                         ifelse(type == "anomaly", "Anomalie", "Cluster"))),
                    caption = "microclimAg") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     aspect.ratio = 1)
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = 8, height = 6, dpi = 150)
    message(sprintf("Carte sauvegardee: %s", output_file))
  }

  print(p)
  invisible(p)
}


#' Générer un rapport HTML automatique
#'
#' @description
#' Génère un rapport HTML ou PDF complet avec les statistiques thermiques,
#' cartes LST, anomalies, comparaison NDVI et interprétation microclimatique.
#'
#' @param lst_raster SpatRaster de LST
#' @param ndvi_raster SpatRaster de NDVI (optionnel)
#' @param weather_data Data.frame de données météo (optionnel)
#' @param output_file Chemin du fichier de sortie (défaut: "rapport_microclimAg.html")
#' @param title Titre du rapport (défaut: "Rapport MicroclimAg")
#' @param open_report Ouvrir automatiquement le rapport (défaut: TRUE)
#'
#' @return Le chemin du rapport généré
#'
#' @examples
#' \dontrun{
#' lst  <- simulate_raster_data(type = "lst",  n_dates = 3)
#' ndvi <- simulate_raster_data(type = "ndvi", n_dates = 3)
#' generate_report(lst, ndvi, output_file = "mon_rapport.html")
#' }
#'
#' @export
generate_report <- function(lst_raster,
                            ndvi_raster  = NULL,
                            weather_data = NULL,
                            output_file  = "rapport_microclimAg.html",
                            title        = "Rapport MicroclimAg",
                            open_report  = TRUE) {

  message("Generation du rapport...")

  # Calculer statistiques LST
  lst_vals <- as.vector(terra::values(lst_raster[[1]]))
  lst_vals <- lst_vals[!is.na(lst_vals)]

  lst_stats <- list(
    mean = round(mean(lst_vals), 2),
    sd   = round(stats::sd(lst_vals), 2),
    min  = round(min(lst_vals), 2),
    max  = round(max(lst_vals), 2),
    n_dates = terra::nlyr(lst_raster)
  )

  # Clustering
  cluster_result <- cluster_thermal_zones(lst_raster, n_clusters = 3)

  # Corrélation NDVI-LST
  cor_result <- NULL
  if (!is.null(ndvi_raster)) {
    cor_result <- compare_lst_vs_ndvi(lst_raster, ndvi_raster, plot = FALSE)
  }

  # Anomalies
  anom_result <- NULL
  if (terra::nlyr(lst_raster) >= 2) {
    anom_result <- calculate_lst_anomaly(lst_raster, method = "zscore")
  }

  # Construire HTML
  html_content <- sprintf('<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8">
  <title>%s</title>
  <style>
    body { font-family: Arial, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; }
    h1 { color: #1D9E75; border-bottom: 3px solid #1D9E75; padding-bottom: 10px; }
    h2 { color: #0F6E56; margin-top: 30px; }
    .stat-grid { display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin: 20px 0; }
    .stat-box { background: #f5f5f5; border-radius: 8px; padding: 15px; text-align: center; }
    .stat-val { font-size: 2em; font-weight: bold; color: #1D9E75; }
    .stat-lab { font-size: 0.85em; color: #666; }
    table { width: 100%%; border-collapse: collapse; }
    th { background: #1D9E75; color: white; padding: 8px 12px; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee; }
    .interp { background: #E1F5EE; border-left: 4px solid #1D9E75; padding: 12px 16px; margin: 15px 0; border-radius: 4px; }
    .footer { margin-top: 50px; color: #999; font-size: 0.85em; border-top: 1px solid #eee; padding-top: 15px; }
  </style>
</head>
<body>
  <h1>%s</h1>
  "<p>Package <strong>microclimAg</strong> &mdash; Rapport g\u00e9n\u00e9r\u00e9 le %s</p>"

  <h2>1. Statistiques LST</h2>
  <div class="stat-grid">
    <div class="stat-box"><div class="stat-val">%.1f°C</div><div class="stat-lab">Température moyenne</div></div>
    <div class="stat-box"><div class="stat-val">%.1f°C</div><div class="stat-lab">Min</div></div>
    <div class="stat-box"><div class="stat-val">%.1f°C</div><div class="stat-lab">Max</div></div>
    <div class="stat-box"><div class="stat-val">%d</div><div class=\"stat-lab\">Dates analys\u00e9es</div>"
  </div>

  <h2>2. Zones thermiques (K-means, k=3)</h2>
  <table>
    <tr><th>Cluster</th><th>Classe</th><th>Temp. moyenne (°C)</th><th>Écart-type</th><th>Pixels</th></tr>
    %s
  </table>

  %s

  %s

  <div class="footer">
    <p>Rapport généré avec <strong>microclimAg</strong> v0.1.0 | IAV Hassan II</p>
  </div>
</body>
</html>',
                          title, title,
                          format(Sys.Date(), "%d/%m/%Y"),
                          lst_stats$mean, lst_stats$min, lst_stats$max, lst_stats$n_dates,
                          paste(apply(cluster_result$stats, 1, function(r) {
                            sprintf('<tr><td>%s</td><td>%s</td><td>%.2f</td><td>%.2f</td><td>%s</td></tr>',
                                    r["cluster"], r["label"], as.numeric(r["temp_mean"]),
                                    as.numeric(r["temp_sd"]), r["n_pixels"])
                          }), collapse = "\n"),
                          if (!is.null(cor_result)) {
                            sprintf('<h2>3. Relation NDVI - LST</h2><div class="interp"><strong>r = %.4f</strong> (p = %.4f)<br>%s</div>',
                                    cor_result$correlation, cor_result$p_value, cor_result$interpretation)
                          } else "",
                          if (!is.null(anom_result)) {
                            anom_vals <- as.vector(terra::values(anom_result[[1]]))
                            anom_vals <- anom_vals[!is.na(anom_vals)]
                            sprintf('<h2>4. Anomalies thermiques</h2><div class="interp">Anomalies Z-score — écart-type: %.2f | min: %.2f | max: %.2f</div>',
                                    stats::sd(anom_vals), min(anom_vals), max(anom_vals))
                          } else ""
  )

  writeLines(html_content, output_file)
  message(sprintf("Rapport HTML genere: %s", output_file))

  if (open_report && interactive()) {
    utils::browseURL(output_file)
  }

  invisible(output_file)
}
