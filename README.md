# microclimAg

> Analyse du microclimat agricole et des îlots de chaleur ruraux — IAV Hassan II

## Présentation

**microclimAg** est un package R développé à l'Institut Agronomique et Vétérinaire Hassan II de Rabat, Maroc. Il permet d'analyser les variations microclimatiques dans les paysages agricoles en combinant télédétection thermique (LST), indices de végétation (NDVI), données météorologiques et machine learning spatial.

## Fonctionnalités

- Import et nettoyage de données LST (MODIS/Landsat)
- Calcul du NDVI et import de l'occupation du sol
- Détection des îlots de chaleur ruraux par K-means
- Analyse de l'autocorrélation spatiale (Moran's I)
- Interpolation spatiale des températures (IDW)
- Modélisation LST ~ NDVI par régression linéaire
- Cartographie thématique et rapport HTML automatique

## Installation

```r
devtools::install_github("rihamdahnoun-droid/microclimAg")
```

## Utilisation rapide

```r
library(microclimAg)

lst    <- download_lst_data("2023-06-01", "2023-08-31",
                             c(-5.5, -5.0, 33.5, 34.0))
lst    <- clean_lst_data(lst)
ndvi   <- calculate_ndvi(n_dates = 3)
zones  <- cluster_thermal_zones(lst, n_clusters = 3)

generate_report(lst, ndvi, output_file = "rapport.html")
```

## Tests

```r
load_all()
test()
# FAIL 0 | WARN 0 | SKIP 0 | PASS 45
```
## Structure

```microclimAg/
├── R/                  # 6 fichiers, 17 fonctions
├── tests/testthat/     # 45 tests unitaires
├── vignettes/          # Workflow complet reproductible
├── man/                # Documentation roxygen2
├── DESCRIPTION
└── NAMESPACE
```
## Auteur

**Riham Dahnoun** — IAV Hassan II, Rabat, Maroc

## Licence

MIT © 2024 — IAV Hassan II
```
