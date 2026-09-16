# ==============================================================================
# Selection of sub-clustering resolution by stability and separation criteria
# Master's Thesis in Bioinformatics | Victor Daniel Velasco Perez
#
# The script sweeps a grid of clustering resolutions within each major lineage
# and quantifies three properties that jointly determine whether a partition is
# defensible, replacing visual inspection of marker plots as the selection
# criterion.
#
#   Stability          Cluster membership is tracked across resolutions.
#   Separation         Mean silhouette width computed over PCA space.
#   Discriminative     Marker genes surviving filters are counted per cluster.
#   Coherence          Screening for simultaneous detection of exclusive markers.
#
# References
#   Zappia L, Oshlack A. Gigascience. 2018;7(7):giy083.
#   Rousseeuw PJ. J Comput Appl Math. 1987;20:53-65.
# ==============================================================================

graphics.off()
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(clustree)
  library(cluster)
})

# ==============================================================================
# 0. SETUP & GLOBAL CONSTANTS
# ==============================================================================
base_dir <- "C:/Master en Bioinformatica/TFM"
setwd(base_dir)

etapa   <- "etapa_2"
in_dir  <- file.path("resultados", etapa, "datos")
out_dir <- file.path("resultados", etapa, "diagnostico_resolucion")

invisible(lapply(c("graficos", "datos"), function(d) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}))

resoluciones             <- seq(0.1, 1.0, by = 0.1)
umbral_deteccion_ambigua <- 0.5

config_linajes <- list(
  list(
    linaje    = "T_cell", 
    etiqueta  = "Célula T", 
    nfeatures = 1500, 
    pcs       = 15,
    pares     = list(c("FOXP3", "CXCL13"), c("CCR7", "GZMB"), c("CD4", "CD8A")),
    panel     = c("CD4", "CD8A", "CCR7", "SELL", "FOXP3", "IL2RA", "CTLA4", "IKZF2",
                  "CXCL13", "PDCD1", "IL21", "GZMB", "PRF1", "IFNG", "GZMK")
  ),
  list(
    linaje    = "B_cell", 
    etiqueta  = "Célula B", 
    nfeatures = 500,  
    pcs       = 8,
    pares     = list(c("IGHD", "PRDM1"), c("MS4A1", "SDC1")),
    panel     = c("MS4A1", "CD19", "IGHD", "SELL", "CR2", "CD27", "CD38", "TBX21",
                  "ITGAX", "SDC1", "XBP1", "PRDM1", "IGJ", "MZB1")
  ),
  list(
    linaje    = "Fibroblast", 
    etiqueta  = "Fibroblasto", 
    nfeatures = 1000, 
    pcs       = 12,
    pares     = list(c("FAP", "PRG4"), c("THY1", "PRG4")),
    panel     = c("PRG4", "CLIC5", "HTRA4", "FAP", "THY1", "PDPN", "CDH11", "PTGFR",
                  "CXCL12", "CCL2", "IL6", "COL1A1", "DKK3", "SFRP1")
  ),
  list(
    linaje    = "Monocyte", 
    etiqueta  = "Macrófago", 
    nfeatures = 800,  
    pcs       = 10,
    pares     = list(c("FCN1", "C1QA"), c("CXCL9", "LYVE1")),
    panel     = c("CD14", "FCN1", "VCAN", "S100A8", "CXCL9", "CXCL10", "IL1B", "SPP1",
                  "C1QA", "CD163", "MARCO", "LYVE1", "MRC1", "TREM2")
  )
)

tema_figura <- theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30")
  )

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

preparar_linaje <- function(seurat_obj, linaje, nfeatures, pcs) {
  suppressWarnings({
    sub <- subset(seurat_obj, type == linaje)
    sub <- NormalizeData(sub, verbose = FALSE) %>%
      FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(seed.use = 42, verbose = FALSE) %>%
      FindNeighbors(dims = 1:pcs, verbose = FALSE)
    
    for (r in resoluciones) {
      sub <- FindClusters(sub, resolution = r, random.seed = 0, verbose = FALSE)
    }
  })
  return(sub)
}

ancho_silueta <- function(obj, columna, pcs) {
  cl <- as.integer(as.character(obj@meta.data[[columna]]))
  if (length(unique(cl)) < 2) return(NA_real_)
  
  emb <- Embeddings(obj, reduction = "pca")[, 1:pcs]
  sil <- silhouette(cl, dist(emb))
  mean(sil[, 3])
}

contar_marcadores <- function(obj, columna) {
  Idents(obj) <- columna
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25, verbose = FALSE)
  
  if (!"avg_log2FC" %in% colnames(markers)) {
    return(c(minimo = 0L, mediana = 0L))
  }
  
  markers <- markers %>% filter(avg_log2FC > 1 & pct.1 > 0.5 & pct.2 < 0.25)
  n_por_cluster <- table(factor(markers$cluster, levels = levels(Idents(obj))))
  c(minimo = min(n_por_cluster), mediana = median(n_por_cluster))
}

clusters_ambiguos <- function(obj, columna, pares, umbral = umbral_deteccion_ambigua) {
  presentes <- Filter(function(p) all(p %in% rownames(obj)), pares)
  if (length(presentes) == 0) return(NA_integer_)
  
  cl   <- obj@meta.data[[columna]]
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  
  ambiguos <- vapply(levels(factor(cl)), function(k) {
    celulas <- which(cl == k)
    any(vapply(presentes, function(p) {
      min(rowMeans(expr[p, celulas, drop = FALSE] > 0)) > umbral
    }, logical(1)))
  }, logical(1))
  
  sum(ambiguos)
}

# ==============================================================================
# 2. RESOLUTION SWEEP & DIAGNOSTICS
# ==============================================================================

amp_seurat <- readRDS(file.path(in_dir, "amp_seurat_normalized.rds"))

cat("\n[+] Sweeping resolutions", min(resoluciones), "to", max(resoluciones),
    "in", length(config_linajes), "lineages...\n")

resumen_global <- list()

for (cfg in config_linajes) {
  cat("    ", cfg$linaje, "\n")
  sub      <- preparar_linaje(amp_seurat, cfg$linaje, cfg$nfeatures, cfg$pcs)
  columnas <- paste0("RNA_snn_res.", resoluciones)
  
  metricas <- lapply(seq_along(resoluciones), function(i) {
    col        <- columnas[i]
    marcadores <- contar_marcadores(sub, col)
    n_cl       <- length(unique(sub@meta.data[[col]]))
    n_amb      <- clusters_ambiguos(sub, col, cfg$pares)
    
    data.frame(
      Linaje             = cfg$etiqueta,
      Resolucion         = resoluciones[i],
      N_Clusters         = n_cl,
      Ancho_Silueta      = round(ancho_silueta(sub, col, cfg$pcs), 4),
      Marcadores_Min     = as.integer(marcadores["minimo"]),
      Marcadores_Mediana = as.integer(marcadores["mediana"]),
      N_Ambiguos         = n_amb,
      Prop_Ambiguos      = round(n_amb / n_cl, 4)
    )
  })
  
  metricas <- bind_rows(metricas)
  resumen_global[[cfg$linaje]] <- metricas
  
  p_arbol <- clustree(sub, prefix = "RNA_snn_res.") +
    labs(title = cfg$etiqueta, subtitle = "Trasvase de células entre particiones consecutivas")
  
  p_silueta <- ggplot(metricas, aes(Resolucion, Ancho_Silueta)) +
    geom_line(colour = "grey50") +
    geom_point(size = 2) +
    tema_figura +
    labs(title = "Separación entre Subclusters",
         subtitle = "Anchura media de silueta sobre el espacio de componentes principales",
         x = "Resolución", y = "Anchura Media de Silueta")
  
  p_marcadores <- ggplot(metricas, aes(Resolucion, Marcadores_Min)) +
    geom_col(fill = "#2c7fb8") +
    tema_figura +
    labs(title = "Capacidad Discriminativa del Subcluster Más Débil",
         subtitle = "Genes marcadores que superan los filtros de la firma",
         x = "Resolución", y = "Número de Genes Marcadores")
  
  p_ambiguos <- ggplot(metricas, aes(Resolucion, Prop_Ambiguos)) +
    geom_col(fill = "#c0392b") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    tema_figura +
    labs(title = "Subclusters con Identidad Ambigua",
         subtitle = paste0("Marcadores excluyentes detectados en más del ",
                           scales::percent(umbral_deteccion_ambigua), " de las células"),
         x = "Resolución", y = "Porcentaje de Subclusters")
  
  p_panel <- p_arbol / (p_silueta | p_marcadores | p_ambiguos) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(title = paste("Diagnóstico de Resolución:", cfg$etiqueta),
                    subtitle = "Estabilidad, separación y capacidad discriminativa",
                    theme = tema_figura)
  
  print(p_panel)
  ggsave(file.path(out_dir, "graficos", paste0("Resolucion_", cfg$linaje, ".jpg")),
         plot = p_panel, device = "jpeg", width = 14, height = 12, dpi = 300)
  
  candidatas <- metricas$Resolucion[which(metricas$Prop_Ambiguos == 0)]
  if (length(candidatas) > 0) {
    r_sugerida <- min(candidatas)
    estado_r   <- "ambigüedad nula"
  } else if (any(!is.na(metricas$Prop_Ambiguos))) {
    r_sugerida <- metricas$Resolucion[which.min(metricas$Prop_Ambiguos)]
    estado_r   <- "mínimo alcanzado, ambigüedad no nula"
  } else {
    r_sugerida <- max(resoluciones)
    estado_r   <- "sin pares evaluables"
  }
  
  Idents(sub) <- paste0("RNA_snn_res.", r_sugerida)
  validos     <- cfg$panel[cfg$panel %in% rownames(sub)]
  
  suppressWarnings({
    p_dot <- DotPlot(sub, features = validos) +
      scale_color_viridis_c(option = "C") +
      tema_figura +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
      labs(title = cfg$etiqueta,
           subtitle = paste0("Resolución = ", format(r_sugerida, nsmall = 1, decimal.mark = "."),
                             " (", estado_r, ")"),
           x = "Gen", y = "Subcluster",
           colour = "Expresión Media", size = "Porcentaje Expresado")
    print(p_dot)
  })
  
  ggsave(file.path(out_dir, "graficos", paste0("PanelAmpliado_", cfg$linaje, ".jpg")),
         plot = p_dot, device = "jpeg", width = 9, height = 5, dpi = 300)
}

tabla_resumen <- bind_rows(resumen_global)
write.csv(tabla_resumen, file.path(out_dir, "datos", "Resolucion_Diagnostico.csv"), row.names = FALSE)
print(tabla_resumen)

cat("\n=================================================================\n")
cat("  Resolution diagnostic completed for", etapa, "\n")
cat("  A defensible resolution maximises silhouette width while keeping the\n")
cat("  weakest cluster above the marker count required by the signature.\n")
cat("=================================================================\n")