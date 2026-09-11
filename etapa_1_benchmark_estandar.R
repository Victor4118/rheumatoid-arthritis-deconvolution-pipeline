# ==============================================================================
# TFM | STAGE I - Benchmarking of deconvolution tools against FACS ground truth
# Author: Victor Daniel Velasco Perez
#
# Bulk RNA-seq and flow cytometry measured on the same validation samples.
# Each algorithm is applied to the bulk matrix and compared against the paired
# FACS fractions. Synovial tissue enters the work in Stage III, not here.
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(immunedeconv)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 0. SETUP & GLOBAL CONSTANTS
# ------------------------------------------------------------------------------
setwd("C:/Master en Bioinformatica/TFM")

data_dir <- "data"
out_dir  <- "resultados/etapa_1"
dir_graf <- file.path(out_dir, "graficos")
dir_tab  <- file.path(out_dir, "tablas")

invisible(lapply(
  c(dir_graf, dir_tab, file.path(out_dir, "datos")),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

CELLTYPES <- c("B cell", "T cell CD4+", "T cell CD8+", "NK cell")

# Output class is an a priori property of each algorithm (documented in literature).
CLASE <- c(
  CIBERSORTx    = "Fraccion", 
  EPIC          = "Fraccion", 
  quanTIseq     = "Fraccion",
  `MCP-counter` = "Puntuacion", 
  xCell         = "Puntuacion"
)

# Only these two accept a user-supplied custom signature (required for Stage III).
FIRMA_PROPIA <- c(
  CIBERSORTx    = TRUE, 
  EPIC          = TRUE, 
  quanTIseq     = FALSE,
  `MCP-counter` = FALSE, 
  xCell         = FALSE
)

PALETA <- c(
  CIBERSORTx    = "#0072B2", 
  EPIC          = "#D55E00", 
  quanTIseq     = "#009E73",
  `MCP-counter` = "#CC79A7", 
  xCell         = "#E69F00"
)

theme_tfm <- theme_bw(11) +
  theme(
    strip.background      = element_rect(fill = "grey92", colour = "grey30"),
    strip.text            = element_text(face = "bold"), 
    panel.grid.minor      = element_blank(),
    axis.title            = element_text(face = "bold"),
    plot.title            = element_text(face = "bold", size = 13),
    plot.subtitle         = element_text(colour = "grey25"),
    plot.caption          = element_text(hjust = 0, size = 8, colour = "grey35"),
    plot.title.position   = "plot", 
    plot.caption.position = "plot",
    legend.position       = "bottom"
  )

# ==============================================================================
# 1. INPUT AND QUALITY CONTROL OF THE EXPRESSION MATRIX
# ==============================================================================
load(file.path(data_dir, "facs_ground_truth.RData"))

facs <- as.data.frame(facs_ground_truth) %>%
  rename(Sample = sample, CellType = cell_type, FACS_Truth = true_fraction)

stopifnot(exists("bulk_expression_matrix"))
bulk <- as.matrix(bulk_expression_matrix)
mode(bulk) <- "numeric"

shared <- intersect(colnames(bulk), unique(facs$Sample))
stopifnot(length(shared) > 0)
N <- length(shared)

bulk_clean <- bulk[!is.na(rownames(bulk)) & rownames(bulk) != "", , drop = FALSE]
if (any(duplicated(rownames(bulk_clean)))) {
  bulk_clean <- rowsum(bulk_clean, rownames(bulk_clean))
}
bulk_clean[!is.finite(bulk_clean)] <- 0

qc_entrada <- tibble(
  comprobacion = c("Genes", "Muestras bulk", "Muestras pareadas", "Valores no finitos",
                   "Simbolos duplicados", "Valor maximo", "Escala inferida"),
  valor = c(
    nrow(bulk), 
    ncol(bulk), 
    N, 
    sum(!is.finite(bulk)),
    sum(duplicated(rownames(bulk))), 
    sprintf("%.1f", max(bulk_clean)),
    if (max(bulk_clean) < 50) "log2" else "lineal"
  )
)

print(as.data.frame(qc_entrada), row.names = FALSE)
write_csv(qc_entrada, file.path(dir_tab, "QC_entrada.csv"))

if (max(bulk_clean) < 50) {
  warning("Matrix appears log2-transformed: reverse the log before CIBERSORTx.")
}

abundancia <- facs %>% 
  filter(CellType %in% CELLTYPES, Sample %in% shared) %>%
  group_by(CellType) %>%
  summarise(
    min_FACS = min(FACS_Truth), 
    mediana  = median(FACS_Truth),
    max_FACS = max(FACS_Truth), 
    .groups  = "drop"
  ) %>% 
  arrange(desc(mediana))

write_csv(abundancia, file.path(dir_tab, "QC_abundancia_FACS.csv"))

# ==============================================================================
# 2. DECONVOLUTION
# ==============================================================================
IDS <- c(
  EPIC          = "epic", 
  quanTIseq     = "quantiseq",
  `MCP-counter` = "mcp_counter", 
  xCell         = "xcell"
)

RAW <- imap(IDS, function(id, nm) {
  tryCatch({
    r <- deconvolute(bulk, method = id, tumor = FALSE)
    colnames(r)[1] <- "cell_type"
    r
  }, error = function(e) {
    message(nm, " skipped: ", conditionMessage(e))
    NULL
  })
})
RAW <- Filter(Negate(is.null), RAW)

write.table(
  data.frame(GeneSymbol = rownames(bulk_clean), bulk_clean, check.names = FALSE),
  file.path(out_dir, "datos", "Bulk_Stage1_CIBERSORTx.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cbx_file <- file.path(data_dir, "CIBERSORTx_Results.csv")
if (file.exists(cbx_file)) {
  cbx <- read.csv(cbx_file, check.names = FALSE)
  colnames(cbx)[1] <- "Sample"
  
  diag_cols <- c("P-value", "Correlation", "RMSE", "Absolute score (sig.score)")
  
  cbx_long <- cbx %>% 
    select(-any_of(diag_cols)) %>%
    pivot_longer(-Sample, names_to = "cell_type", values_to = "Score") %>%
    mutate(Method = "CIBERSORTx")
  
  sum_cols <- function(cols) {
    cols <- intersect(cols, colnames(cbx))
    if (!length(cols)) rep(NA_real_, nrow(cbx))
    else rowSums(cbx[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  pred_cbx <- cbx %>% 
    transmute(
      Sample,
      `B cell`     = sum_cols(c("B cells naive", "B cells memory")),
      `T cell CD4+` = sum_cols(c("T cells CD4 naive", "T cells CD4 memory resting",
                                 "T cells CD4 memory activated",
                                 "T cells follicular helper", "T cells regulatory (Tregs)")),
      `T cell CD8+` = sum_cols("T cells CD8"),
      `NK cell`     = sum_cols(c("NK cells resting", "NK cells activated"))
    ) %>%
    pivot_longer(-Sample, names_to = "CellType", values_to = "Score") %>%
    mutate(Method = "CIBERSORTx")
} else {
  cbx_long <- tibble()
  pred_cbx <- tibble()
}

# ==============================================================================
# 3. OUTPUT SCALE
# ==============================================================================
completo <- bind_rows(
  imap_dfr(RAW, ~ pivot_longer(.x, -cell_type, names_to = "Sample", values_to = "Score") %>% 
             mutate(Method = .y)),
  cbx_long
)

qc_escala <- completo %>%
  group_by(Method, Sample) %>% 
  summarise(s = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(Method) %>% 
  summarise(suma_mediana = round(median(s), 3), .groups = "drop") %>%
  mutate(
    clase  = unname(CLASE[Method]),
    escala = case_when(
      clase == "Puntuacion"              ~ "puntuacion (u.a.)",
      between(suma_mediana, 0.80, 1.20)  ~ "fraccion, composicion cerrada",
      between(suma_mediana, 80, 120)     ~ "PORCENTAJE 0-100",
      suma_mediana > 0                   ~ "fraccion, composicion abierta",
      TRUE                               ~ "ANOMALA"
    )
  )

print(as.data.frame(qc_escala), row.names = FALSE)
write_csv(qc_escala, file.path(dir_tab, "QC_escala_salida.csv"))

pct <- qc_escala$Method[qc_escala$escala == "PORCENTAJE 0-100"]
if (any(qc_escala$escala == "ANOMALA")) {
  warning("Anomalous output scale: check for absolute mode before interpreting.")
}

# ==============================================================================
# 4. AGGREGATION TO THE FACS PANEL
# ==============================================================================
MAPA <- list(
  "B cell" = c("B cell", "B cell naive", "B cell memory", "B cell plasma",
               "B lineage", "Class-switched memory B cell"),
  "T cell CD4+" = c("T cell CD4+", "T cell CD4+ (non-regulatory)",
                    "T cell regulatory (Tregs)", "T cell CD4+ naive",
                    "T cell CD4+ memory", "T cell CD4+ central memory",
                    "T cell CD4+ effector memory", "T cell CD4+ Th1",
                    "T cell CD4+ Th2", "T cell follicular helper"),
  "T cell CD8+" = c("T cell CD8+", "T cell CD8+ naive", "T cell CD8+ central memory",
                    "T cell CD8+ effector memory", "Cytotoxic lymphocytes"),
  "NK cell" = c("NK cell", "NK cell resting", "NK cell activated")
)

agregar <- function(res, id, nm) {
  agg <- tryCatch(
    immunedeconv::map_result_to_celltypes(res, CELLTYPES, id) %>%
      rownames_to_column("CellType") %>%
      pivot_longer(-CellType, names_to = "Sample", values_to = "Score"),
    error = function(e) {
      largo <- pivot_longer(res, -cell_type, names_to = "Sample", values_to = "Score")
      imap_dfr(MAPA, function(sin, lin) {
        sub <- filter(largo, cell_type %in% sin)
        if (!nrow(sub)) {
          tibble(CellType = lin, Sample = unique(largo$Sample), Score = NA_real_)
        } else {
          sub %>% 
            group_by(Sample) %>%
            summarise(Score = sum(Score, na.rm = TRUE), .groups = "drop") %>%
            mutate(CellType = lin)
        }
      })
    }
  )
  mutate(agg, Method = nm)
}

trazabilidad <- imap_dfr(RAW, function(res, nm) {
  imap_dfr(MAPA, ~ tibble(
    Method        = nm, 
    CellType      = .y,
    tipos_nativos = paste(intersect(res$cell_type, .x), collapse = " + ")
  ))
})
write_csv(trazabilidad, file.path(dir_tab, "QC_agregacion_tipos.csv"))

predictions <- bind_rows(imap_dfr(RAW, ~ agregar(.x, IDS[[.y]], .y)), pred_cbx) %>%
  filter(CellType %in% CELLTYPES, Sample %in% shared) %>%
  mutate(
        Score       = if_else(Method %in% pct, Score / 100, Score),
    Tipo_Salida = unname(CLASE[Method])
  ) %>%
  inner_join(facs, by = c("Sample", "CellType"))

# ==============================================================================
# 5. DENOMINATOR AGREEMENT
# ==============================================================================
suma_facs <- facs %>% 
  filter(CellType %in% CELLTYPES, Sample %in% shared) %>%
  group_by(Sample) %>% 
  summarise(suma_FACS = sum(FACS_Truth), .groups = "drop")

qc_denom <- predictions %>% 
  filter(Tipo_Salida == "Fraccion") %>%
  group_by(Method, Sample) %>% 
  summarise(suma_metodo = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  left_join(suma_facs, by = "Sample") %>% 
  mutate(cociente = suma_metodo / suma_FACS)

write_csv(qc_denom, file.path(dir_tab, "QC_denominadores.csv"))

print(
  qc_denom %>% 
    group_by(Method) %>%
    summarise(cociente_mediano = round(median(cociente), 2), .groups = "drop") %>%
    as.data.frame(), 
  row.names = FALSE
)

# ==============================================================================
# 6. AVAILABILITY AND METRICS
# ==============================================================================
disponibilidad <- expand_grid(Method = names(CLASE), CellType = CELLTYPES) %>%
  filter(Method %in% unique(predictions$Method)) %>%
  left_join(
    predictions %>% 
      group_by(Method, CellType) %>%
      summarise(
        n_obs = sum(!is.na(Score)), 
        sd    = sd(Score, na.rm = TRUE),
        media = mean(Score, na.rm = TRUE), 
        .groups = "drop"
      ),
    by = c("Method", "CellType")
  ) %>%
  mutate(estado = case_when(
    is.na(n_obs) | n_obs == 0 ~ "no definido por el metodo",
    media == 0                ~ "estimacion nula",
    is.na(sd) | sd == 0       ~ "constante entre muestras",
    TRUE                      ~ "evaluable"
  ))

write_csv(disponibilidad, file.path(dir_tab, "QC_disponibilidad.csv"))

P_MIN <- 2 / factorial(N)
perms <- function(x) {
  if (length(x) <= 1) list(x) 
  else do.call(c, lapply(seq_along(x), function(i) lapply(perms(x[-i]), function(p) c(x[i], p))))
}

p_exacto <- function(a, b) {
  if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  if (length(a) > 8) return(suppressWarnings(cor.test(a, b, method = "spearman")$p.value))
  obs <- cor(a, b, method = "spearman")
  mean(abs(vapply(perms(seq_along(a)), function(p) cor(a, b[p], method = "spearman"), numeric(1))) >= abs(obs) - 1e-9)
}
message(sprintf("n = %d | minimum attainable two-sided p = %.3f", N, P_MIN))

metrics <- predictions %>%
  group_by(Method, CellType, Tipo_Salida) %>%
  summarise(
    n = sum(!is.na(Score)), 
    v = sd(Score, na.rm = TRUE),
    Spearman = if (n >= 3 && !is.na(v) && v > 0) round(cor(FACS_Truth, Score, method = "spearman"), 3) else NA_real_,
    p_exacto = if (n >= 3 && !is.na(v) && v > 0) round(p_exacto(FACS_Truth, Score), 3) else NA_real_,
    Pearson  = if (n >= 3 && !is.na(v) && v > 0) round(cor(FACS_Truth, Score), 3) else NA_real_,
    RMSE     = if (first(Tipo_Salida) == "Fraccion") round(sqrt(mean((Score - FACS_Truth)^2)), 4) else NA_real_,
    NRMSE    = if (first(Tipo_Salida) == "Fraccion") round(sqrt(mean((Score - FACS_Truth)^2)) / diff(range(FACS_Truth)), 3) else NA_real_,
    Sesgo    = if (first(Tipo_Salida) == "Fraccion") round(mean(Score - FACS_Truth), 4) else NA_real_,
    .groups  = "drop"
  ) %>%
  mutate(supera_rango = if_else(!is.na(NRMSE) & NRMSE > 1, "SI", "no")) %>%
  select(-v) %>% 
  arrange(Tipo_Salida, Method, CellType)

print(as.data.frame(metrics), row.names = FALSE)
write_csv(metrics, file.path(dir_tab, "T_metricas_benchmark.csv"))

elegibilidad <- tibble(
  Method       = names(CLASE), 
  clase_salida = unname(CLASE),
  firma_propia = unname(FIRMA_PROPIA)
) %>%
  filter(Method %in% unique(predictions$Method)) %>%
  mutate(decision = if_else(firma_propia, "PASA A ETAPA III", "excluido"))

write_csv(elegibilidad, file.path(dir_tab, "T_elegibilidad_etapaIII.csv"))

# ==============================================================================
# 7. FIGURES
# ==============================================================================
pie <- sprintf(
  paste(
    "n = %d muestras pareadas. Metricas descriptivas: el valor p",
    "bilateral minimo alcanzable es %.3f, por lo que ningun\ncontraste",
    "puede alcanzar significacion estadistica."
  ), 
  N, P_MIN
)

df_frac <- filter(predictions, Tipo_Salida == "Fraccion")
lims <- df_frac %>% 
  group_by(CellType) %>%
  summarise(hi = max(c(FACS_Truth, Score)) * 1.08, .groups = "drop")

blank <- bind_rows(
  transmute(lims, CellType, FACS_Truth = 0, Score = 0),
  transmute(lims, CellType, FACS_Truth = hi, Score = hi)
)

p_exactitud <- ggplot(df_frac, aes(FACS_Truth, Score, colour = Method)) +
  geom_blank(data = blank, inherit.aes = FALSE, aes(FACS_Truth, Score)) +
  geom_abline(linetype = "dashed", colour = "grey35") +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_wrap(~ CellType, scales = "free") +
  scale_colour_manual(values = PALETA, name = "Metodo:") +
  labs(
    title    = "Exactitud de los metodos que estiman fraccion celular",
    subtitle = "Estimacion frente a la fraccion medida por citometria de flujo",
    x        = "Fraccion medida por citometria", 
    y        = "Fraccion estimada",
    caption  = paste(
      "La linea discontinua es la identidad. Los ejes comparten rango",
      "dentro de cada panel, de modo que la distancia vertical a la",
      "recta es\nel error. Por encima: sobreestimacion; por debajo:",
      "subestimacion.\n", pie
    )
  ) +
  theme_tfm + 
  theme(aspect.ratio = 1)

heat <- disponibilidad %>% 
  select(Method, CellType, estado) %>%
  left_join(select(metrics, Method, CellType, Spearman), by = c("Method", "CellType")) %>%
  mutate(
    lab = case_when(
      !is.na(Spearman)                  ~ sprintf("%.2f", Spearman),
      estado == "no definido por el metodo" ~ "n/d",
      estado == "estimacion nula"       ~ "0", 
      TRUE                              ~ "cte"
    ),
    Method = factor(Method, rev(names(CLASE)))
  )

p_heatmap <- ggplot(heat, aes(CellType, Method, fill = Spearman)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = lab), size = 3.9) +
  scale_fill_gradient2(
    limits   = c(-1, 1), 
    low      = "#0072B2", 
    mid      = "grey97",
    high     = "#D55E00", 
    na.value = "grey90", 
    name     = "Rho de\nSpearman"
  ) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title    = "Concordancia de ordenamiento con la citometria de flujo",
    subtitle = "Correlacion de Spearman entre la estimacion y la fraccion medida",
    x        = "Linaje celular", 
    y        = NULL,
    caption  = paste(
      "Metrica aplicable a ambas clases de salida: mide asociacion",
      "monotona y no exige unidades homogeneas. n/d indica que el\n",
      "metodo no define ese linaje, una ausencia estructural.\n", pie
    )
  ) +
  theme_tfm + 
  theme(
    panel.grid     = element_blank(), 
    axis.ticks     = element_blank(),
    legend.position = "right"
  )

p_denom <- ggplot(qc_denom, aes(Sample, suma_metodo, colour = Method)) +
  geom_point(
    data         = suma_facs, 
    inherit.aes  = FALSE, 
    aes(Sample, suma_FACS),
    shape        = 4, 
    size         = 3.4, 
    stroke       = 1.1, 
    colour       = "grey15"
  ) +
  geom_point(size = 3, alpha = 0.9, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = PALETA, name = "Metodo:") +
  ylim(0, NA) +
  labs(
    title    = "Cobertura del panel linfoide y comparabilidad de denominadores",
    subtitle = "Suma de los cuatro linajes en cada muestra, por metodo y por citometria",
    x        = "Muestra", 
    y        = "Suma de los cuatro linajes",
    caption  = paste(
      "Las aspas corresponden a la citometria. Un metodo situado",
      "sistematicamente por encima normaliza sus fracciones frente a",
      "una\npoblacion mas estrecha que la acotada por el citometro:",
      "parte del sesgo positivo procede de esa diferencia de",
      "denominador y no\nde una perdida de exactitud.\n", pie
    )
  ) +
  theme_tfm

p_punt <- predictions %>% 
  filter(Tipo_Salida == "Puntuacion") %>%
  group_by(Method, CellType) %>% 
  filter(sd(Score, na.rm = TRUE) > 0) %>% 
  ungroup() %>%
  ggplot(aes(FACS_Truth, Score, colour = Method)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_point(size = 2.6, alpha = 0.9) +
  facet_grid(Method ~ CellType, scales = "free") +
  scale_colour_manual(values = PALETA, guide = "none") +
  labs(
    title    = "Metodos basados en puntuacion frente a la citometria de flujo",
    subtitle = "La puntuacion es comparable entre muestras, nunca entre tipos celulares",
    x        = "Fraccion medida por citometria", 
    y        = "Puntuacion en unidades arbitrarias",
    caption  = paste(
      "Sin linea de identidad: una puntuacion en unidades arbitrarias",
      "no es dimensionalmente comparable con una fraccion, por lo que",
      "el\nerror cuadratico medio y el sesgo no estan definidos.\n", pie
    )
  ) +
  theme_tfm

p_nrmse <- metrics %>% 
  filter(Tipo_Salida == "Fraccion", !is.na(NRMSE)) %>%
  ggplot(aes(CellType, NRMSE, fill = Method)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey20") +
  geom_col(
    position = position_dodge2(preserve = "single"), 
    width    = 0.72,
    colour   = "grey20", 
    linewidth = 0.3
  ) +
  scale_fill_manual(values = PALETA, name = "Metodo:") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Error normalizado y suelo de deteccion por linaje",
    subtitle = "Error cuadratico medio dividido por el rango dinamico de la fraccion medida",
    x        = "Linaje celular", 
    y        = "Error normalizado",
    caption  = paste(
      "La normalizacion permite comparar linajes cuya abundancia",
      "difiere en dos ordenes de magnitud. Por encima de la linea de",
      "puntos\nel error supera el rango completo de variacion medido:",
      "el metodo no discrimina entre muestras para ese linaje.\n", pie
    )
  ) +
  theme_tfm

# Save figures loop
figure_list <- list(
  list(p_exactitud, "exactitud_fracciones", 9, 8),
  list(p_heatmap,   "concordancia_spearman", 8.5, 5),
  list(p_denom,     "cobertura_denominadores", 8.5, 5.5),
  list(p_punt,      "metodos_puntuacion", 10, 7),
  list(p_nrmse,     "suelo_deteccion", 8.5, 5.5)
)

for (f in figure_list) {
  suppressWarnings(
    ggsave(
      filename = file.path(dir_graf, paste0(f[[2]], ".jpg")), 
      plot     = f[[1]],
      width    = f[[3]], 
      height   = f[[4]], 
      dpi      = 300, 
      bg       = "white"
    )
  )
}

writeLines(capture.output(sessionInfo()), file.path(dir_tab, "sessionInfo.txt"))
message("Stage I complete: ", out_dir)
