# Diagrama actualizado del flujo completo del TFM.
# Genera versiones PNG (para Word) y PDF (vectorial).

suppressPackageStartupMessages(library(grid))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

project_root <- dirname(script_dir)
out_dir <- file.path(project_root, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

colours <- list(
  config = "#C93A2C",
  analysis = "#7D3C98",
  biology = "#216A8C",
  export = "#238B57",
  final = "#087F73",
  arrow = "#4D4D4D",
  text = "#FFFFFF",
  side = "#424242"
)

nodes <- data.frame(
  y = c(0.925, 0.830, 0.735, 0.630, 0.525,
        0.420, 0.310, 0.200, 0.105, 0.035),
  h = c(0.052, 0.064, 0.064, 0.080, 0.068,
        0.068, 0.082, 0.068, 0.060, 0.050),
  label = c(
    "M0 - Configuraci\u00f3n global\nk = 5 - semilla - modo all",
    "M1 - Preprocesamiento\nFiltro - varianza - VST - selecci\u00f3n combinada - z-score",
    "M2 - Reducci\u00f3n dimensional\nPCA - mapas de difusi\u00f3n",
    "M3 - Clustering\nK-means - jer\u00e1rquico - espectral - consenso\nLeiden - MCL - NMF - GMM",
    "M4 - Validaci\u00f3n y selecci\u00f3n\nSilueta - estabilidad por remuestreo - comparaci\u00f3n histol\u00f3gica posterior",
    "M5 - Expresi\u00f3n diferencial\nDESeq2 one-versus-rest - FDR - log2 fold-change",
    "M6 - Caracterizaci\u00f3n biol\u00f3gica posterior\nComposici\u00f3n histol\u00f3gica - firmas - genes marcadores\nEnriquecimiento Reactome - prueba hipergeom\u00e9trica - correcci\u00f3n BH",
    "M7 - Visualizaci\u00f3n\nPCA - heatmap final - volcanes - matriz porcentual - Sankey",
    "M8 - Exportaci\u00f3n\nFiguras - tablas - asignaciones - par\u00e1metros",
    "Clusters transcript\u00f3micos"
  ),
  fill = c(
    colours$config,
    rep(colours$analysis, 4),
    rep(colours$biology, 3),
    colours$export,
    colours$final
  ),
  stringsAsFactors = FALSE
)

draw_workflow <- function() {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = "sans")))

  grid.text(
    "Pipeline de an\u00e1lisis transcript\u00f3mico",
    x = 0.46, y = 0.985,
    gp = gpar(fontsize = 17, fontface = "bold", col = "#222222")
  )

  # Flechas primero para que queden por debajo de las cajas.
  for (i in seq_len(nrow(nodes) - 1)) {
    y_start <- nodes$y[i] - nodes$h[i] / 2
    y_end <- nodes$y[i + 1] + nodes$h[i + 1] / 2
    grid.lines(
      x = unit(c(0.46, 0.46), "npc"),
      y = unit(c(y_start, y_end), "npc"),
      arrow = arrow(type = "closed", length = unit(2.2, "mm")),
      gp = gpar(col = colours$arrow, lwd = 1.4)
    )
  }

  for (i in seq_len(nrow(nodes))) {
    is_final <- i == nrow(nodes)
    width <- if (is_final) 0.43 else 0.60
    radius <- if (is_final) unit(0.025, "npc") else unit(0.012, "npc")

    grid.roundrect(
      x = 0.46, y = nodes$y[i],
      width = width, height = nodes$h[i],
      r = radius,
      gp = gpar(fill = nodes$fill[i], col = NA)
    )
    grid.text(
      nodes$label[i], x = 0.46, y = nodes$y[i],
      gp = gpar(
        col = colours$text,
        fontsize = if (i == 4) 9.0 else if (i %in% c(2, 5, 7, 8)) 9.2 else 9.7,
        fontface = "bold",
        lineheight = 1.05
      )
    )
  }

  # Etiquetas laterales que resumen las fases del trabajo.
  grid.text("PROCESAMIENTO\nY AN\u00c1LISIS", x = 0.84, y = 0.72,
            gp = gpar(col = colours$analysis, fontsize = 10, fontface = "bold"))
  grid.text("INTERPRETACI\u00d3N\nPOSTERIOR", x = 0.84, y = 0.31,
            gp = gpar(col = colours$biology, fontsize = 10, fontface = "bold"))
  grid.text("COMUNICACI\u00d3N\nY EXPORTACI\u00d3N", x = 0.84, y = 0.12,
            gp = gpar(col = colours$export, fontsize = 10, fontface = "bold"))

  popViewport()
}

png(
  filename = file.path(out_dir, "pipeline_workflow_updated.png"),
  width = 2400, height = 3300, res = 300, bg = "white"
)
draw_workflow()
dev.off()

pdf(
  file = file.path(out_dir, "pipeline_workflow_updated.pdf"),
  width = 8, height = 11, family = "Helvetica", useDingbats = FALSE
)
draw_workflow()
dev.off()

cat("Workflow created in:", out_dir, "\n")
