---
title: "PEC1 ANÁLISIS DE DATO ÓMICOS - JULIA ESCUDERO FELIU"
output:
  word_document: default
  html_document: default
  pdf_document: default
---


knitr::opts_chunk$set(echo = TRUE)

options(repos = c(CRAN = "https://cloud.r-project.org/"))


## Selección y carga del dataset 

# Cargamos los datos descargados en en R
data <- read.table("/Users/iei/NASH-ALL-Count.txt", header = TRUE, sep = "\t", row.names = 1)

# Verificamos los datos cargados
head(data)

```
## Carga de metadatos y creación del objeto SummarizedExperiment

# Cargamos los metadatos de las muestras
sample_metadata <- read.table("/Users/iei/RNAseq_codes.txt", header = TRUE, sep = "\t")

# Verificamos, igual que antes, los datos cargados
head(sample_metadata)

# Creación del objeto: 

# Cargamos el paquete SummarizedExperiment
library(SummarizedExperiment)

# Configuramos los nombres de las filas en los metadatos
rownames(sample_metadata) <- sample_metadata$Sample_ID

# Filtramos columnas de 'data' y filas de 'sample_metadata' para eliminar elementos "NA"
data <- data[, !grepl("^NA", colnames(data))]
sample_metadata <- sample_metadata[!grepl("^NA", rownames(sample_metadata)), ]

# Filtramos 'data' y 'sample_metadata' para que solo contengan las muestras en común
data <- data[, colnames(data) %in% rownames(sample_metadata)]
sample_metadata <- sample_metadata[rownames(sample_metadata) %in% colnames(data), ]

# Antes de crear el objeto, nos aseguramos de que las muestras estén en el mismo orden
sample_metadata <- sample_metadata[match(colnames(data), rownames(sample_metadata)), ]

# Creamos el objeto SummarizedExperiment
se <- SummarizedExperiment(
    assays = list(counts = as.matrix(data)),
    colData = sample_metadata
)

# Verificamos el objeto creado
print(se)
summary(se)

## Análisis exploratorio del dataset

# Resumen estadístico básico de los datos de conteo
counts <- assay(se)  # Extraer la matriz de conteos

# Calculamos estadísticas descriptivas por transcrito (filas)
row_stats <- data.frame(
    mean = rowMeans(counts),
    median = apply(counts, 1, median),
    sd = apply(counts, 1, sd),
    min = apply(counts, 1, min),
    max = apply(counts, 1, max)
)

# Mostramos las primeras filas del resumen estadístico
head(row_stats)

# Histograma de la media de expresión por transcrito
hist(row_stats$mean, breaks = 50, main = "Distribución de la media de expresión por transcrito",
     xlab = "Media de expresión", ylab = "Frecuencia")

# Boxplot de los valores de expresión por muestra
boxplot(counts, outline = FALSE, main = "Boxplot de expresión por muestra",
        ylab = "Nivel de expresión", xlab = "Muestras")

# Calculamos la matriz de correlación entre muestras
cor_matrix <- cor(counts)

# Visualizamos la matriz de correlación (se requiere el paquete pheatmap para la visualización así que primero lo instalamos)
if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
}
library(pheatmap)
pheatmap(cor_matrix, main = "Matriz de correlación entre muestras")


save(se, file = "SummarizedExperiment.Rda")



