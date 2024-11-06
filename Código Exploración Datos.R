options(repos = c(CRAN = "https://cran.rstudio.com/"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
install.packages("dplyr")

library(SummarizedExperiment)
library(dplyr)

# Establecemos nuestro directorio de trabajo
setwd("~/Omiques/metaboData/Datasets/2018-MetabotypingPaper")


# Cargamos los archivos
data_info <- read.csv("DataInfo_S013.csv", stringsAsFactors = FALSE)
data_values <- read.csv("DataValues_S013.csv", stringsAsFactors = FALSE)
descripcion <- readLines("description.md")

# Arreglamos la estructura
data_values <- data_values %>% select(-X.1)
rownames(data_values) <- data_values$SUBJECTS
data_values <- data_values %>% select(-SUBJECTS)
data_matrix <- as.matrix(data_values)
col_data <- data_info %>% select(VarName, varTpe, Description)

# Aseguramos que las columnas de data_matrix coincidan con los nombres en col_data
common_vars <- intersect(colnames(data_matrix), col_data$VarName)

# Filtramos data_matrix y col_data para mantener solo las variables comunes
data_matrix <- data_matrix[, common_vars, drop = FALSE]
col_data <- col_data[col_data$VarName %in% common_vars, ]

# Crear el objeto SummarizedExperiment
sumexp <- SummarizedExperiment(
  assays = SimpleList(counts = data_matrix),
  colData = col_data
)


# Comparamos la primera columna de data_matrix con la primera fila de col_data para comprobar que nuestros datos cuadran.

first_match <- colnames(data_matrix)[1] == col_data$VarName[1]
last_match <- colnames(data_matrix)[ncol(data_matrix)] == col_data$VarName[nrow(col_data)]


# Agregamos nuestro archivo description como nuestra metadata
metadata(sumexp)$description <- paste(descripcion, collapse = "\n")

# Revisar el objeto
print(sumexp)

# Dimensiones del objeto SummarizedExperiment
dim(sumexp)

# Nombres de las filas y columnas para una visión general
head(rownames(sumexp))
head(colnames(sumexp))

#Exploramos diferentes componentes de nuestro Summarized Experiment
colData(sumexp)
rowData(sumexp)
assayNames(sumexp)
summary(colData(sumexp))

# Convertir la matriz assay(sumexp) a un data frame temporal
assay_df <- as.data.frame(assay(sumexp))

# Intentar convertir cada columna de assay_df a numérica
assay_df_numeric <- as.data.frame(lapply(assay_df, function(x) as.numeric(as.character(x))))

# Filtrar solo las columnas numéricas
numeric_data <- assay_df_numeric[, sapply(assay_df_numeric, function(x) !all(is.na(x)))]

# Calcular la media de las columnas numéricas, ignorando NA
meancol <- apply(numeric_data, 2, mean, na.rm = TRUE)

# Mostrar los primeros resultados
head(meancol)

# Filtramos las columnas con varianza mayor que cero
numeric_data_nonconstant <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) != 0]

# Escalamos los datos y reemplazamos los valores NA por la mediana de esa columna 
numeric_data_scaled <- scale(numeric_data_nonconstant, center = TRUE, scale = TRUE)
numeric_data_scaled[is.na(numeric_data_scaled)] <- apply(numeric_data_scaled, 2, function(x) median(x, na.rm = TRUE))

# Realizamos el PCA con los datos ya limpios y escalados
pca_result <- prcomp(numeric_data_scaled, center = TRUE, scale. = TRUE)

# Extraemos las tres primeras componentes
pca_scores <- pca_result$x[, 1:3]  # Accede a las tres primeras columnas de las puntuaciones

head(pca_scores)

library(plotly)

# Creamos gráfico 3D para las 3 primeras componentes principales
plot_ly(x = pca_scores[, 1], y = pca_scores[, 2], z = pca_scores[, 3], 
        type = "scatter3d", mode = "markers") %>%
  layout(title = "3D PCA",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))
