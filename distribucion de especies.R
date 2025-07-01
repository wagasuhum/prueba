library(SpaDES)
library(raster)
library(xfun)
library(ReIns)
library(dplyr)


# load("rep_distspp_cod/rep_distspp_objetos_operativo.RData")

# 1. Cargar insumos

# 1.1 Cargar rutas de los modelos Nivel1 BioModelos

list_raster <-list.files(path = "C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/BioModelos_N1", pattern="*.tif$", full.names = T)

# 1.2 Leer archivos RUNAP, en formato raster (.tif).

RUNAP_1990 <- raster("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_rasters/RUNAP_1990_pnnid.tif")
RUNAP_2010 <- raster("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_rasters/RUNAP_2010_pnnid.tif")


# 2. Funciones

# 2.1 Representatividad de especies por periodo de SINAP
#
# raster_spec: raster stack, de modelos de distribución de especies  
# SINAP_shp: shape, de periodo de SINAP
# name_table: vector character, para darle nombre a la matriz de riqueza, describe el contenido y temporalidad del SINAP
#
# return matriz de riqueza de especies en el conjunto de Areas Protegidas (AP) del SINAP
# por especie y vector de representatividad total 

rep_spp_SINAP <- function(raster_spp, SINAP_shp, name_table){
  # De Raster Stack binario a data.frame binaria
  x <- as.data.frame(raster_spp, xy=TRUE)
  
  # Espacializar el data frame
  coordinates(x) <- ~ x + y
  
  # Homologar sistemas de coordenadas de insumos
  SINAP_shp <- SINAP_shp  %>% spTransform(crs(raster_spp))
  crs(x) <- crs(raster_spp)
  
  # Cruzar datos geográficos de los insumos 
  temp <- over(SINAP_shp, x)
  
  #generar matriz presencia-ausencia, numero de Areas Protegidas por especie
  temp[is.na(temp)] = 0
  temp <- t(temp)
  
  # Cantidad de Areas Protegidas en las cuales esta presente cada especie
  ap_spec.sum <-apply(temp,1,sum)
  ap.riqueza.1 <-cbind(temp,ap_spec.sum)
  
  # presencia-ausencia de especies en el conjunto de AP dentro del SINAP
  nc <- ncol(ap.riqueza.1) # columna en la que se ubica el conteo
  pres_espec_all_ap <- replace(ap.riqueza.1[,nc], ap.riqueza.1[,nc] > 0 , 1)
  ap.riqueza.1 <-cbind(ap.riqueza.1,pres_espec_all_ap)
  
  # Porcentaje de la representatividad de especies en las AP, en otras palabras,
  # porcion de las especies usadas en el calculo que estan presentes en el 
  # conjunto del AP's
  sumPresAus <- sum(ap.riqueza.1[,ncol(ap.riqueza.1)])
  repre_riqueza <-sumPresAus/(dim(espec_mod)[[3]])*100
  list_data <- list(ap.riqueza.1,repre_riqueza)
  names(list_data) <- c(paste0("Matriz riqueza"," ", name_table), 
                        "Representatividad total")
  return(list_data)
}
