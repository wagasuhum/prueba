# Bien Representado Ecológicamente
#
# Indicador: Porcentaje y cambio en el porcentaje de la representatividad de la riqueza de especies del SINAP
#
# Porcentaje y cambio en el porcentaje de la representatividad de la riqueza de especies corresponde al porcentaje 
# y cambio en el porcentaje de especies protegidas dentro de las unidades del patrimonio natural y cultural del país 
# incluidas dentro del Sistema Nacional de Áreas Protegidas con relación al número de especies disponibles para el 
# análisis

# librerias necesaria
library("dismo")
library("sf")
library("rgdal")
library("raster")
library("qpcR")
library("dplyr")

# Configurando carpetas

dir.create("productos", showWarnings = F)
dir.create("productos/rep_ri", showWarnings = F)
dir.create("productos/rep_ri/", showWarnings = F)
