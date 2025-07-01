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
library("raster")
library("qpcR")
library("dplyr")

# Configurando carpetas

dir.create("productos", showWarnings = F)
dir.create("productos/rep_ri", showWarnings = F)
dir.create("productos/rep_ri/", showWarnings = F)


# 1.1 modelos de especies
espec_mod <- stack(list.files("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/BioModelos_N1","tif",full.names=TRUE)) 

# 1.2 Cargar SINAP alias historico RUNAP a WGS84.

ap1990 <- shapefile("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_shapefiles/RUNAP_1990.shp")
ap2010 <- shapefile("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_shapefiles/RUNAP_2010.shp")

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


# 2.2 Calculo del aporte de cada Area Protegida (AP) a la representatividad 
#de la riqueza de especies.
# matriz_rique: matriz, de riqueza de especies calculada con rep_spp_SINAP()  
# tiempo: vector character, de temporalidad de las unidades de analisis, por 
#ejemplo el año el que representa las AP's 
#
# return data.frame, porcion de representatividad por cada AP

repre_ind_AP <- function(matriz_rique = matriz_rique, tiempo){
  
  # extraer matriz binaria calculada con la funcion rep_spp_SINAP
  mat_bin <-matriz_rique[[1]]
  temp <- mat_bin[,1:(ncol(mat_bin)-2)]/mat_bin[,(ncol(mat_bin)-1)]
  temp[is.na(temp)] <- 0
  temp <-apply(temp,2,sum)
  
  # Porcion de representatividad de cada especie en cada UAE
  repre_sp_AP <- 100*(as.data.frame(temp))/nrow(mat_bin)
  names(repre_sp_AP) <- paste0("RRi", tiempo)
  return(repre_sp_AP)
}


# 2.3 Cambio en la representatividad de especies en SINAP
#
# cambio_rep_anual: vector numeric, de la representatividad de especies total 
# del conjunto de Areas Protegidas por periodo

cambio_rep_anual <- function(rep_aps = Rep_aps){
  cambio.rep.anual<- abs(diff(rep_aps))
  colnames(cambio.rep.anual) <- "cambio.rep.anual"
  return(cambio.rep.anual)
}

# 2.4 Calculo del aporte de cada Area Protegida (AP) al cambio de la representatividad 
#de la riqueza de especies.
# Establece el cambio en la representatividad de la riqueza de especies entre 
# un par de años consecutivos (t2 y t1). En terminos simples, calcula t2 - t1 por
# area protegida.
# shp1: shapefile, datos vectoriales del SINAP para el periodo 1
# rique_ind1: data.frame, aporte de riqueza calculada con la funcion rep_ind_AP
# en el periodo 1
# tiempo1: vector character, de temporalidad del SINAP en el periodo 1. Por 
#ejemplo el año el que representa las AP's 
# rique_ind2: matriz, de riqueza de especies calculada con rep_spp_SINAP() 
# para el periodo 2 
# tiempo1: vector character, de temporalidad del SINAP en el periodo 2. Por 
#ejemplo el año el que representa las AP's 
#
# return data.frame, porcion de aporte en el cambio de la representatividad por cada AP

repre_deltaind_AP <- function(shp1, rique_ind1, tiempo1, shp2, rique_ind2, tiempo2 ){
  
  # Areas protegidas que comparten los periodos
  index <-  shp2@data$IDPNN %in% shp1@data$IDPNN
  
  # Restar el valor del aporte de la riqueza por AP en los dos periodos
  diff_21 <-  rique_ind2[which(index == T) , 1] - rique_ind1[, 1]
  
  # crear un vector en donde se guarde la diferencia de los dos vectores
  deltaRRi  <- rep(NA, nrow(shp2@data))
  deltaRRi[which(index == T)] <- diff_21
  
  # Porcion del cambio de representatividad de cada especie por AP
  repre_delta_AP <- as.data.frame(deltaRRi) 
  names(repre_delta_AP) <- paste0("dRRi_", tiempo1, "_", tiempo2)
  return(repre_delta_AP)
}


# 3. Aplicación 

# 3.1 Indicador
# Ic: Porcentaje de la representatividad de la riqueza de especies del SINAP (%RRi)

# 3.1.1 Aplicar la funcion rep_spp_SINAP() para toda la serie de shapefiles 
# (historico RUNAP)

repres_90 <- rep_spp_SINAP(espec_mod, ap1990,"RUNAP1990")
repres_10 <- rep_spp_SINAP(espec_mod, ap2010,"RUNAP2010")

# extraer representatividad por cada temporalidad
Rep_aps <- rbind(
  repres_90$`Representatividad total`,
  repres_10$`Representatividad total`
)

colnames(Rep_aps) <- "Repre_APs"
rownames(Rep_aps) <- c("RUNAP 1990",
                       "RUNAP 2010"
)

# 3.2 Indicador
# Id: Cambio en el porcentaje de la representatividad de la riqueza de especies
# del SINAP (dRRi)
# Incremento tiempo 1 a tiempo 2 en la representatividad riqueza total

tasa.increm.anual <- cambio_rep_anual(rep_aps = Rep_aps)

# 3.3 Aporte de cada Area Protegida

# 3.3.1 Aporte de cada AP a la representatividad de la riqueza de especies 
# (RRi). Se evalua el % de representatividad de cada AP al total

RP_AP1990_AP <- repre_ind_AP(repres_90, "1990")
RP_AP2010_AP <- repre_ind_AP(repres_10, "2010")

# 3.3.2 Aporte de cada AP al cambio de la representatividad de la riqueza de 
# especies (dRRi). Se evalua el delta del % de representatividad de cada AP
# total

deltaRP_AP1990_AP <- repre_deltaind_AP(ap1990, RP_AP1990_AP, "90", ap1990, RP_AP1990_AP, "90")
deltaRP_AP2010_AP <- repre_deltaind_AP(ap1990, RP_AP1990_AP, "90", ap2010, RP_AP2010_AP, "10")

# 3.3.3 Agregar columnas %RRi y deltaRRi al shapefile original

ap1990@data <- cbind(ap1990@data, round(RP_AP1990_AP, 3), round(deltaRP_AP1990_AP, 3))
ap2010@data <- cbind(ap2010@data, round(RP_AP2010_AP, 3), round(deltaRP_AP2010_AP, 3))

# 4. Se escribe un shapefile nuevo que incluya la representatividad y el cambio en
# la representatividad contenidas en las columnas %RRi(año) y deltaRRi.



# Convertir de sp a sf
ap1990_sf <- st_as_sf(ap1990)

# Reproyectar
ap1990_proj <- st_transform(ap1990_sf, crs = 4326)

st_write(ap1990_proj, "productos/rep_ri/RUNAPRRi/RUNAP_1990_rri.shp", delete_layer = TRUE)


# 5 Representatividad riqueza por territorial

# 5.1 Preparación y calculo de representatividad por territorial

# 5.1.1. Se aplica la herramienta "intersect" de mapas en ArcGis 10.5 entre shapefile territoriales con multitemporal
# RUNAP. El shapefile de direcciones territoriales se encuentra en capas_base_ejemplos/Territorial.
# Una vez realizada la interseccion guardar el shp dentro de la carpeta del proyecto.

# 5.1.2. Cargar todos los dbf de RUNAP-territorialb
mapas_runap_ter <- list.files("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_rasters", pattern="terr.dbf$", full.names=TRUE)
mapas_runap_ter <- lapply(mapas_runap_ter, foreign::read.dbf)

# 5.1.3. Se agregan los porcentajes por Territorial y se hace una sumatoria.

terr_2010 <- aggregate(RRi2010 ~ nombre_1, data = mapas_runap_ter[[1]], sum)
colnames(terr_2010) <- c("Terr", "RRi2010")


terr_1990 <- aggregate(RRi1990 ~ nombre_1, data = mapas_runap_ter[[1]], sum)
colnames(terr_1990) <- c("Terr", "RRi1990")


# 5.1.4. Se unen todos los valores de cada RUNAP en una sola tabla.

dat_rep_ri_terr <- as.data.frame(t(cbind(terr_1990[,2], terr_2010[,2])))

# 5.2 Calculo de la diferencia en la representatividad entre territoriales, dRRIterr 

tasa.incremento.territorial <- apply(dat_rep_ri_terr, 2, diff) %>% as.data.frame() # al trabajar varios años, se hace necesario transponer la matriz inicial
row.names(tasa.incremento.territorial) <- as.vector(terr_1990[,1])
colnames(tasa.incremento.territorial) <- "1990-2010"

# save(espec_mod, ap1990, ap2010, rep_spp_SINAP,
#      repre_ind_AP, cambio_rep_anual, repre_deltaind_AP, Rep_aps, tasa.increm.anual, RP_AP1990_AP, 
#      RP_AP2010_AP, deltaRP_AP1990_AP, deltaRP_AP2010_AP, mapas_runap_ter, terr_1990, 
#      terr_2010, dat_rep_ri_terr, tasa.incremento.territorial, file = "productos/rep_ri/rep_ri_objetos_operativo.RData"
#  )