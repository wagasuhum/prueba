library("exactextractr")
library("sf")
library("terra")
library("dplyr")

# 1. Insumos

# 1.1 modelos de especies
# Primero descromprimir el archivo Nivel1.7z en su respectiva carpeta
espec_mod <- list.files("D:/otros_procesos/biomodelos/Presente_SIMSINAP2023/","tif",full.names=TRUE)  #----- funcionara en debian?

# 1.2 Cargar SINAP alias historico RUNAP a WGS84.

# Función para leer shapefiles desde archivos ZIP
leer_shp_desde_zip <- function(zip_file,epsg = "EPSG:4326") {
  # Intenta leer el shapefile directamente
  tryCatch(
    {
      shp <- read_sf(paste0("/vsizip/", zip_file), quiet = T) %>%  #----- funcionara en debian?
        st_transform(crs = st_crs(epsg)) %>% 
        st_zm()
      # Asegurarse de que las geometrías sean 2D
      return(shp)
    },
    error = function(e) {
      # Si hay un error, listar archivos dentro del ZIP
      archivos_en_zip <- unzip(zip_file, list = TRUE)$Name
      # Buscar archivos .shp
      shp_files <- archivos_en_zip[grep("\\.shp$", archivos_en_zip)]  #----- funcionara en debian?
      # Leer el primer shapefile encontrado en el ZIP
      shp <- read_sf(paste0("/vsizip/", zip_file, "/", shp_files[1]), quiet = T) %>%   #----- funcionara en debian?
        st_transform(crs = st_crs(epsg)) %>% 
        st_zm()
      return(shp)
    }
  )
}

# 10 minutes to load all ap db/918 mb ram
path_ap <- "D:/otros_procesos/ap_sinap/" #----- funcionara en debian?
nms <- list.dirs(path_ap, recursive = F, full.names = F)

shps <- list.dirs(path_ap, recursive = F, full.names = T) %>% 
  lapply(X = ., function(X){
    a <- list.files(X, pattern = "*.zip$", full.names = T, recursive = T )
    return(a)
  }) %>% lapply(X = ., function(X){
    tmpShp <- list()
    for(i in 1:length(X)){
      tmpShp[[i]] <- leer_shp_desde_zip(X[i])
    }
    tmpShp <- bind_rows(tmpShp)
    
    # remover duplicados (algunas AP estan duplicados sus poligonos)
    nm <- sub(".*[^0-9]([0-9]+)$", "\\1", tmpShp$URL)
    tmpShp$id <- nm
    tmpShp <- tmpShp[!duplicated(nm),] 
    
    return(tmpShp)
  })
names(shps) <- nms

# 1.3 Cargar territoriales
territoriales <- read_sf("rep_rri_otros/Territoriales/Territorales.shp") %>% 
  st_transform(st_crs("EPSG:4326"))

# 2. Funciones

# 2.0 Generar Matrices de Presencia-Ausencia
#
# Esta función toma una lista de rutas a archivos TIFF de rasters binarios y un shapefile de polígonos,
# y genera una matriz de presencia-ausencia indicando en qué polígonos está presente cada capa de raster.
#
# rasters: Un objeto de archivo raster (spatRaster) con las capas de especies.
# poly: Un objeto de archivo shapefile con los polígonos de áreas protegidas.
#
# return: Una matriz de presencia-ausencia donde las filas corresponden a las capas de raster
# y las columnas a los polígonos. Los valores son 1 si la especie está presente en el polígono,
# y 0 en caso contrario.

generar_matrices <- function(rasters, poly) {
  
  wCol <- terra::ext(-79.0083333333333, -66.85, -4.23333333333333, 12.4583333333333)
  
  r <- lapply(X = rasters, FUN = function(X){
      terra::crop(terra::rast(X), wCol)
    })
  
  # Convertir la lista de archivos raster en un objeto 'rast' usando el paquete 'terra'
  r <- terra::rast(r)

  # Obtener los nombres de las capas raster
  nm <- names(r)
  
  # Realizar una extracción precisa de los valores de raster dentro de cada polígono
  # del shapefile utilizando el paquete 'exactextractr' y convertir el shapefile a 'sf'
  # La función 'exact_extract' con la opción 'max' devuelve el valor máximo por polígono
  r <- exactextractr::exact_extract(r, sf::st_as_sf(poly), 'max', progress = FALSE)
  
  # Generar la matriz de presencia-ausencia
  # Transponer la matriz para que las filas correspondan a las capas raster (especies)
  # y las columnas a las áreas protegidas (polígonos)
  mat <- t(r)
  
  # Asignar los nombres de las filas a la matriz, que corresponden a los nombres de las capas raster
  rownames(mat) <- nm
  
  # Liberar memoria eliminando objetos no necesarios y ejecutando la recolección de basura
  rm(r); gc()
  
  # Devolver la matriz de presencia-ausencia
  return(mat)
}

# 2.1 Representatividad de especies por periodo de SINAP
#
# raster_spec: raster stack, de modelos de distribución de especies  
# SINAP_shp: shape, de periodo de SINAP
# name_table: vector character, para darle nombre a la matriz de riqueza, describe el contenido y temporalidad del SINAP
#
# return matriz de riqueza de especies en el conjunto de Areas Protegidas (AP) del SINAP
# por especie y vector de representatividad total 

rep_spp_SINAP <- function(raster_spp = espec_mod, SINAP_shp = shps[[56]], name_table = nms[56], gruposSize = 1100){
  
  matTempL <- list()
  
  grupos <- ceiling(length(raster_spp) / gruposSize)
  
  # Función para mostrar la barra de progreso
  progress_bar <- function(current, total) {
    progress <- (current / total) * 100
    cat(sprintf("\rProgress: [%-50s] %d%%", 
                paste(rep("=", progress / 2), collapse = ""), round(progress)))
    flush.console()
  }
  
  for(i in 1:grupos) {
    #i <- 1
    progress_bar(i, grupos)
    
    if(i == 1) {
      r0 <- 1
      r1 <- min(gruposSize, length(raster_spp))  # Asegurarse de no exceder la longitud
    } else if(i == grupos) {
      r0 <- (i - 1) * gruposSize + 1
      r1 <- length(raster_spp)
    } else {
      r0 <- (i - 1) * gruposSize + 1
      r1 <- i * gruposSize
    }
    
    matTempL[[i]] <- generar_matrices(rasters = espec_mod[r0:r1], poly = SINAP_shp)
    
  }
  
  # Completar la barra de progreso al 100%
  progress_bar(grupos, grupos)
  cat("\n")  # Nueva línea después de completar la barra de progreso
  
  temp <- do.call("rbind", matTempL)
  rm(matTempL);gc()

  # Cantidad de Areas Protegidas en las cuales esta presente cada especie
  ap_spec.sum <-apply(temp,1,sum, na.rm = T)
  ap.riqueza.1 <-cbind(temp,ap_spec.sum)
  
  # presencia-ausencia de especies en el conjunto de AP dentro del SINAP
  nc <- ncol(ap.riqueza.1) # columna en la que se ubica el conteo
  pres_espec_all_ap <- replace(ap.riqueza.1[,nc], ap.riqueza.1[,nc] > 0 , 1)
  ap.riqueza.1 <-cbind(ap.riqueza.1,pres_espec_all_ap)
  
  # Porcentaje de la representatividad de especies en las AP, en otras palabras,
  # porcion de las especies usadas en el calculo que estan presentes en el 
  # conjunto del AP's
  sumPresAus <- sum(ap.riqueza.1[,ncol(ap.riqueza.1)], na.rm = T)
  repre_riqueza <-sumPresAus/(length(raster_spp))*100
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

repre_deltaind_AP <- function(shp1 = shps[[2]], rique_ind1 = RP_AP[[2]], tiempo1 = nms[2], shp2 = shps[[3]], rique_ind2 = RP_AP[[3]], tiempo2 = nms[3] ){
  
  nm1 <- shp1$id %>% 
     as.numeric()
  nm2 <- shp2$id %>% 
    as.numeric()
  
  # Areas protegidas que comparten los periodos
  index <-  nm2 %in% nm1
  
  # Restar el valor del aporte de la riqueza por AP en los dos periodos
  diff_21 <-  (rique_ind2[which(index == T) , 1] - rique_ind1[, 1]) %>% abs()
  
  # crear un vector en donde se guarde la diferencia de los dos vectores
  deltaRRi  <- rep(NA, nrow(shp2))
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

repres_ap <- mapply(FUN = rep_spp_SINAP, 
                    SINAP_shp = shps, 
                    name_table = nms, 
                    gruposSize = 1100, 
                    SIMPLIFY = FALSE)
# > tm2-tm1
# Time difference of 5.942718 hours

# extraer representatividad por cada temporalidad
Rep_aps <- lapply(X = repres_ap, FUN = function(X){
    X["Representatividad total"]
  }) %>% 
  bind_rows() %>% 
  as.data.frame()

colnames(Rep_aps) <- "Repre_APs"
rownames(Rep_aps) <- nms

# 3.2 Indicador
# Id: Cambio en el porcentaje de la representatividad de la riqueza de especies
# del SINAP (dRRi)
# Incremento tiempo 1 a tiempo 2 en la representatividad riqueza total

tasa.increm.anual <- cambio_rep_anual(rep_aps = Rep_aps$Repre_APs)

# 3.3 Aporte de cada Area Protegida

# 3.3.1 Aporte de cada AP a la representatividad de la riqueza de especies 
# (RRi). Se evalua el % de representatividad de cada AP al total

RP_AP <- mapply(FUN = repre_ind_AP, 
                    matriz_rique = repres_ap, 
                    tiempo = nms, 
                    SIMPLIFY = FALSE)
# tm2-tm1
# Time difference of 10.91409 secs

# RP_AP1990_AP <- repre_ind_AP(repres_90, "1990")

# 3.3.2 Aporte de cada AP al cambio de la representatividad de la riqueza de 
# especies (dRRi). Se evalua el delta del % de representatividad de cada AP
# total

deltaRP_AP <- mapply(FUN = repre_deltaind_AP, 
                shp1 = shps[1:length(nms)-1], 
                rique_ind1 = RP_AP[1:length(nms)-1],
                tiempo1 = nms[1:length(nms)-1], 
                shp2 = shps[2:length(nms)], 
                rique_ind2 = RP_AP[2:length(nms)],
                tiempo2 = nms[2:length(nms)],
                SIMPLIFY = FALSE)
# tm2-tm1
# Time difference of 3.415342 secs

# 3.3.3 Agregar columnas %RRi y deltaRRi al shapefile original

for (i in 1:length(RP_AP)) {
  #i <- 1
  shps[[i]] <- shps[[i]] %>%
    cbind(round(RP_AP[[i]], 3))
}

for (i in 1:length(deltaRP_AP)) {
  #i <- 1
  a <- i+1
  shps[[a]] <- shps[[a]] %>%
    cbind(round(deltaRP_AP[[i]], 3))
}

# 5 Representatividad riqueza por territorial

# 5.1 Preparación y calculo de representatividad por territorial

# 5.1.1. Intersectar territorial y shapefiles: warnings are ok (equator+small areas)

tm1 <- Sys.time()
shps_territoriales <- lapply(X = shps, FUN = function(X){
  sf::sf_use_s2(FALSE)
  st_intersection(X, territoriales) %>% st_drop_geometry()
})
tm2 <- Sys.time()
# tm2-tm1
# Time difference of 30 min

# 5.1.2. Se agregan los porcentajes por Territorial y se hace una sumatoria.

terr <- list()

for(i in 1:length(shps_territoriales)){
  # i <- 1
  datai <- shps_territoriales[[i]]
  colRRi <- grep(pattern = "^RRi\\d{4}$", x = colnames(datai))
  terri <- aggregate(datai[,colRRi] ~ datai[,"nombre"], data=datai, FUN=sum)
  colnames(terri) <- c("Terr", nms[i])
  terr[[i]] <- terri
  
  rm(datai, terri)
}

# 5.1.3. Se unen todos los valores de cada RUNAP en una sola tabla.

dat_rep_ri_terr <- Reduce(function(x, y) merge(x, y, by = "Terr", all = TRUE), terr)
dat_rep_ri_terr_anual <- as.data.frame(t(dat_rep_ri_terr[-1]))  # Transponer sin la columna "Terr"

# 5.2 Calculo de la diferencia en la representatividad entre territoriales, dRRIterr 

# Calcular la diferencia entre las columnas, manejando NA
tasa.incremento.territorial <- t(apply(dat_rep_ri_terr_anual, 2, function(x) abs(diff(x, na.rm = TRUE))))
row.names(tasa.incremento.territorial) <- dat_rep_ri_terr$Terr

