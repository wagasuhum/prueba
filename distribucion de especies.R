library(SpaDES)
library(raster)
library(xfun)
library(ReIns)
library(dplyr) # para manipular data


# load("rep_distspp_cod/rep_distspp_objetos_operativo.RData")

# 1. Cargar insumos

# 1.1 Cargar rutas de los modelos Nivel1 BioModelos

list_raster <-list.files(path = "C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/BioModelos_N1", pattern="*.tif$", full.names = T)

# 1.2 Leer archivos RUNAP, en formato raster (.tif).

RUNAP_1990 <- raster("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_rasters/RUNAP_1990_pnnid.tif")
RUNAP_2010 <- raster("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/RUNAP_rasters/RUNAP_2010_pnnid.tif")


# 1.3 Cargar archivo shapefiile Territorial

# Cargar shapefile de territoriales con sf
Territoriales <- st_read("C:/Users/walter.garcia/Documents/capas_base_ejemplos/capas_base_ejemplos/Territorial/Territoriales_final.shp")

# 2. Funciones
#
# 2.1 Distribución de especies en el RUNAP (DistSpRUNAP)
#
# path.Spraster: vector character, directorio en donde se encuentra el raster de distribucion de las
# especies
# ras.RUNAP: raster, raster del sistema nacional de areas protegidas, puede ser cualquier 
# capa de poligonos rasterizada
# do.equal_proj: logico, se prende o apaga la reproyección de los BioModelos al mismo CRS de los raster del RUNAP.
# El indicador necesaita que los raster de los BioModelos este en el mismo CRS que los raster del RUNAP.
# Este proceso consume tiempo y recursos, es posible que el proceso se mas rapido dentro de un SIG y posterior
# cargarlo en R. Por defecto TRUE.
#
#
# Return: 
# data.frame con las siguientes columnas
#
# RUNAP: Area Protegida ID
# freq: numero de pixeles de la distribucion de la especie dentro de cada parque
# species: nombre del archivo raster de la distribucion de la especie
# in_RUNAP: ¿existen pixeles dentro del RUNAP?
# per_protected: porcentaje de la distribución de la especie en Colombia que es protegida
# needed_to_protect: porcentaje de la distribución de la especie a proteger
# total_area_km2: area de distribución de la especie en km2
# prop_achieved: proporcion del area de la distribucion a ser protegida que 
# esta siendo realmente protegida
# achieved: 0 si la proporcion de area protegida real (prop_achieved) es < 0.9, y 1 
# si es > 0.9
#
# Details: Rodrigues, A.S.L., et al, (2004a) Effectiveness of the global protected 
# area network in representing species diversity. Nature 428(6983), 640-643.
# Rodrigues, A.S.L., et al. (2004b) Global Gap Analysis: Priority Regions for 
# Expanding the Global Protected-Area Network. BioScience 54(12), 1092-1100.

DistSpRUNAP <- function(path.Spraster, ras.RUNAP, do.equal_proj = T){
  
  print(path.Spraster)
  
  #leer raster por especie
  focus_sp <- raster::raster(path.Spraster)
  
  if(do.equal_proj == T){
    
    # cambiar proyeccion
    template <- projectRaster(from = focus_sp, to= ras.RUNAP)
    
    #template is an empty raster that has the projected extent of r2 but is aligned with r1 (i.e. same resolution, origin, and crs of r1)
    focus_sp <- projectRaster(from = focus_sp, to= template)
    
  }
  
  # Revisar si el raster del BioModelo intercepta con el raster de las Areas Protegidas
  # del RUNAP 
  
  if(!is.null(intersect(extent(ras.RUNAP), focus_sp))){
    
    # Convertir raster a Data.frame
    cells <- rasterToPoints(focus_sp)
    cells <- as.data.frame(cells)
    # Elegir celdas de presencia de la especie
    cells <- subset(cells,cells[3]==1)
    
    # Extraer valores del RUNAP que intersectan con el area de distribución
    # de la especie
    RUNAP_data <- raster::extract(ras.RUNAP, cells[1:2])
    
    # Agregar la intersección del RUNAP (Areas Protegidas) y BioModelo al data.frame
    # cells 
    cells$RUNAP <- RUNAP_data
    
    if(nrow(cells) > 0){
      # contar los pixeles que estan dentro del RUNAP
      InterSpAP <- plyr::count(cells[4])
      InterSpAP$species <- names(focus_sp)
      InterSpAP$year_RUNAP <- strsplit(names(ras.RUNAP),"_")[[1]][[2]]
    }else{
      stop("sin datos en las celdas")
    }
    
    
    # Asignar valores dependiendo si los pixeles estan dentro o fuera del RUNAP, guardar
    # en columna "in_PNN"
    
    InterSpAP[which(is.na(InterSpAP$RUNAP)),"in_RUNAP"] <- 0
    InterSpAP[which(!is.na(InterSpAP$RUNAP)),"in_RUNAP"] <- 1
    
    # Suma de pixeles fuera (0) y dentro (1) del conjunto de Areas Protegidas
    count_RUNAP <- aggregate(InterSpAP$freq, list(RUNAP = InterSpAP$in_RUNAP),"sum")
    
    # Area total de la distribución de la especie
    total_area <- sum(count_RUNAP$x)
    
    # Kilometro cuadrado a la resolución dada del raster y con la proyección 
    # +proj=tmerc +lat_0=4.59620041666667 +lon_0=-74.0775079166667 +k=1 +x_0=1000000 
    # +y_0=1000000 +ellps=GRS80 +units=m +no_defs 
    
    km2ResRas <- ((res(focus_sp)[[1]]/1000)^2)
    
    # Area total a kilometros (metros a kilometros)
    total_area_km2 <- total_area * km2ResRas
    
    # Area total protegida
    total_protected <- count_RUNAP[which(count_RUNAP$RUNAP == "1"), 2] #porque no se paso a km 2?
    
    if(length(total_protected) > 0){
      InterSpAP$per_protected <- (total_protected/total_area) * 100
    }else{
      InterSpAP$per_protected<-0
    }
    
    # Analisis de Brecha
    # Decisiones de area necesaria para proteger
    # area menor a 1000 km2
    if(total_area_km2 <= 1000){
      InterSpAP$needed_to_protect <- 100
    }else{
      # area mayor de 250000 km2
      if(total_area_km2 >= 250000){
        InterSpAP$needed_to_protect <- 10 
      }else{
        # area menor de 250000 km2 pero mayor a 1000 km2
        InterSpAP$needed_to_protect<-lm.func(total_area_km2)
      }
    }
    
    # agregar columna de area de distribución total de la especie
    InterSpAP$total_area_km2 <- total_area_km2
    
    # porcion del area exitosamente protegida de la necesaria
    InterSpAP$prop_achieved <- InterSpAP$per_protected / InterSpAP$needed_to_protect
    
    # se ha logrado proteger mas del 0.9 de la distribución necesaria a proteger?
    if(InterSpAP$prop_achieved[[1]] > 0.9){
      InterSpAP$achieved<-1
    }else{
      InterSpAP$achieved<-0
    }
    
  }else{
    
    # en caso de no intersectar
    # crear una data frame sin datos con todas las columnas
    InterSpAP <- as.data.frame(matrix(nrow = 1,ncol = 10))
    colnames(InterSpAP) <- c("RUNAP", "freq", "species", "year_RUNAP", "in_RUNAP", 
                             "per_protected", "needed_to_protect", "total_area_km2",   
                             "prop_achieved", "achieved")
  }
  
  return(InterSpAP)
  
}

# 2.2 Calculo del area a proteger cuando el area de distribución se encuentra entre
# 1000 km2 y 250000 km2
# 
# x = vector numeric, area de distribucion total de la especie
# 
# return:
# Area minima a protegerpara especies con distribuciones mayores a 1000 km2 y 
# menores a 250000 km2
#
# Details: Rodrigues, A.S.L., et al, (2004a) Effectiveness of the global protected 
# area network in representing species diversity. Nature 428(6983), 640-643.
# Rodrigues, A.S.L., et al. (2004b) Global Gap Analysis: Priority Regions for 
# Expanding the Global Protected-Area Network. BioScience 54(12), 1092-1100.


lm.func<-function(x){
  
  # El area a proteger es calculada como una función del área de distribución 
  # de la especie y se escala como una función lineal > 10% para especies 
  # con distribuciones menores a 250000 km2, hasta < 100% para especies 
  # con distribuciones mayores a 1000 km2
  
  areas <- c(1000, 2.5e5)
  proporcion <- c(100, 10)
  TempDF <- data.frame("xTemp" = areas, "yTemp" = proporcion)
  m <- diff(TempDF$yTemp) / diff(TempDF$xTemp)
  b <- ((TempDF$yTemp[[1]]*TempDF$xTemp[[2]]) - (TempDF$yTemp[[2]]*TempDF$xTemp[[1]])) /
    (TempDF$xTemp[[2]] - TempDF$xTemp[[1]])
  
  y <- m*x + b
  
  return(y)
}

# 2.3 Calcula el valor de la representatividad de la distribución de especies en un sistema de referencia 
# (nacional o territorial) y por area protegida dentro de ese sistema de referencia 
# 
# stats_periodo = estadisticas de la distribución por especie y por area en un sistema de referencia,
# objeto creado con la función DistSpRunap
#
# return: lista con 3 objetos,
#
# media_rep_distSpp = media de la representatividad de la distribución de especes en el sistema de referencia
# media_rep_distSpp_AP = media de la representitividad de la distribucion de especies por Area Protegida 
# id_AP = Identificador unico de cada area protegida

rep_distSpp <- function(stats_periodo = sp_stats_1990, raster_runap = RUNAP_1990){
  
  #stats_allsp es la tabla generada por el loop.
  
  ind_PNN <- stats_periodo %>% count(RUNAP, achieved)
  
  #contar número de especies que están en un área
  PNN_sp <- aggregate(ind_PNN$n,list(RUNAP = ind_PNN$RUNAP), "sum")
  
  #contar cuántas de esas especies alcanzan su target
  PNN_achieved <- subset(ind_PNN, ind_PNN$achieved == 1)
  
  PNN_stats <- merge(PNN_sp, PNN_achieved, "RUNAP", all.x = T)
  
  colnames(PNN_stats)<-c("RUNAP","no_species","ac","sp_ach")
  
  PNN_stats[which(is.na(PNN_stats$ac)),"sp_ach"] <- 0
  
  #calcular indicador por área
  PNN_stats$ac <- (PNN_stats$sp_ach/PNN_stats$no_species)*100
  
  #media indicador:
  
  media_all <-  mean(PNN_stats$ac)
  
  unique_ID_AP <- unique(na.omit(as.data.frame(raster_runap)))
  
  colnames(unique_ID_AP) <- "RUNAP"
  
  no_distsspp <- unique_ID_AP$RUNAP[which(unique_ID_AP$RUNAP %in% PNN_stats$RUNAP == F)]
  
  for(j in 1:nrow(PNN_stats)){
    raster_runap[ raster_runap[] == PNN_stats$RUNAP[j] ] <- PNN_stats$ac[j]  
  }
  
  for(i in 1:length(no_distsspp)){
    raster_runap[ raster_runap[] == no_distsspp[i] ] <- NA 
  }
  
  # objetos para resultados
  id_AP <- PNN_stats$RUNAP
  medias_rep_distrspp_AP <- PNN_stats$ac
  
  
  return(list(media_rep_distSpp = media_all, media_rep_distSpp_AP = medias_rep_distrspp_AP, 
              id_AP = id_AP, raster_rep_distrSpp = raster_runap, AP_stats = PNN_stats))
  
}


# 3. Aplicacion

# 3.1 Nacional

# 3.1.1 calculo de la representatividad de la distribución de especies en SINAP

# 1990
names(RUNAP_1990) <- "RUNAP_1990"

sp_stats_1990 <- lapply(X = list_raster, function(X) DistSpRUNAP(path.Spraster = X, ras.RUNAP = RUNAP_1990)) 
sp_stats_1990 <- do.call(rbind, sp_stats_1990)
rep_distSpp_1990 <- rep_distSpp(stats_periodo = sp_stats_1990, raster_runap = RUNAP_1990)

#2010
names(RUNAP_2010) <- "RUNAP_2010"
sp_stats_2010 <- lapply(X = list_raster, function(X) DistSpRUNAP(path.Spraster = X, ras.RUNAP = RUNAP_2010)) 
sp_stats_2010 <- do.call(rbind, sp_stats_2010)
rep_distSpp_2010 <- rep_distSpp(stats_periodo = sp_stats_2010, raster_runap = RUNAP_2010)

# vector de la media de representatividad de distribucion para cada periodo a nivel nacional
rep_distSpp_Nal <- c(rep_distSpp_1990$media_rep_distSpp, rep_distSpp_2010$media_rep_distSpp) %>% 
  round(3)


# 3,1,2 calculo de la diferencia de integridad entre años (indicador)
Delta_rep_distSpp_Nal <- c(0, diff(rep_distSpp_Nal)) %>% round(3)

# 3.2 Territorial

# 3.2.1 calculo de la representatividad de la distribución de especies en SIRAP

# 1990
rep_distSpp_Terr_1990 <- lapply(1:nrow(Territoriales), function(x) { 
  
  # cortar las areas protegidas a la extension de cada territorial
  raster.Areas.terrx <-  raster::crop(rep_distSpp_1990$raster_rep_distrSpp, Territoriales[x, ]) %>%
    raster::mask(Territoriales[x, ]) 
  
  
  df.Areas.terrx <- as.data.frame(raster.Areas.terrx) %>% na.omit() %>% unique()
  
  media_territorialx <- mean(df.Areas.terrx[ , 1], na.rm = T)
  
  return(list(raster.terrx = raster.Areas.terrx, media_terrx = media_territorialx))
  
}
)

# 2010
rep_distSpp_Terr_2010 <- lapply(1:nrow(Territoriales), function(x) { 
  
  # cortar las areas protegidas a la extension de cada territorial
  raster.Areas.terrx <-  raster::crop(rep_distSpp_2010$raster_rep_distrSpp, Territoriales[x, ]) %>%
    raster::mask(Territoriales[x, ]) 
  
  df.Areas.terrx <- as.data.frame(raster.Areas.terrx) %>% na.omit() %>% unique()
  
  media_territorialx <- mean(df.Areas.terrx[ , 1], na.rm = T)
  
  return(list(raster.terrx = raster.Areas.terrx, media_terrx = media_territorialx))
  
}
)

# 3.2.2 Calcular delta de la la representatividad de la distribución de especies en SIRAP

# Cada una de las listas de calculo la representatividad de la distribución de especies por territorial 
# tiene el mismo orden que el shapefile de las territoriales, que se puede obtener con Territoriales$nombre

rep_distSpp_Terr_res <- list()

for(i in 1:length(Territoriales)){
  
  rep_distSpp_terri <- c(rep_distSpp_Terr_1990[[i]]$media_terrx, 
                         rep_distSpp_Terr_2010[[i]]$media_terrx)
  
  Delta_rep_distSpp_terri <- c(0, diff(rep_distSpp_terri))
  
  rep_df <- data.frame(c(rep_distSpp_terri, Delta_rep_distSpp_terri))
  
  rep_distSpp_Terr_res[[i]] <- rep_df
}

rep_distSpp_Terr_res <- do.call(cbind.data.frame, rep_distSpp_Terr_res) %>% t() %>% round(3)

tiempos <- c("1990", "2010")
colnames(rep_distSpp_Terr_res) <- c(paste0("rep_distSpp_", tiempos), paste0("Delta_rep_distSpp_", tiempos))
rownames(rep_distSpp_Terr_res) <- Territoriales$nombre

# 3.3 Area protegida

# 3.3.1 calculo de la representatividad de la distribución de especies por AP

# se desarrollo en el 3.1.1

# 3.3.2 Calcular delta de la la representatividad de la distribución de especies por AP

drep_distSpp_AP_1990_2010 <- rep_distSpp_2010$raster_rep_distrSpp - rep_distSpp_1990$raster_rep_distrSpp
names(drep_distSpp_AP_1990_2010) <- "drep_distSpp_AP_1990_2010"
