####### AGREGANDO ESQUEMA DE MUESTREO

################# Cargamos las biblioecas necesarias ####################3
library(questionr)
library(survey)
library(dplyr)
library(stargazer)

#--------------------------------------------------------------------
## Cargamos datos
concentrado_extendido <- read.csv("/home/milo/Documents/egtp/BID/olg_egtp/src/judy/output/enigh2022_ns_concentradohogar_extendido.csv")

concentrado_extendido <- subset(concentrado_extendido, salud >0)

## Generamos gasto catastrófico y empobrecedor
# Calculamos ingreso disponible
concentrado_extendido$ing_disp <- concentrado_extendido$ing_cor - concentrado_extendido$alimentos

# Calculamos gasto catastrófico
concentrado_extendido$hog_gto_catas30 <- ifelse((concentrado_extendido$salud/concentrado_extendido$ing_disp)>=0.30, 1,0)

# Calculamos gasto empobrecedor
concentrado_extendido$hog_empobrecido <- ifelse(concentrado_extendido$ing_cor>=concentrado_extendido$alimentos & (concentrado_extendido$ing_cor-concentrado_extendido$salud) <=concentrado_extendido$alimentos,1,0)

## Definimos el esquema de muestreo
enigh_design<-svydesign(id=~upm, strata=~est_dis, weight=~factor, data=concentrado_extendido, nest=TRUE)
options(survey.lonely.psu="adjust")

#### Generamos tablas 

## Tabla de gasto promedio por claves J 
gasto_tri_clases_j <-svyby(~gasto_tri, ~clave,design = enigh_design , svymean, vartype=c("se","cv"))
gasto_tri_clases_j

## Tabla de gasto promedio por claves J y por seguridad social
gasto_tri_clases_j_ss <-svyby(~gasto_tri, ~clave + seguridad_social,design = enigh_design , svymean, vartype=c("se","cv"))
gasto_tri_clases_j_ss

####  Hacemos regresiones logísticas
reg_1 <- svyglm(hog_gto_catas30 ~ seguridad_social,  design = enigh_design, family=quasibinomial)
reg_2 <- svyglm(hog_gto_catas30 ~ seguridad_social + ing_cor,  design = enigh_design, family=quasibinomial)

stargazer(reg_1, reg_2, type = "text")

