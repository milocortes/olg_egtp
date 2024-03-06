####### AGREGANDO ESQUEMA DE MUESTREO

################# Cargamos las biblioecas necesarias ####################3
library(questionr)
library(survey)
library(ggplot2)

###############################################
### Definimos una ruta temporal para guardar los microdatos
temporal <- tempfile()
download.file("https://www.inegi.org.mx/contenidos/programas/enoe/15ymas/microdatos/2015trim4_csv.zip",temporal)
files = unzip(temporal, list=TRUE)$Name
unzip(temporal, files=files[grepl("csv",files)])
sdemt <- read.csv("ENOEN_SDEMT421.csv")


# Se selecciona la población de referencia que es: población ocupada mayor de 15 años con entrevista completa y condición de residencia válida.
sdemt <- subset(sdemt, sdemt$clase2 == 1 & sdemt$eda>=20 & sdemt$eda<=65 & sdemt$r_def==0 & (sdemt$c_res==1 | sdemt$c_res==3))


##### Formal
sdemt$cohortes_formal<-0

sdemt$cohortes_formal[sdemt$eda >=20 | sdemt$eda < 25 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "20-25-Formal"
sdemt$cohortes_formal[sdemt$eda >=25 | sdemt$eda < 30 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "25-30-Formal"
sdemt$cohortes_formal[sdemt$eda >=30 | sdemt$eda < 35 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "30-35-Formal"
sdemt$cohortes_formal[sdemt$eda >=35 | sdemt$eda < 40 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "35-40-Formal"
sdemt$cohortes_formal[sdemt$eda >=40 | sdemt$eda < 45 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "40-45-Formal"
sdemt$cohortes_formal[sdemt$eda >=45 | sdemt$eda < 50 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "45-50-Formal"
sdemt$cohortes_formal[sdemt$eda >=50 | sdemt$eda < 55 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "50-55-Formal"
sdemt$cohortes_formal[sdemt$eda >=55 | sdemt$eda < 60 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "55-60-Formal"
sdemt$cohortes_formal[sdemt$eda >=60 | sdemt$eda < 65 & (sdemt$mh_col==2 | sdemt$mh_col==4 | sdemt$mh_col==6 | sdemt$mh_col==8 | sdemt$mh_col==10)] <- "60-65-Formal"

#### Informal
sdemt$cohortes_informal<-0

sdemt$cohortes_informal[sdemt$eda >=20 | sdemt$eda < 25 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "20-25-Informal"
sdemt$cohortes_informal[sdemt$eda >=25 | sdemt$eda < 30 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "25-30-Informal"
sdemt$cohortes_informal[sdemt$eda >=30 | sdemt$eda < 35 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "30-35-Informal"
sdemt$cohortes_informal[sdemt$eda >=35 | sdemt$eda < 40 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "35-40-Informal"
sdemt$cohortes_informal[sdemt$eda >=40 | sdemt$eda < 45 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "40-45-Informal"
sdemt$cohortes_informal[sdemt$eda >=45 | sdemt$eda < 50 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "45-50-Informal"
sdemt$cohortes_informal[sdemt$eda >=50 | sdemt$eda < 55 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "50-55-Informal"
sdemt$cohortes_informal[sdemt$eda >=55 | sdemt$eda < 60 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "55-60-Informal"
sdemt$cohortes_informal[sdemt$eda >=60 | sdemt$eda < 65 & (sdemt$mh_col==1 | sdemt$mh_col==3 | sdemt$mh_col==5 | sdemt$mh_col==7 | sdemt$mh_col==9)] <- "60-65-Informal"


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------

##### Formal
sdemt$cohortes_formal<-0

sdemt$cohortes_formal[sdemt$eda >=20 | sdemt$eda < 25 & sdemt$cs_p13_1<5] <- "20-25"
sdemt$cohortes_formal[sdemt$eda >=25 | sdemt$eda < 30 & sdemt$cs_p13_1<5] <- "25-30"
sdemt$cohortes_formal[sdemt$eda >=30 | sdemt$eda < 35 & sdemt$cs_p13_1<5] <- "30-35"
sdemt$cohortes_formal[sdemt$eda >=35 | sdemt$eda < 40 & sdemt$cs_p13_1<5] <- "35-40"
sdemt$cohortes_formal[sdemt$eda >=40 | sdemt$eda < 45 & sdemt$cs_p13_1<5] <- "40-45"
sdemt$cohortes_formal[sdemt$eda >=45 | sdemt$eda < 50 & sdemt$cs_p13_1<5] <- "45-50"
sdemt$cohortes_formal[sdemt$eda >=50 | sdemt$eda < 55 & sdemt$cs_p13_1<5] <- "50-55"
sdemt$cohortes_formal[sdemt$eda >=55 | sdemt$eda < 60 & sdemt$cs_p13_1<5] <- "55-60"
sdemt$cohortes_formal[sdemt$eda >=60 | sdemt$eda < 65 & sdemt$cs_p13_1<5] <- "60-65"

#### Informal
sdemt$cohortes_informal<-0

sdemt$cohortes_informal[sdemt$eda >=20 | sdemt$eda < 25 & sdemt$cs_p13_1>=5] <- "20-25"
sdemt$cohortes_informal[sdemt$eda >=25 | sdemt$eda < 30 & sdemt$cs_p13_1>=5] <- "25-30"
sdemt$cohortes_informal[sdemt$eda >=30 | sdemt$eda < 35 & sdemt$cs_p13_1>=5] <- "30-35"
sdemt$cohortes_informal[sdemt$eda >=35 | sdemt$eda < 40 & sdemt$cs_p13_1>=5] <- "35-40"
sdemt$cohortes_informal[sdemt$eda >=40 | sdemt$eda < 45 & sdemt$cs_p13_1>=5] <- "40-45"
sdemt$cohortes_informal[sdemt$eda >=45 | sdemt$eda < 50 & sdemt$cs_p13_1>=5] <- "45-50"
sdemt$cohortes_informal[sdemt$eda >=50 | sdemt$eda < 55 & sdemt$cs_p13_1>=5] <- "50-55"
sdemt$cohortes_informal[sdemt$eda >=55 | sdemt$eda < 60 & sdemt$cs_p13_1>=5] <- "55-60"
sdemt$cohortes_informal[sdemt$eda >=60 | sdemt$eda < 65 & sdemt$cs_p13_1>=5] <- "60-65"

# Definimos el esquema de muestreo
sdemtdesign<-svydesign(id=~upm, strata=~est_d_tri, weight=~fac_tri, data=sdemt, nest=TRUE)
options(survey.lonely.psu="adjust")

ingocup_cohortes_formal <- svyby(~ingocup, ~cohortes_formal, sdemtdesign, svymean, vartype=c("se","cv"))

ingocup_cohortes_informal <- svyby(~ingocup, ~cohortes_informal, sdemtdesign, svymean, vartype=c("se","cv"))


# Tomamos como base el ingreso mensual del primer cohorte del grupo de trabajo calificado
base <- ingocup_cohortes_formal$ingocup[1]

ingocup_cohortes_formal$ingocup <- ingocup_cohortes_formal$ingocup/base
ingocup_cohortes_informal$ingocup <- ingocup_cohortes_informal$ingocup/base

ingocup_cohortes_formal$tipo <- "formal"
ingocup_cohortes_informal$tipo <- "informal"

names(ingocup_cohortes_formal)[1] <- "cohorte"
names(ingocup_cohortes_informal)[1] <- "cohorte"

ingocup_cohortes <- rbind.data.frame(ingocup_cohortes_formal, ingocup_cohortes_informal)
rownames(ingocup_cohortes) <- NULL
